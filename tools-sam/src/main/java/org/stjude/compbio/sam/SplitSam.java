package org.stjude.compbio.sam;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SAMTag;

import org.apache.commons.cli.Options;
import org.stjude.compbio.common.RefNameMap;
import org.stjude.compbio.util.CountingMap;
import org.stjude.compbio.util.DelimitedListBuilder;
import org.stjude.compbio.util.MainHelper;

/**
 * Splits a sam/bam file into multiple smaller pieces
 */
public class SplitSam {
	private static final String READ_CTR = "records read";
	private static final String WRITE_CTR = "records written";

	/**
	 * Default reference sequence name (not included in output filename)
	 */
	private static final String DEFAULT_REF = "default";
	/**
	 * Reference sequence name for unmapped reads with no reference
	 */
	private static final String NOREF_REF = "noref";
	
	/**
	 * Pattern for stripping mate suffixes
	 */
	private static final Pattern MATE_SUFFIX_PATTERN =
			Pattern.compile("/[12]$");
	
	private boolean filterNonPrimary;
	private boolean byRefEnabled;
	private boolean byReadEnabled;
	private boolean byReadGroupEnabled;
	private boolean addReadGroupIdToRecordsEnabled;
	private String readGroupIdToAddToRecords = null;
	private boolean mated;
	private boolean stripMateSuffixEnabled;
	private int readsPerFile;
	private int suffixLen;
	private int numBuckets;
	private Integer lenDivisor;
	private Map<String, Integer> refNumBuckets;
	private File outDir;
	private String outPrefix;
	private String outExt;
	private SamReader reader;
	private SAMFileHeader header;
	private SAMFileWriterFactory factory;
	private CountingMap<String> counters;
	
	public SplitSam() {
		this.factory = new SAMFileWriterFactory();
		this.counters = CountingMap.newInitialized(0, READ_CTR, WRITE_CTR);
	}

	/**
	 * Enable/disable the non-primary filter
	 */
	public void setNonPrimaryFilterEnabled(boolean filter) {
		this.filterNonPrimary = filter;
	}
	
	/**
	 * Enabled mated processing with the given number of buckets
	 * 
	 * @param numBuckets number of buckets
	 * @param lenDivisor length divisor (optional)
	 */
	public void enableMated(int numBuckets, Integer lenDivisor) {
		this.mated = true;
		this.numBuckets = numBuckets;
		this.refNumBuckets = new RefNameMap<Integer>();
		this.lenDivisor = lenDivisor;
	}
	
	/**
	 * Enables by-read processing with the given number of reads per file
	 * 
	 * @param readsPerFile number of reads per file
	 */
	public void enableByRead(int readsPerFile) {
		this.byReadEnabled = true;
		this.readsPerFile = readsPerFile;
	}

	/**
	 * @param byRefEnabled the byRefEnabled to set
	 */
	public void setByRefEnabled(boolean byRefEnabled) {
		this.byRefEnabled = byRefEnabled;
	}

	/**
	 * @param byReadGroupEnabled the byReadGroupEnabled to set
	 */
	public void setByReadGroupEnabled(boolean byReadGroupEnabled) {
		this.byReadGroupEnabled = byReadGroupEnabled;
	}

	/**
	 * @param addReadGroupIdToRecordsEnabled the addReadGroupIdToRecordsEnabled to set
	 */
	public void setAddReadGroupIdToRecordsEnabled(
			boolean addReadGroupIdToRecordsEnabled) {
		this.addReadGroupIdToRecordsEnabled = addReadGroupIdToRecordsEnabled;
	}

	/**
	 * @param stripMateSuffixEnabled the stripMateSuffixEnabled to set
	 */
	public void setStripMateSuffixEnabled(boolean stripMateSuffixEnabled) {
		this.stripMateSuffixEnabled = stripMateSuffixEnabled;
	}

	/**
	 * @param suffixLen the suffixLen to set
	 */
	public void setSuffixLen(int suffixLen) {
		this.suffixLen = suffixLen;
	}

	/**
	 * Sets the output directory, base prefix, and extension
	 * @param dir output directory
	 * @param prefix basename prefix
	 * @param ext extension, without the . (should be sam or bam)
	 */
	public void setOutput(File dir, String prefix, String ext) {
		this.outDir = dir;
		this.outPrefix = prefix;
		this.outExt = ext;
	}

	/**
	 * @param header the header to set
	 */
	public void setHeader(SAMFileHeader header) {
		this.header = header;
		
		// Set reference-specific numBuckets values if applicable
		if(lenDivisor != null) {
			for(SAMSequenceRecord seqRec:
					header.getSequenceDictionary().getSequences()) {
				refNumBuckets.put(seqRec.getSequenceName(),
						seqRec.getSequenceLength() / lenDivisor + 1);
			}
		}
		
		// Remember single RGID if there is one
		if(addReadGroupIdToRecordsEnabled) {
			List<SAMReadGroupRecord> readGroups = header.getReadGroups();
			if(readGroups != null && readGroups.size() == 1) {
				readGroupIdToAddToRecords = readGroups.get(0).getId();
			}
			else {
				System.err.println(
						"Cannot assign read group IDs because there are " +
						readGroups.size() + " in the header (requires " +
						"exactly 1 read group)");
				readGroupIdToAddToRecords = null;
			}
		}
	}

	/**
	 * @param reader the reader to set
	 */
	public void setReader(SamReader reader) {
		this.reader = reader;
		this.setHeader(reader.getFileHeader());
	}

	/**
	 * @return the counters
	 */
	public CountingMap<String> getCounters() {
		return counters;
	}

	/**
	 * Helper for the main method
	 */
	private static SamMainHelper main = new SamMainHelper(getMainOptions());

	private static void usage(int exitcode) {
		main.usage(
				"java " + SplitSam.class.getName() + " [OPTION]... PREFIX EXT\n" +
				"Splits a sam/bam file into smaller pieces.  Output is\n" +
				"named with prefix PREFIX and extension EXT, with\n" +
				"distinguishing parts in between.  The extension\n" +
				"determines the output format.  PREFIX may contain path.\n\n" +
				"You may split by reference name, read group, and/or by\n" +
				"number of reads.\n\n" +
				"If your data is paired-end, and you would like to preserve\n" +
				"sort order while keeping mates in the same file, then you\n" +
				"can use -m with -b (and optionally -l).  In this mode, all\n" +
				"inputs are divided by read name into B buckets.  You can\n" +
				"specify the number of buckets using -b, or compute by taking" +
				"chromosome length divided by the -l value.  If yo use -l,\n" +
				"you still need -b for the no-ref case.\n\n");
		System.exit(exitcode);
	}

	public static Options getMainOptions() {
		Options options = new Options();
		options.addOption("i", true, "input sam/bam file if not stdin");
		options.addOption("c", false, "split by ref name (usually chromosome)");
		options.addOption("n", true, "split every 'arg' records");
		options.addOption("a", true, "suffix length for -n default 3");
		options.addOption("r", false, "split by read group");
		MainHelper.addLongOpt(options, "add-unique-rgid-to-records", null, 
				"add RGID to records where missing, if and only if there is " +
				"exactly one RG in header");
		options.addOption("m", false, "keep mates together");
		options.addOption("M", false, "strips mate suffixes /1 and /2");
		options.addOption("b", true, "number of buckets for -m splitting");
		options.addOption("l", true, "seq len divisor to compute b on the fly");
		options.addOption("p", false, "filter out non-primary records");
		return options;
	}

	/**
	 * @param args
	 * @throws FileNotFoundException
	 */
	public static void main(String[] args) throws IOException {
		main.addValidationStringencyOption();
		if(args.length == 0) usage(0);
		
		// Create worker
		SplitSam worker = new SplitSam();
		
		// Parse command line directly into worker
		main.setArgs(args);
		File outPrefixFile = main.getOutputFile(0);
		String outExt = main.getArgString(1);
		main.setStrict(false);
		worker.setNonPrimaryFilterEnabled(main.hasArg("p"));
		boolean byRef = main.hasArg("c");
		boolean byRead = main.hasArg("n");
		boolean mated = main.hasArg("m");
		boolean stripMateSuffix = main.hasArg("M");
		boolean byRg = main.hasArg("r");
		if(mated) {
			if(byRead) {
				throw new MainHelper.CommandLineException("Cannot do both -n and -m");
			}
			else if(!main.hasArg("b")) {
				throw new MainHelper.CommandLineException("-m requires -b");
			}
			else if(main.hasArg("l") && !byRef) {
				System.err.println("WARNING: -l is ignored without -c");
			}
			int numBuckets = main.getArgInt("b", 0);
			worker.enableMated(numBuckets, main.getArgInt("l"));
		}
		else if(byRead) {
			worker.enableByRead(main.getArgInt("n"));
		}
		worker.setAddReadGroupIdToRecordsEnabled(
				main.hasArg("add-unique-rgid-to-records"));
		worker.setByRefEnabled(byRef);
		worker.setByReadGroupEnabled(byRg);
		worker.setStripMateSuffixEnabled(stripMateSuffix);
		worker.setSuffixLen(main.getArgInt("a", 3));
		
		// Get and validate output file prefix and extension
		worker.setOutput(
				outPrefixFile.getParentFile(),
				outPrefixFile.getName(),
				outExt);
		
		// Open input
		SamReader reader = main.openSamReader("i");
		worker.setReader(reader);
		
		// Read the header
		SAMFileHeader header = reader.getFileHeader();
		worker.setHeader(header);
		
		// Do the split
		try {
			worker.split();
		}
		finally {
			reader.close();
			System.err.println("Stats:\n" + worker.getCounters());
		}
	}
	
	public void split() {		
		// Initialize RefWriter map
		Map<String, RefWriter> writers = new HashMap<String, RefWriter>();
		
		// Do the split
		try {
			for(SAMRecord record: reader) {
				// Count it
				counters.increment(READ_CTR);

				// Check filter
				if(filterNonPrimary && record.getNotPrimaryAlignmentFlag()) {
					continue;
				}
				
				// Strip off /1 /2 suffixes if appropriate
				// Note: this has to happen before getting the output key, or
				// there can be problems with -m
				if(stripMateSuffixEnabled) {
					Matcher m = MATE_SUFFIX_PATTERN.matcher(
							record.getReadName());
					record.setReadName(m.replaceFirst(""));
				}
				
				// Add RGID if appropriate
				if(readGroupIdToAddToRecords != null &&
						record.getReadGroup() == null) {
					record.setAttribute(
							SAMTag.RG.name(),
							readGroupIdToAddToRecords);
				}
				
				// Write it to the appropriate writer
				
				// Get string key identifying the output file to which to write
				String outputKey = getOutputKey(record);
				
				// Get the RefWriter to use
				if(!writers.containsKey(outputKey)) {
					writers.put(outputKey, new RefWriter(outputKey));
				}
				RefWriter writer = writers.get(outputKey);
				
				// If it is not open, then open it
				if(!writer.isOpen()) writer.open();
				
				// Write the record
				writer.addAlignment(record);
				counters.increment(WRITE_CTR);
			}
		}
		finally {
			for(RefWriter writer: writers.values()) writer.close();
		}
	}

	/**
	 * Gets a string output key identifying the writer to which to write
	 * 
	 * @param record record to test
	 * @return output key
	 */
	private String getOutputKey(SAMRecord record) {
		// Get reference name, casting * to null
		String refName = record.getReferenceName();
		if("*".equals(refName)) refName = null;

		// Build up a list of key parts
		List<String> parts = new ArrayList<String>();
		if(byReadGroupEnabled) {
			SAMReadGroupRecord readGroup = record.getReadGroup();
			if(readGroup != null) parts.add(readGroup.getId());
		}
		if(byRefEnabled) {
			parts.add((refName == null) ? NOREF_REF : refName);
		}
		if(mated) {
			// Determine number of buckets
			int numBuckets;
			if(refName == null || lenDivisor == null) {
				numBuckets = this.numBuckets;
			}
			else {
				numBuckets = refNumBuckets.get(refName);
			}
			// Determine and add bucket
			int bucket = Math.abs(record.getReadName().hashCode()) % numBuckets;
			parts.add(String.format("%0" + suffixLen + "d", bucket));
		}
		
		// If we got no parts, then return default ref
		if(parts.isEmpty()) return DEFAULT_REF;

		// Otherwise, put parts together
		DelimitedListBuilder dlb = new DelimitedListBuilder("-");
		dlb.appendAll(parts);
		return dlb.toString();
	}

	private class RefWriter {
		/**
		 * String key identifying the output file
		 */
		String outputKey;
		/**
		 * The open SAMFileWriter (may be null)
		 */
		SAMFileWriter writer;
		/**
		 * The next file number to use
		 */
		int next;
		/**
		 * The number written on this writer
		 */
		int written;
		
		public RefWriter(String outputKey) {
			this.outputKey = outputKey;
			writer = null;
			next = 0;
			written = 0;
		}
		
		/**
		 * Returns true if writer is ready to accept alignments
		 */
		public boolean isOpen() {
			return writer != null;
		}
		
		/**
		 * Adds an alignment using open writer
		 * 
		 * @param record record to write
		 * @throws NullPointerException if writer is not open
		 */
		public void addAlignment(SAMRecord record) {
			writer.addAlignment(record);
			++written;
			if(byReadEnabled && written >= readsPerFile) {
				writer.close();
				writer = null;
			}
		}
		
		/**
		 * Opens the writer, creating a new file if necessary
		 */
		public void open() {
			// Build up filename
			DelimitedListBuilder suffix = new DelimitedListBuilder("-");
			String fileName = outPrefix;
			if(!outputKey.equals("") && outputKey != DEFAULT_REF) {
				suffix.append(outputKey);
			}
			if(byReadEnabled) {
				suffix.append(String.format("%0" + suffixLen + "d", next++));
			}
			fileName = fileName + suffix.toString() + "." + outExt;
			File file = new File(outDir, fileName);

			// Create writer
			writer = factory.makeSAMOrBAMWriter(header, true, file);
			written = 0;
		}
		
		/**
		 * Closes the writer (safe to use if not open)
		 */
		public void close() {
			if(writer != null) writer.close();
		}
	}
}
