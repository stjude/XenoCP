package org.stjude.compbio.sam;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;

import org.apache.commons.cli.Options;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.sam.block.Block;
import org.stjude.compbio.sam.block.Blocklist;
import org.stjude.compbio.sam.block.BlocklistFactory;
import org.stjude.compbio.util.CountingMap;
import org.stjude.compbio.util.MainHelper;
import org.stjude.compbio.util.MainHelper.CommandLineException;
import org.stjude.compbio.util.RegexMap;
import org.stjude.compbio.util.io.IoRecord;
import org.stjude.compbio.util.io.TabularTextReader;
import org.stjude.compbio.util.io.TabularTextReader.Header;
import org.stjude.compbio.common.formats.seqdict.SeqDictRefNameMapFactory;
import org.stjude.compbio.common.RefNameMap;

/**
 * Copies a sam/bam file, performing any of a set of simple tweaks, such as
 * sam/bam format conversion, setting sorting, fixing sorting if the file is
 * almost sorted, etc.
 * 
 * @author mrusch
 */
public class TweakSam implements Closeable {
	private static final String READ_CTR = "records read";
	private static final String WRITE_CTR = "records written";
	private static final String DUP_CTR =
			"records filtered due to read id duplication";
	private static final String FAILED_RG_CTR =
			"records had no match for RG assignment";
	private static final BlocklistFactory BLOCKLIST_FACTORY =
			new BlocklistFactory(false, null);
	
	/**
	 * Input reader
	 */
	private SamReader reader;
	/**
	 * Input iterator (opened on init)
	 */
	private SAMRecordIterator iterator = null;
	
	/**
	 * Output writer
	 */
	private SAMFileWriter writer = null;
	
	/**
	 * Read group assignments, null if disabled
	 */
	private RegexMap<String> rgAssignments = null;
	/**
	 * ReadUnmapper, null if disabled
	 */
	private ReadUnmapper readUnmapper = null;
	/*
	 * Locus options
	 */
	private Integer pos;
	private int left;
	private int right;
	/*
	 * Various options
	 */
	private char allele;
	private boolean crossRef;
	private boolean hasCrossRefFilter;
	private ReadLenChooser rlc;
	private RandomReadChooser rrc;
	private double minMatchRate;
	private int inclFlags;
	private int exclFlags;
	private boolean clearDupEnabled;
	private boolean filterConsecDupsEnabled;
	private boolean filterAllDupsEnabled;
	private boolean normalizePairsEnabled;
	private boolean makeUbamEnabled;
	private long maxRecords = Long.MAX_VALUE;
	
	/**
	 * Read counters
	 */
	private final CountingMap<String> counters;
	
	/**
	 * Instantiates a TweakSam worker from a SamReader
	 * 
	 * @param reader open reader
	 */
	public TweakSam(SamReader reader) {
		this.reader = reader;
		this.counters = CountingMap.newInitialized(0,
				READ_CTR, WRITE_CTR, FAILED_RG_CTR);
	}
	
	/**
	 * Gets the read counters
	 */
	public CountingMap<String> getCounters() {
		return counters;
	}

	/**
	 * Opens the output writer
	 * 
	 * @param outputFile output file, or null for System.out
	 * @param headerFile header file, or null to use input header
	 * @param sortOrder sort order for sorting or declaration, null for input
	 * @param setSortOrderEnabled if true, then sort order is a declaration
	 * @param md5 whether to generate md5 file
	 */
	public void openOutputWriter(
			File outputFile,
			File headerFile,
			SortOrder sortOrder,
			boolean setSortOrderEnabled,
			boolean md5) {
		// Get header
		SAMFileHeader header;
		if(headerFile != null) {
			header = SamHeaderHelper.ingestHeader(headerFile).clone();
			System.err.println("Read header with " +
					header.getSequenceDictionary().size() + "SQs, " +
					header.getReadGroups().size() + "RGs, and " +
					header.getProgramRecords().size() + "PGs");
		}
		else {
			header = reader.getFileHeader();
		}
		
		// If a sort order was set, then use it
		if(sortOrder != null) header.setSortOrder(sortOrder);

		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		factory.setCreateMd5File(md5);
		boolean presorted = setSortOrderEnabled || sortOrder == null;
		if(outputFile != null) {
			writer = factory.makeSAMOrBAMWriter(header, presorted, outputFile);
		}
		else {
			writer = factory.makeSAMWriter(header, presorted, System.out);
		}
		
		// Wrap it in a pair normalize wrapper if appropriate
		if(normalizePairsEnabled) {
			writer = new PairNormalizingSAMFileWriter(writer);
		}
	}
	
	/**
	 * Initializes to read the entire input file
	 */
	public void initWholeInput() {
		if(iterator != null) {
			throw new IllegalStateException("Already initialized.");
		}
		iterator = reader.iterator();
	}
		
	/**
	 * Initializes to read all reads that have no reference name.
	 * 
	 * Requires indexed input file
	 * 
	 * @throws UnsupportedOperationException if reader is not indexed
	 */
	public void initNoRefName() throws UnsupportedOperationException {
		if(iterator != null) {
			throw new IllegalStateException("Already initialized.");
		}
		if(!reader.hasIndex()) {
			throw new UnsupportedOperationException(
					"Must use indexed BAM when processing reference-less " +
					"records");
		}
		iterator = reader.queryUnmapped();
	}
	
	/**
	 * Initializes to read a particular reference or reference/position
	 * 
	 * @param refName reference name (required)
	 * @param pos position (or null/0 for whole reference)
	 * @param left left margin (ignored if pos is null)
	 * @param right right margin (ignored if pos is null)
	 */
	public void initLocus(String refName, Integer pos, int left, int right) {
		if(iterator != null) {
			throw new IllegalStateException("Already initialized.");
		}
		if(!reader.hasIndex()) {
			throw new UnsupportedOperationException(
					"Must use indexed BAM when processing reference-less " +
					"records");
		}
		RefNameMap<String> refNameMap = SeqDictRefNameMapFactory.newInstance(reader.getFileHeader().getSequenceDictionary());
		refName = refNameMap.get(refName);
		if(pos == null || pos == 0) {
			iterator = reader.queryContained(refName, 0, 0);
		}
		else {
			this.pos = pos;
			this.left = left;
			this.right = right;
			iterator = reader.queryOverlapping(
					refName, pos - left, pos + right);
		}

	}
	
	/**
	 * Enables assigning read groups based on definitions in the file
	 * 
	 * @param file tab-delimited text file containing read name regex and RG
	 * @throws IOException
	 */
	public void enableRgAssignment(File file) throws IOException {
		rgAssignments = new RegexMap<String>();
		TabularTextReader<?> tr =
				TabularTextReader.open(file, "\t", Header.forbidden());
		for(IoRecord<?> record: tr) {
			rgAssignments.put(
					Pattern.compile(record.valueForOrdinal(0)),
					record.valueForOrdinal(1));
		}
	}
	
	/**
	 * Enables unmapping using a names file
	 * 
	 * @param file read names file
	 */
	public void enableUnmapping(File file) throws IOException {
		readUnmapper = new ReadUnmapper();
		readUnmapper.readNamesFile(file);
	}

	/**
	 * Requires a particular allele at the position
	 */
	public void setAllele(char allele) {
		this.allele = allele;
	}
	
	/**
	 * Enables the cross-ref filter
	 * 
	 * @param crossRef true to require cross-ref, false to forbid
	 */
	public void enableCrossRefFilter(boolean crossRef) {
		this.hasCrossRefFilter = true;
		this.crossRef = crossRef;
	}

	/**
	 * Sets a record chooser that uses read length
	 */
	public void setReadLenChooser(ReadLenChooser rlc) {
		this.rlc = rlc;
	}
	
	public void setRandomReadChooser(RandomReadChooser rrc) {
		this.rrc = rrc;
	}
	
	public void setMinMatchRate(double minMatchRate) {
		this.minMatchRate = minMatchRate;
		if(minMatchRate != 0.0) BLOCKLIST_FACTORY.setMdParsingEnabled(true);
	}
	public void setInclFlags(int inclFlags) {
		this.inclFlags = inclFlags;
	}
	public void setExclFlags(int exclFlags) {
		this.exclFlags = exclFlags;
	}
	public void setClearDupEnabled(boolean clearDupEnabled) {
		this.clearDupEnabled = clearDupEnabled;
	}
	public void setFilterConsecDupsEnabled(boolean filterConsecDupsEnabled) {
		this.filterConsecDupsEnabled = filterConsecDupsEnabled;
	}
	public void setFilterAllDupsEnabled(boolean filterAllDupsEnabled) {
		this.filterAllDupsEnabled = filterAllDupsEnabled;
	}
	public void setNormalizePairsEnabled(boolean normalizePairsEnabled) {
		this.normalizePairsEnabled = normalizePairsEnabled;
	}
	public void setMakeUbamEnabled(boolean makeUbam) {
		this.makeUbamEnabled = makeUbam;
	}
	public void setMaxRecords(long maxRecords) {
		this.maxRecords = maxRecords;
	}

	/**
	 * Helper for the main method
	 */
	private static SamMainHelper main = new SamMainHelper(getMainOptions());

	private static void usage(int exitcode) {
		main.usage(
				"java -jar TweakSam.jar [OPTION]...\n" +
				"Copies a sam/bam file while performing any/all of the\n" +
				"following tweaks among others:\n" +
				"1. Format conversion.  Output format is determined by\n" +
				"   extension.  Stdout output is always SAM.\n" +
				"2. Sorting for almost sorted data.  If the data is almost\n" +
				"   sorted, then sorting may be fixed using the -f option.\n" +
				"   A buffer size may be specified with -b.  If any record\n" +
				"   is further away in the input than (buffer size), then\n" +
				"   the command will fail.  The sort order to use is\n" +
				"   given by -O.\n" +
				"3. Sort order declaration.  If the data is already sorted\n" +
				"   but the order is not declared, or if you are fixing\n" +
				"   sorting using -f, you may use -O to declare the\n" +
				"   correct sort order.  This will be set in the output\n" +
				"   header.\n" +
				"4. Extraction by sequence name.  Specify a reference\n" +
				"   sequence name using -c or specify -n to extract records\n" +
				"   with no reference sequence name.  These are only valid\n" +
				"   for coordinate-sorted, indexed bam files" +
				"5. CIGAR overflow detection.  Detects overflow of the\n" +
				"   numeric part of a CIGAR element and soft clips the\n" +
				"   alignment at that position.\n" +
				"6. Random sampling.  Specify a probability using -r to\n" +
				"   randomly sample from the input bam.  Mate pairs should\n" +
				"   NOT be broken up, but sorting is preserved.\n" +
				"7. Alignment and allele selection.  Specify a position\n" +
				"   using -c and -p to extract reads aligned at that\n" +
				"   position.  Use -a to specify the allele at that\n" +
				"   position, or -L and/or -R to specify an interval over\n" +
				"   which there must be contiguous alignment.  Use -M to\n" +
				"   require a certain match rate on either side of the -p\n" +
				"   position.\n" +
				"8. Flag filtering using -g and -G like samtools view -f\n" +
				"   and -F\n" +
				"9. Others (see arg list)");
		System.exit(exitcode);
	}

	public static Options getMainOptions() {
		Options options = new Options();
		options.addOption("i", true, "input sam/bam file if not stdin");
		options.addOption("o", true, "output sam/bam file if not stdout");
		options.addOption("5", false, "create MD5 file while writing output");
		MainHelper.addLongOpt(options, "set-sort-order", "ORDER",
				"output sort order to declare: queryname, coordinate, or " +
				"unsorted.  Input must already be sorted in the proper " +
				"order; if it is not, then use -O");
		options.addOption("c", true, "only include records with the given\n" +
				"reference name (requires indexed bam input)");
		options.addOption("n", false, "only include records with no\n" +
				"reference name (requires indexed bam input)");
		options.addOption("x", false,
				"only include paired records whose mates are mapped to a\n" +
				"different chromosome");
		options.addOption("X", false,
				"exclude paired records whose mates are mapped to a\n" +
				"different chromosome");
		options.addOption("l", true,
				"only include records with read length in range; formats\n" +
				"are N >N <N [M,N] (M,N) [M,N) (M,N]");
		options.addOption("r", true,
				"randomly sample records with the given probability of\n" +
				"including any one record.  Range is (0.0,1.0).\n");
		options.addOption("p", true,
				"show only records that have a base aligned to the given\n" +
				"position (requires -c)");
		options.addOption("a", true,
				"show only records with the given allele at the position\n" +
				"specified with -p (valid values: C, G, T, A; requires -p");
		options.addOption("L", true,
				"requires -p; indicates that the read must also align for\n" +
				"n bases to the left of p, for a total of n+1 aligned bases");
		options.addOption("R", true,
				"requires -p; indicates that the read must also align for\n" +
				"n bases to the right of p, for a total of n+1 aligned bases");
		options.addOption("m", true,
				"maximum number of reads to return");
		options.addOption("M", true,
				"minimum match rate in the alignment block on either side\n" +
				"of the position indicated by -p.  The block is split at\n" +
				"the -p position, and each part of the block is tested\n" +
				"separately; both parts must pass for the read to be\n" +
				"included.");
		options.addOption("g", true,
				"show only records with the specified flag(s) set");
		options.addOption("G", true,
				"show only records with the specified flag(s) unset");
		options.addOption("D", false, "exclude marked duplicates");
		options.addOption("d", false, "clear duplicate flag");
		options.addOption("C", false, "exclude consecutive duplicates in " +
				"queryname-sorted input (these are not marked dups, they " +
				"are dups in read name/pairing flags)");
		MainHelper.addLongOpt(options, "filter-read-id-dups", null,
				"filter out records that are duplicates based on read id " +
				"(read name + pairing flags).  When there are multiple " +
				"records for the same read, only the first is written. " +
				"Note: significantly increases memory requirements.");
		MainHelper.addLongOpt(options, "normalize-pairs", null,
				"Aggressively works to normalize read pairs, regardless " +
				"of input sorting.  May dramatically increase memory " +
				"requirements.");
		options.addOption("h", true,
				"replace header with one from the specified file");
		options.addOption("y", true,
				"assign read groups according to mapping in file; the file " +
				"is tab-delimited text with a regex and an arbitrary " +
				"identifier.  If the read name matches the regex, then the " +
				"read group ID is set to the ID indicated in the file");
		options.addOption("q", false,
				"filter out records marked as 'duplicate' or 'fails vendor " +
				"quality checks'");
		options.addOption("u", true,
				"set reads named in the file to unmapped; add XU tag");
		options.addOption("U", false,
				"strip all alignment info; make it a UBAM");
		return options;
	}

	public static void main(String[] args) throws CommandLineException, IllegalStateException, IOException {
		main.addValidationStringencyOption();
		main.addSortOrderOption();
		if(args.length == 0) usage(0);
		
		// Parse command line
		main.setArgs(args);
		main.setStrict(false);
		String refName = main.getArgString("c");
		boolean noRefName = main.hasArg("n");
		boolean hasRefNameFilter = noRefName || refName != null;
		boolean crossRef = main.hasArg("x");
		boolean hasCrossRefFilter = crossRef || main.hasArg("X");
		ReadLenChooser rlc = null;
		if(main.hasArg("l")) {
			rlc = ReadLenChooser.parse(main.getArgString("l"));
		}
		RandomReadChooser rrc = null;
		if(main.hasArg("r")) {
			double prob = Double.parseDouble(main.getArgString("r"));
			if(prob <= 0 || prob >= 1) {
				throw new MainHelper.CommandLineException(
						"Invalid probability: " + prob +
						", must be in the range (0.0,1.0)");
			}
			rrc = new RandomReadChooser(prob);
		}
		int pos = main.getArgInt("p", 0);
		char allele = 0;
		if(main.hasArg("a")) {
			String alleleStr = main.getArgString("a");
			if(alleleStr.length() != 1) {
				throw new MainHelper.CommandLineException(
						"Invalid allele: " + alleleStr);
			}
			allele = alleleStr.toUpperCase().charAt(0);
			if(allele != 'C' && allele != 'G' && 
					allele != 'T' && allele != 'A') {
				throw new MainHelper.CommandLineException(
						"Invalid allele: " + alleleStr);
			}
		}
		int left = main.getArgInt("L", 0);
		int right = main.getArgInt("R", 0);
		double minMatchRate = 0.0;
		if(main.hasArg("M")) {
			minMatchRate = Double.parseDouble(main.getArgString("M"));
		}
		int inclFlags = main.getArgInt("g", 0);
		int exclFlags = main.getArgInt("G", 0);
		if(main.hasArg("D")) {
			exclFlags |= 0x400; // SAMRecord.DUPLICATE_READ_FLAG;
		}
		if(main.hasArg("q")) {
			exclFlags |= 0x400; // SAMRecord.DUPLICATE_READ_FLAG;
			exclFlags |= 0x200; // SAMRecord.READ_FAILS_VENDOR_QUALITY_CHECK_FLAG
		}
		boolean clearDup = main.hasArg("d");
		boolean filterConsecDups = main.hasArg("C");
		boolean filterAllDups = main.hasArg("filter-read-id-dups");
		boolean normalizePairs = main.hasArg("normalize-pairs");
		File headerFile = main.getInputFile("h");
		File rgFile = main.getInputFile("y");
		File suNamesFile = main.getInputFile("u");
		boolean makeUbam = main.hasArg("U");
		long maxRecords = (long)main.getArgDouble("m", 0);
		boolean md5 = main.hasArg("5");
		SortOrder sortOrder = main.getSortOrder();
		boolean setSortOrderEnabled = main.hasArg("set-sort-order");
		if(setSortOrderEnabled) {
			if(sortOrder != null) {
				throw new MainHelper.CommandLineException(
						"Cannot provide both -O and set-sort-order");
			}
			sortOrder = SortOrder.valueOf(
					main.getArgString("set-sort-order"));
			if(filterConsecDups && sortOrder != SortOrder.queryname) {
				throw new MainHelper.CommandLineException(
						"-C requires queryname-sorted input, but your " +
						"--set-sort-order was not queryname");
			}
		}
		
		// Validate ref name filtering options
		if(noRefName && refName != null) {
			throw new MainHelper.CommandLineException("Cannot specify both -c and -n");
		}
		else if(pos == 0) {
			if(allele != 0) {
				throw new MainHelper.CommandLineException("-a requires -p");
			}
			else if(left != 0) {
				throw new MainHelper.CommandLineException("-L requires -p");
			}
			else if(right != 0) {
				throw new MainHelper.CommandLineException("-R requires -p");
			}
		}
		// -p without -c is warning only to allow pre-filtering
		if(refName == null && pos != 0 && !main.hasArg("i")) {
			System.err.println("WARNING -p without -c requires pre-filtering");
		}
		
		// Get and validate output file
		File outputFile = main.getOutputFile("o");
		
		// Open input
		SamReader reader = main.openSamReader("i");
		if(filterConsecDups && !setSortOrderEnabled &&
				reader.getFileHeader().getSortOrder() != SortOrder.queryname) {
			reader.close();
			throw new MainHelper.CommandLineException(
					"-C can only be used with queryname-sorted input");
		}
		
		// Make sure reader is sorted and indexed if doing ref name filtering
		if(hasRefNameFilter && !reader.hasIndex()) {
			reader.close();
			throw new MainHelper.CommandLineException(
					"Input must be indexed bam if using -c, or -n");
		}
		
		// Make worker
		TweakSam worker = new TweakSam(reader);
		
		// Read in read group assignments, if we're doing that
		if(rgFile != null) {
			try {
				worker.enableRgAssignment(rgFile);
			} catch (IOException e) {
				worker.close();
				throw new MainHelper.CommandLineException(
						"Could not read read group assignments file: " + rgFile,
						e);
			}
		}
		
		// Read in names of reads to set to unmapped, if we're doing that
		if(suNamesFile != null) {
			try {
				worker.enableUnmapping(suNamesFile);
			} catch (IOException e) {
				worker.close();
				throw new MainHelper.CommandLineException(
						"Could not read set-unmapped read names file: " +
						rgFile, e);
			}
		}
		
		// Set other options
		worker.setAllele(allele);
		if(hasCrossRefFilter) worker.enableCrossRefFilter(crossRef);
		worker.setReadLenChooser(rlc);
		worker.setRandomReadChooser(rrc);
		worker.setMinMatchRate(minMatchRate);
		worker.setInclFlags(inclFlags);
		worker.setExclFlags(exclFlags);
		worker.setClearDupEnabled(clearDup);
		worker.setFilterConsecDupsEnabled(filterConsecDups);
		worker.setFilterAllDupsEnabled(filterAllDups);
		worker.setNormalizePairsEnabled(normalizePairs);
		worker.setMakeUbamEnabled(makeUbam);
		if(maxRecords > 0) worker.setMaxRecords(maxRecords);
		
		// Open output, with updated header
		worker.openOutputWriter(
				outputFile, headerFile, sortOrder, setSortOrderEnabled, md5);
		
		
		// Initialize the worker, applying any ref/pos criteria
		if(noRefName) {
			worker.initNoRefName();
		}
		else if(refName != null) {
			worker.initLocus(refName, pos, left, right);
		}
		else {
			worker.initWholeInput();
		}
		
		// Do the tweak
		try {
			worker.tweak();
		}
		finally {
			reader.close();
			worker.close();
			System.err.println("Stats:\n" + worker.getCounters());
		}
	}
	
	public void tweak() {
		ReadId prevReadId = null;
		
		Set<ReadId> readIds = new HashSet<ReadId>();
		
		while(iterator.hasNext()) {
			// Get record
			SAMRecord record = iterator.next();
			// Count it
			counters.increment(READ_CTR);
			
			// If doing flag filtering, then see if it passes filters
			if(inclFlags > 0 && (record.getFlags() & inclFlags) == 0) {
				continue;
			}
			if(exclFlags > 0 && (record.getFlags() & exclFlags) > 0) {
				continue;
			}
			
			// Perform duplicate read filtering
			if(filterAllDupsEnabled) {
				ReadId readId = new ReadId(record);
				if(readIds.contains(readId)) {
					counters.increment(DUP_CTR);
					continue;
				}
				else {
					readIds.add(readId);
				}
			}
			else if(filterConsecDupsEnabled) {
				ReadId readId = new ReadId(record);
				if(readId.equals(prevReadId)) {
					counters.increment(DUP_CTR);
					continue;
				}
				else {
					prevReadId = readId;
				}
			}
			
			// If doing cross-reference filtering, then check
			if(hasCrossRefFilter && crossRef == isRecordOneRef(record)) {
				continue;
			}

			// If filtering by read length, then check
			if(rlc != null && !rlc.check(record)) continue;
			
			// If randomly sampling, see if it should be excluded
			if(rrc != null && !rrc.check(record)) continue;
			
			// Assign read group
			if((rgAssignments != null)) {
				String rgid = rgAssignments.get(record.getReadName());
				if(rgid != null) record.setAttribute("RG", rgid);
				else counters.increment(FAILED_RG_CTR);
			}
			
			// If selecting by position, then get the read position, and
			// filter out reads that are not aligned at that position and
			// reads with the wrong allele, if allele filtering is in use
			if(pos != null && pos != 0) {
				// Get the read position (and do left, right, mismatch 
				// filtering)
				int readPos = getReadPos(record);
				
				// Filter out reads not aligned at the position
				if(readPos == 0) continue;
				
				// Filter by allele if specified
				if(allele != 0 &&
						record.getReadString().charAt(readPos - 1) 
						!= allele) continue;
				
				// Add genomic pos and read pos to record
				record.setAttribute("YG", pos);
				record.setAttribute("YR", readPos);
			}
			
			// Clear dup if appropriate
			if(clearDupEnabled && record.getDuplicateReadFlag()) {
				record.setDuplicateReadFlag(false);
			}
			
			// Set unmapped if appropriate
			if((readUnmapper != null) && readUnmapper.hits(record)) {
				ReadUnmapper.setUnmapped(record);
			}
			
			// Strip to unaligned, if appropriate
			if(makeUbamEnabled) {
				record = ReadUnmapper.strip(record, true);
			}
			
			// Write it
			writer.addAlignment(record);
			if(counters.increment(WRITE_CTR) >= maxRecords) break;
		}
		
		if(normalizePairsEnabled) {
			counters.putAll(
					((PairNormalizingSAMFileWriter)writer).getCounters());
		}
	}
	
	/**
	 * Closes the underlying writer
	 */
	public void close() {
		if(writer != null) writer.close();
	}
	
	/**
	 * Class that allows selecting reads by length
	 */
	public static class ReadLenChooser {
		private int min;
		private int max;
		public ReadLenChooser(int min, int max) {
			this.min = min;
			this.max = max;
		}
		
		public boolean check(SAMRecord record) {
			int len = record.getReadLength();
			return min <= len && len <= max;
		}
		
		public static ReadLenChooser parse(String opt) {
			int min = 0;
			int max = Integer.MAX_VALUE;
			boolean open = false;
			
			switch(opt.charAt(0)) {
			case '<':
				max = Integer.parseInt(opt.substring(1)) - 1;
				break;
				
			case '>':
				min = Integer.parseInt(opt.substring(1)) + 1;
				break;
				
			case '[':
				open = true;
				// fall through
				
			case '(':
				// Find comma
				int comma = opt.indexOf(',');
				
				// Get min
				min = Integer.parseInt(opt.substring(1, comma));
				if(open) ++min;
				
				// Get max
				int last = opt.length() - 1;
				open = opt.charAt(last) == ')';
				max = Integer.parseInt(opt.substring(comma + 1, last));
				if(open) --max;
				
				break;
				
			default:
				min = max = Integer.parseInt(opt);
			}
			
			return new ReadLenChooser(min, max);
		}
	}
	
	/**
	 * Class that allows randomly selecting reads in a smart way.
	 * 
	 * You choose initially what proportion of reads should be chosen, and then
	 * you ask about each read using check().  check() should return true for
	 * about the correct proportion of reads, for a sufficiently large set of
	 * reads.  It should also return the same value for two records that have
	 * the same read name (i.e. for both reads in a pair).  Finally, it should
	 * produce different results for different instances.
	 * 
	 * Note that the precision is hard-coded to 997 (the value of PRECISION).
	 */
	public static class RandomReadChooser {
		/**
		 * Hash multiplier to use for read name hashing.
		 */
		public static int MULTIPLIER = 127;
		/**
		 * Number of buckets for hashing
		 */
		public static int PRECISION = 997;
		
		/**
		 * Random shift amount to create different results for every instance.
		 */
		private int shift;
		/**
		 * Cutoff
		 */
		private int cutoff;
		
		/**
		 * Constructs a new random read chooser with the given probability.
		 * 
		 * Note that the effective precision of probability is the 10th of a
		 * percent.
		 * 
		 * @param prob probability in the range (0.0, 1.0)
		 */
		public RandomReadChooser(double prob) {
			this.cutoff = new Double(
					prob * PRECISION).intValue();
			this.shift = new Double(
					Math.floor(Math.random() * PRECISION)).intValue();
		}
		
		/**
		 * Returns true if a read should be chosen
		 * 
		 * @param record the record to check
		 * @return true if the record should be used
		 */
		public boolean check(SAMRecord record) {
			// Get the portion of the read name to hash by trimming off the last
			// character, which sometimes indicates read 1/2
			String toHash = record.getReadName();
			if(toHash.length() > 1) {
				toHash = toHash.substring(0, toHash.length() - 1);
			}
			
			// Hash it to an integer in the range (-inf,inf)
			int rawHash = hash(toHash);
			
			// Map to range [0,PRECISION)
			/*
			 * Note: this does bias against value 0, when there is int overflow
			 * in hash()
			 * 
			 * Here are the ranges, FYI:
			 * 
			 * Value        Range
			 * -----------  ------------
			 * rawHash      (-inf,inf)
			 * % PRECISION  [-PRECISION + 1, PRECISION - 1]
			 * + shift      [-PRECISION + 1, 2 * PRECISION - 1)
			 * + PRECISION  [1, 3 * PRECISION - 1)
			 * % PRECISION  [0, PRECISION)
			 */
			int rangedHash =
				((rawHash % PRECISION) + shift + PRECISION) % PRECISION;
			
			// Check against cutoff
			return rangedHash < cutoff;
		}
		
		/**
		 * Custom hash function that uses MULTIPLIER as the multiplier and 
		 * which should still distribute well for strings that have high
		 * similarity at the beginning
		 * 
		 * @param str string to hash
		 * @return hash value
		 */
		private int hash(String str) {
			int hash = 0;
			byte[] bytes = str.getBytes();
			for(int i = bytes.length - 1; i >= 0; i--) {
				hash = MULTIPLIER * hash + bytes[i];
			}
			return hash;
		}
	}
		
	/**
	 * Returns true if there is a single reference that contains both reads
	 * in the pair (which is true for single-end reads)
	 * 
	 * @param record record to test
	 * @return true unless paired and mate's ref is different
	 */
	private static boolean isRecordOneRef(SAMRecord record) {
		if(!record.getReadPairedFlag()) return true;
		String self = record.getReferenceName();
		String mate = record.getMateReferenceName();
		if(self == mate) return true;
		if(self == null) {
			return mate == null || mate.equals(
					SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
		}
		else if(mate == null) {
			return false;
		}
		else {
			return self.equals(mate) || mate.equals("=");
		}
	}
	
	/**
	 * Returns the position in the read that aligns to the configured
	 * reference position.
	 * 
	 * Returns 0 if the read does not align to the reference at that position
	 * or if it does not satisfy left/right/minMatchRate
	 * 
	 * @param record the record to test
	 * @return 1-based position in the read of the reference position
	 */
	private int getReadPos(SAMRecord record) {
		// Loop through aligned regions
		Blocklist blocklist = BLOCKLIST_FACTORY.newBlocklist(record);
		for(Block block: blocklist) {
			RefInterval refInterval = block.getRefInterval();

			// Only consider aligned blocks 
			if(!block.isAligned()) continue;
				
			// If we got past the position, then the base must not be aligned
			if(refInterval.getLeftInBase() > pos) return 0;
						
			// If we aren't at the position yet, then keep looking
			if(refInterval.getRight() < pos) continue;
			
			// We found the block
			
			// Do left/right checks
			// Note: left/right amount is num bases not including the one at the
			// position.  That, along with the fact that the interval is
			// interbase, and the pos is base, accounts for the particular
			// inequalities below
			if(pos - refInterval.getLeft() <= left) return 0;
			if(refInterval.getRight() - pos < right) return 0;
			
			// Calculate read position
			int readPos = block.getReadInterval().getLeftInt() -
					refInterval.getLeftInt() + pos;
			
			// Do match rate checks
			if(minMatchRate > 0) {
				// Check left side, not including position
				if(!checkMatchRate(
						block,
						0,
						readPos - block.getReadInterval().getLeftInt() - 1)) {
					return 0;
				}
				// Check right side, not including position
				if(!checkMatchRate(
						block,
						readPos - block.getReadInterval().getLeftInt(),
						block.getReadInterval().getRightInt())) {
					return 0;
				}
			}
			
			// Return read position
			return readPos;
		}
		// Never found the position
		return 0;
	}

	/**
	 * Checks the match rate in a read interval
	 * 
	 * @param cigarMd CigarMd built off the SAMRecord
	 * @param left left end of interval to check (interbase)
	 * @param right right end of interval to check (interbase)
	 * @return true if match rate is >= minMatchRate
	 */
	private boolean checkMatchRate(Block block, int left, int right) {
		if(left >= right) return true;
		byte[] readBases = block.getReadBases();
		byte[] refBases = block.getRefBases();
		int matches = 0;
		for(int i = left; i < right; i++) {
			if(readBases[i] == refBases[i]) ++matches;
		}
		int length = right - left + 1;
		double matchRate = ((double)matches) / length; 
		return matchRate >= minMatchRate;
	}
}
