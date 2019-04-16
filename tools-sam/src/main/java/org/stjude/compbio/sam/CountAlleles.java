package org.stjude.compbio.sam;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;

import org.apache.commons.cli.Options;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.sam.block.Block;
import org.stjude.compbio.sam.block.Blocklist;
import org.stjude.compbio.sam.block.BlocklistFactory;
import org.stjude.compbio.util.CountingMap;
import org.stjude.compbio.util.MainHelper;
import org.stjude.compbio.util.MainHelper.CommandLineException;
import org.stjude.compbio.util.io.IoField;
import org.stjude.compbio.util.io.IoFormat;
import org.stjude.compbio.util.io.IoRecord;
import org.stjude.compbio.util.io.IoRecordBuilder;
import org.stjude.compbio.util.io.TabularTextReader;
import org.stjude.compbio.util.io.TabularTextReader.Header;
import org.stjude.compbio.util.io.TabularTextWriter;

/**
 * Counts the number of reads in a bam file that support the various alleles at 
 * positions given in an input file
 * 
 * @author mrusch
 */
public class CountAlleles {
	private static final BlocklistFactory BLOCKLIST_FACTORY =
			new BlocklistFactory(false, null);
	
	/**
	 * Input BAM reader
	 */
	private SamReader reader;
	/**
	 * Whether to include dups
	 */
	private boolean includeDupsEnabled = false;
	
	/**
	 * Constructor using a SamReader
	 * 
	 * @param reader open reader
	 */
	public CountAlleles(SamReader reader) {
		if(!reader.hasIndex()) {
			throw new UnsupportedOperationException("Must use indexed BAM");
		}
		this.reader = reader;
	}
	
	public boolean isIncludeDupsEnabled() {
		return includeDupsEnabled;
	}

	public void setIncludeDupsEnabled(boolean includeDupsEnabled) {
		this.includeDupsEnabled = includeDupsEnabled;
	}

	public static final String REF_NAME_KEY = "refName";
	public static final String POS_KEY = "pos";
		
	/**
	 * Pre-processes an input file for processing and returns an output IoFormat.
	 * 
	 * This assumes that any column with a header that starts with "allele"
	 * (case-insensitive) holds an allele to be counted.  The output has
	 * all of the original input columns, followed by specified allele counts,
	 * followed by columns for each possible allele
	 * 
	 * @param reader input reader
	 * @param writer output BufferedWriter
	 * @param refNameField name or 1-based position of ref name field (null for 1)
	 * @param posField name or 1-based position of pos field (null for 2)
	 */
	public static IoFormat<String> preProcessFile(TabularTextReader<String> reader,
			String refNameField, String posField) {
	
		Pattern custPattern = Pattern.compile("^[Aa]llele");
		
		IoFormat<String> inFormat = reader.getIoFormat();
		
		IoFormat<String> outFormat = new IoFormat<String>();
		List<IoField<String>> custFields = new LinkedList<IoField<String>>();
		
		// Loop over input fields
		for(IoField<String> ioField: inFormat.getFields()) {
			// All input fields copied directly to output
			outFormat.addField(new IoField<String>(ioField));
			
			// For custom allele fields, add key to input, and add custom field
			// to be appended to output
			String name = ioField.getName();
			Matcher matcher = custPattern.matcher(name);
			if(matcher.find()) {
				inFormat.addKey(ioField, name);
				custFields.add(new IoField<String>(
						name + "-count", null, name));
			}
		}
		
		// Add fields to output
		for(IoField<String> ioField: custFields) {
			outFormat.addField(ioField);
		}
		for(String alleleStr: new String[] { "C", "G", "T", "A", "N", "-" }) {
			String name = alleleStr;
			if(name.equals("-")) name = "unaligned";
			outFormat.addField(new IoField<String>(
					name + "-count", null, alleleStr));
		}
		
		// Add input special keys
		addSpecialKey(inFormat, REF_NAME_KEY, refNameField);
		addSpecialKey(inFormat, POS_KEY, posField);

		return outFormat;
	}
	
	private static void addSpecialKey(
			IoFormat<String> format, String key, String field) {
		
		IoField<String> ioField;
		try {
			int col = Integer.parseInt(field);
			ioField = format.getOrdinalMap().get(col - 1);
		}
		catch(NumberFormatException e) {
			ioField = format.getNameMap().get(field);
		}
		format.addKey(ioField, key);
	}

	/**
	 * Processes a tab-delimited text file of positions and writes out a new
	 * version annotated with allele counts.
	 * 
	 * This makes use of the keys in the format to determine what to do.  Keys
	 * are as follows:
	 * 
	 * Input:
	 * - REF_NAME_KEY (required) which column holds ref name
	 * - POS_KEY (required) which column holds pos
	 * - all other non-null keys: columns with alleles
	 * 
	 * Output:
	 * - single character keys: allele counts for that allele
	 * - non-null keys matching input: corresponding allele counts
	 *  
	 * @param reader input reader
	 * @param writer output writer
	 * @throws IOException 
	 */
	public void processFile(
			TabularTextReader<String> reader,
			TabularTextWriter<String> writer) throws IOException {
		
		// Get custom allele keys in input
		List<String> custKeys = new LinkedList<String>();
		for(IoField<String> field: reader.getIoFormat().getFields()) {
			String key = field.getKey();
			if(key == null) continue;
			if(key.equals(REF_NAME_KEY) || key.equals(POS_KEY)) continue;
			custKeys.add(key);
		}
		
		// Get specific allele keys in output
		List<String> alleleKeys = new ArrayList<String>(8);
		for(IoField<String> field: writer.getIoFormat().getFields()) {
			String key = field.getKey();
			if(key == null) continue;
			if(key.length() == 1) alleleKeys.add(key);
		}
		
		// Output record builder
		IoRecordBuilder<String> builder =
				new IoRecordBuilder<String>(writer.getIoFormat());
		
		// Process file
		for(IoRecord<String> record: reader) {
			String refName = record.valueForKey(REF_NAME_KEY);
			int pos = Integer.parseInt(record.valueForKey(POS_KEY));
			
			CountingMap<Character> counts = count(refName, pos);
			
			// Initialize output record with the input record
			builder.init(record);
			
			// For custom keys, look up allele for current record
			for(String key: custKeys) {
				String alleleStr = record.valueForKey(key);
				if(alleleStr == null || alleleStr.length() != 1) {
					builder.setValueForKey(key, "");
				}
				else {
					builder.setValueForKey(
							key, getAlleleCountStr(counts, alleleStr));
				}
			}
			
			// For allele keys, the key is the allele
			for(String key: alleleKeys) {
				builder.setValueForKey(
						key, getAlleleCountStr(counts, key));
			}
			
			writer.writeIoRecord(builder.getRecord());
		}
	}
	
	/**
	 * Gets an allele count as a string, where the allele is a string
	 */
	private static String getAlleleCountStr(CountingMap<Character> counts, String alleleStr) {
		return Integer.toString(counts.getCount(alleleStr.charAt(0)));
	}

	/**
	 * Counts alleles at a position
	 * 
	 * @param refName reference name
	 * @param pos position
	 * @return CountingMap for alleles
	 */
	public CountingMap<Character> count(String refName, int pos) {
		CountingMap<Character> out = new CountingMap<Character>();
		
		SAMRecordIterator iterator = reader.queryOverlapping(
				refName, pos, pos);
	
		try {
			while(iterator.hasNext()) {
				SAMRecord record = iterator.next();
				if(!includeDupsEnabled && record.getDuplicateReadFlag()) continue;
				
				int readPos = getReadPos(record, pos);
				
				char allele =
					(readPos == 0)
					? '-'
					: record.getReadString().charAt(readPos - 1);
				
				out.increment(allele);
			}
		}
		finally {
			iterator.close();
		}
		
		return out;
	}

	/**
	 * Returns the position in the read that aligns to the configured
	 * reference position.
	 * 
	 * Returns 0 if the read does not align to the reference at that position
	 * or if it does not satisfy left/right/minMatchRate
	 * 
	 * @param record the record to test
	 * @param pos position
	 * @return 1-based position in the read of the reference position
	 */
	private int getReadPos(SAMRecord record, int pos) {
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
			
			// We found the block; calculate and return read position
			return block.getReadInterval().getLeftInt() -
					refInterval.getLeftInt() + pos;
		}
		// Never found the position
		return 0;
	}



	/**
	 * Helper for the main method
	 */
	private static SamMainHelper main = new SamMainHelper(getMainOptions());

	private static void usage(int exitcode) {
		main.usage(
				"java " + CountAlleles.class.getName() + " [OPTION]... BAM\n" +
				"Reads a BAM file at various positions and reports the counts "
				+ "of reads having various alleles at those positions.\n\n"
				+ "Spec file is tab-delimited text containing columns for:\n"
				+ "1. Reference name (default first column, see --refName)\n"
				+ "2. Position (default first column, see --pos)\n"
				+ "3. Specific alleles you want counted (optional, headers "
				+ "must start with the word 'allele'.\n\n"
				+ "Output file will be the same as input, but with columns "
				+ "appended containing the counts.  First, any custom "
				+ "allele counts (#3 above) are appended, and then allele "
				+ "counts for C, G, T, A, N, and '-' (skips)"
				);
		System.exit(exitcode);
	}

	public static Options getMainOptions() {
		Options options = new Options();
		options.addOption("i", true, "input spec file if not stdin");
		options.addOption("o", true, "output result file if not stdout");
		MainHelper.addLongOpt(options, "ref-name", "NAME",
				"Name or one-based position of column containing ref name");
		MainHelper.addLongOpt(options, "pos", "POS",
				"Name or one-based position of column containing position");
		MainHelper.addLongOpt(options, "include-dups", null,
				"Include marked dups in the counting");
		return options;
	}

	public static void main(String[] args) throws CommandLineException, IllegalStateException, IOException {
		main.addValidationStringencyOption();
		if(args.length == 0) usage(0);
		
		// Parse command line
		main.setArgs(args);
		File bamFile = main.getInputFile(0);
		main.setStrict(false);
		File inputFile = main.getInputFile("i");
		File outputFile = main.getOutputFile("o");
		String refNameField = main.getArgString("ref-name", "1");
		String posField = main.getArgString("pos", "2");
		boolean includeDups = main.hasArg("include-dups");
		
		// Open input
		TabularTextReader<String> reader;
		if(inputFile == null) {
			reader = new TabularTextReader<String>(
					System.in, "\t", Header.required(null));
		}
		else {
			reader = new TabularTextReader<String>(
					inputFile, "\t", Header.required(null));
		}
		
		// Preprocess it
		IoFormat<String> format = CountAlleles.preProcessFile(reader, refNameField, posField);
		
		// Open output
		TabularTextWriter<String> writer;
		if(outputFile == null) {
			writer = new TabularTextWriter<String>(System.out, "\t", format, true); 
		}
		else {
			writer = new TabularTextWriter<String>(outputFile, "\t", format, true); 
		}
		
		// Do the counting
		SamReader bamReader = main.openSamReader(bamFile);
		CountAlleles worker = new CountAlleles(bamReader);
		worker.setIncludeDupsEnabled(includeDups);
		worker.processFile(reader, writer);
		
		// Close everything
		writer.close();
		reader.close();
		bamReader.close();
	}

}
