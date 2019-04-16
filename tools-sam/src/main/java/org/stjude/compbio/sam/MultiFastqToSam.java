package org.stjude.compbio.sam;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.ReservedTagConstants;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;

import org.apache.commons.cli.Options;
import org.stjude.compbio.common.formats.fastq.FastqHelper;
import org.stjude.compbio.common.formats.fastq.SortedFastqMultiReader;
import org.stjude.compbio.common.formats.fastq.FastqReader;
import org.stjude.compbio.util.MainHelper;
import org.stjude.compbio.util.io.FormatException;

/**
 * Reads in multiple sorted FASTQ files and converts to unaligned BAM files,
 * populating the XX tag with extra related sequences.
 * 
 * IMPORTANT: the input files must have records sorted by query name.
 * 
 * IMPORTANT: input FASTQs must have quality encoded using standard Phred+33
 * encoding, as there is no quality encoding checking in this implementation.
 * 
 * This is analagous to FastqToSam except that it allows specifying multiple
 * FASTQs with related sequences to the primary 1-2 FASTQs, and the matching
 * sequences go in the XX tag.
 */
public class MultiFastqToSam {
	/**
	 * Input multi-reader
	 */
	private SortedFastqMultiReader reader;
	/**
	 * Number of input reads (1 for SE, 2 for PE)
	 */
	private int numInputReads;
	/**
	 * Input labels.
	 * 
	 * Indexes are the same as the matching reader within the multi-reader.
	 * The first 1 (for SE) or 2 (for PE) elements will be ignored
	 */
	private String[] labels;
	/**
	 * Output writer
	 */
	private SAMFileWriter writer;
	
	// Derived from the writer's header
	/**
	 * The header
	 */
	private SAMFileHeader header;
	/**
	 * Read group ID
	 */
	private String rgId;
	
	/**
	 * Constructs a MultiFastqToSam without input/output (they need to be
	 * specified later).
	 */
	public MultiFastqToSam() {
		this.reader = null;
		this.numInputReads = 0;
		this.labels = null;
		this.writer = null;
	}

	/**
	 * Sets up the object to use the given inputs
	 * 
	 * @param file file with SE sequencing data
	 * @param extraFiles map from extra input file to label
	 * @return the multireader
	 * @throws FormatException if the input is not properly formatted
	 * @throws IOException on read error
	 */
	public void setInputsSe(File file, Map<File, String> extraFiles)
	throws IOException, FormatException {
		List<File> files = new ArrayList<File>(1);
		files.add(file);
		setInputs(files, extraFiles);
	}

	/**
	 * Sets up the object to use the given inputs
	 * 
	 * @param file1 file with PE read #1 sequencing data
	 * @param file2 file with PE read #2 sequencing data
	 * @param extraFiles map from extra input file to label
	 * @return the multireader
	 * @throws FormatException if the input is not properly formatted
	 * @throws IOException on read error
	 */
	public void setInputsPe(File file1, File file2,
			Map<File, String> extraFiles)
	throws IOException, FormatException {
		List<File> files = new ArrayList<File>(1);
		files.add(file1);
		files.add(file2);
		setInputs(files, extraFiles);
	}

	/**
	 * Sets up the object to use the given inputs
	 * 
	 * @param file file with SE sequencing data
	 * @param extraFiles map from extra input file to label
	 * @return the multireader
	 * @throws FormatException if the input is not properly formatted
	 * @throws IOException on read error
	 */
	private void setInputs(List<File> files, Map<File, String> extraFiles)
	throws IOException, FormatException {
		this.numInputReads = files.size();
		
		// Initialize labels and build reader/opt indexes list
		int numInputFiles = extraFiles.size() + numInputReads;
		// Init labels
		this.labels = new String[numInputFiles];
		int labelIndex = numInputReads;
		// Init readers
		List<FastqReader> readers = new ArrayList<FastqReader>(numInputFiles);
		for(File file: files) readers.add(new FastqReader(file));
		// Init opt indexes
		List<Integer> optIndexes = new ArrayList<Integer>(extraFiles.size());
		// Fill
		for(Entry<File, String> entry: extraFiles.entrySet()) {
			labels[labelIndex] = entry.getValue();
			readers.add(new FastqReader(entry.getKey()));
			optIndexes.add(labelIndex);
			++labelIndex;
		}
		
		// Construct reader
		this.reader = new SortedFastqMultiReader(readers, optIndexes);
	}
	
	/**
	 * Sets the writer
	 */
	public void setWriter(SAMFileWriter writer) {
		this.writer = writer;
		this.header = writer.getFileHeader();
		List<SAMReadGroupRecord> rgs = header.getReadGroups();
		this.rgId = rgs.isEmpty() ? null : rgs.get(0).getId();
	}

	public void writeSam() {
		// Read and write
		for(FastqRecord[] records: reader) {
			// Build up XX tag
			// Iterate over the extra records
			ExtraSeqs extra = new ExtraSeqs();
			for(int i = numInputReads; i < records.length; i++) {
				// If there was no extra record for this index, then skip it
				if(records[i] == null) continue;
				// Add
				extra.addIfNotEmpty(labels[i], records[i]);
			}
			
			// Write
			for(int i = 0; i < numInputReads; i++) {
				writer.addAlignment(createSamRecord(records[i], i == 1, extra));
			}
		}
	}

	/**
	 * Creates a SAMRecord for a FASTQ record.
	 * 
	 * This method was originally taken from Picard's FastqToSam.
	 * 
	 * @param fqRec FASTQ record
	 * @param second whether this is the second in a pair
	 * @param extra ExtraSeqs for the extra sequences tag, or null
	 * @return new SAMRecord
	 */
	private SAMRecord createSamRecord(FastqRecord fqRec, boolean second,
			ExtraSeqs extraSeqs) {
	    final SAMRecord samRec = new SAMRecord(header);
	    samRec.setReadName(FastqHelper.getName(fqRec));
	    samRec.setReadString(fqRec.getReadString());
	    samRec.setReadUnmappedFlag(true);
	    samRec.setAttribute(ReservedTagConstants.READ_GROUP_ID, rgId);
	    samRec.setBaseQualityString(fqRec.getBaseQualityString());
	    if(numInputReads == 2) {
	        samRec.setReadPairedFlag(true);
	        samRec.setMateUnmappedFlag(true);
	        samRec.setFirstOfPairFlag(!second);
	        samRec.setSecondOfPairFlag(second);
	    }
	    extraSeqs.setSAMRecordAttribute(samRec);
	    return samRec;
	}

	/**
	 * Sets up input information on a MultiFastqToSam object using information
	 * read from a descriptor file.
	 * 
	 * @param worker uninitialized MultiFastqToSam object
	 * @param descriptorReader reader open on the descriptor file
	 * @throws IOException on read error or invalid input file
	 * @throws FormatException if descriptor file or FASTQs are bad.
	 */
	public static void setupFromDescriptor(MultiFastqToSam worker,
			BufferedReader descriptorReader) throws IOException, FormatException {
		Map<String, File> mainFiles = new HashMap<String, File>();
		Map<File, String> extraFiles = new HashMap<File, String>();
		
		// Loop over the descriptor file
		String line;
		while((line = descriptorReader.readLine()) != null) {
			// Parse the line
			String[] parts = line.split("\t");
			File file = new File(parts[0]);
			String label = parts[1];
			
			// Validate the file
			if(!file.isFile()) {
				throw new FileNotFoundException("Not a file: " + file);
			}
			
			// Check for special case labels
			if(label.equals("")) label = "0";
			if(label.matches("[012]")) {
				mainFiles.put(label, file);
			}
			else {
				extraFiles.put(file, label);
			}
		}
		
		// Check for inconsistent se/pe
		boolean isSe = false;
		switch(mainFiles.size()) {
		case 0:
			throw new FormatException(
					"Descriptor must specify main input file(s) " +
					"(label 0 or 1/2)");
		case 1:
			isSe = true;
			if(!mainFiles.containsKey("0")) {
				throw new FormatException(
						"Only one paired-end input file specified");
			}
			break;
		case 2:
			isSe = false;
			if(!mainFiles.containsKey("0")) break;
		case 3:
			throw new FormatException("Both SE and PE input files specified.");
		}
		
		// Do the setup
		if(isSe) {
			worker.setInputsSe(
					mainFiles.get("0"), extraFiles);
		}
		else {
			worker.setInputsPe(
					mainFiles.get("1"), mainFiles.get("2"), extraFiles);
		}
	}

	/**
	 * Builds a header given a list of TAG=VALUE strings
	 * 
	 * @param tagAndValueStrs list of TAG=VALUE strings
	 * @return new header, with a read group if list is non-null/empty
	 */
	public static SAMFileHeader makeHeaderFromTagAndValues(
			List<String> tagAndValueStrs) {
		// Instantiate the header
		SAMFileHeader header = new SAMFileHeader();
		// Always queryname sorted for this class, since inputs must be sorted
		header.setSortOrder(SortOrder.queryname);
		// Make and add the RG if tags/values were given
		SAMReadGroupRecord rg = makeReadGroupFromTagAndValues(tagAndValueStrs);
		if(rg != null) header.addReadGroup(rg);
		// Return
		return header;
	}

	/**
	 * Builds a read gruop record given a list of TAG=VALUE strings
	 * 
	 * @param tagAndValueStrs list of TAG=VALUE strings
	 * @return new read group record, null if list is null/empty
	 */
	public static SAMReadGroupRecord makeReadGroupFromTagAndValues(
			List<String> tagAndValueStrs) {
		// Return null if there are no tags/values
		if(tagAndValueStrs == null || tagAndValueStrs.isEmpty()) return null;
		
		// Parse into map
		Map<String, String> tagToValueMap = new HashMap<String, String>();
		for(String tagAndValueStr: tagAndValueStrs) {
			String[] parts = tagAndValueStr.split("=", 2);
			tagToValueMap.put(parts[0], parts[1]);
		}
		// Handle id
		String id = null;
		if(tagToValueMap.containsKey("RG")) {
			if(tagToValueMap.containsKey("ID")) {
				throw new IllegalArgumentException(
						"Cannot specify both RG and ID tags because RG is an " +
						"alias for ID");
			}
			else {
				id = tagToValueMap.remove("RG");
			}
		}
		else if(tagToValueMap.containsKey("ID")) {
			id = tagToValueMap.remove("ID");
		}
		else {
			tagToValueMap.put("ID", "0");
		}
		
		// Build the read group
		SAMReadGroupRecord rg = new SAMReadGroupRecord(id);
		for(Entry<String, String> entry: tagToValueMap.entrySet()) {
			rg.setAttribute(entry.getKey(), entry.getValue());
		}
		return rg;
	}

	/**
	 * Helper for the main method
	 */
	private static MainHelper main = new MainHelper(getMainOptions());

	private static void usage(int exitcode) {
		main.usage(
				"java " + MultiFastqToSam.class.getName() +
					" [OPTIONS] [TAGS]\n\n"
				+
				"This reads input files from a descriptor file, which may " +
				"also be piped into stdin.  The descriptor file is a two-" +
				"column tabular file with columns:\n" +
				"1. Input file path\n" +
				"2. Label, which is one of:\n" +
				"   * 0/blank for single-end read\n" +
				"   * 1 for paired-end read 1\n" +
				"   * 2 for paired-end read 2\n" +
				"   * <label> to put in XX tag with the given label.\n\n" +
				"The TAGS arguments are all of the form XY=Val, where " +
				"XY is the two-letter tag name, and Val is the value. " +
				"Either RG or ID may be used for the read group's ID. " +
				"These are ignored if -u is used.");
		System.exit(exitcode);
	}

	public static Options getMainOptions() {
		Options options = new Options();
		options.addOption("i", true, "input descriptor file if not stdin");
		options.addOption("o", true, "output BAM file (req'd)");
		options.addOption("u", true, "input ubam from which to read header");
		return options;
	}

	public static void main(String[] args) throws IOException, FormatException {
		if(args.length == 0) usage(0);
		main.setArgs(args);
		
		// Parse command line
		File outputFile = main.getOutputFile("o");
		main.setStrict(false);
		File descriptorFile = main.getInputFile("i");
		File ubamTplFile = main.getInputFile("u");
		List<String> tagAndValueStrs = main.getArgs();
		
		// Open descriptor file/stream
		BufferedReader descriptorReader;
		if(descriptorFile == null) {
			descriptorReader =
					new BufferedReader(new InputStreamReader(System.in));
		}
		else {
			descriptorReader =
					new BufferedReader(new FileReader(descriptorFile));
		}
		
		// Create worker object
		MultiFastqToSam worker = new MultiFastqToSam();
		
		// Read descriptor file and setup inputs
		setupFromDescriptor(worker, descriptorReader);

		// Make header using tag=value strings
		SAMFileHeader header;
		if(ubamTplFile == null) {
			header = makeHeaderFromTagAndValues(tagAndValueStrs);
		}
		else {
			header = SamHeaderHelper.ingestHeader(ubamTplFile);
		}
		
		// Open writer and supply to the worker
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		SAMFileWriter writer =
				factory.makeSAMOrBAMWriter(header, true, outputFile);
		worker.setWriter(writer);
		
		// Do it
		worker.writeSam();
		
		// Close the writer
		writer.close();
	}
}
