package org.stjude.compbio.mapping.standard;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;

import org.apache.commons.cli.Options;
import org.stjude.compbio.common.RefNameMap;
import org.stjude.compbio.common.formats.seqdict.SeqDictRefNameMapFactory;
import org.stjude.compbio.common.seq.SequenceHelper;
import org.stjude.compbio.sam.ReadId;
import org.stjude.compbio.sam.SamMainHelper;
import org.stjude.compbio.util.CountingMap;
import org.stjude.compbio.util.MainHelper;

/**
 * Simplified version of Picard's MergeBamAlignment.
 * 
 * This version requires that the BAMs have the same set of records in the same
 * order.
 * 
 * It uses all of the information from the aligned BAM, and copies read groups,
 * and unused attributes from the unaligned BAM.
 */
public class SimpleMergeBamAlignment {
	private static final String ALN_READ_CTR = "aligned records read";
	private static final String UBAM_READ_CTR = "unaligned records read";
	private static final String WRITE_CTR = "records written";

	public static final ExtraOption DEFAULT_EXTRA_OPTION = ExtraOption.ERROR;
	public static enum ExtraOption { ERROR, TRIM, ALL }
	
	/**
	 * Aligned BAM reader
	 */
	private final SamReader alignedReader;
	/**
	 * Unaligned BAM reader
	 */
	private final SamReader ubamReader;
	/**
	 * Sequence dictionary reader (may be null)
	 */
	private final SamReader seqDictReader;
	/**
	 * PG reader (may be null)
	 */
	private final SamReader pgReader;
	/**
	 * PGID to set on RGs (may be null)
	 */
	private final String pgId;
	/**
	 * ExtraOption
	 */
	private final ExtraOption extraOption;
	/**
	 * BLAT mode (assumes input is SAM from psl2sam.pl)
	 */
	private boolean blatMode;
	/**
	 * Ignore the non-primary flag
	 */
	private boolean ignoreNonPrimary;
	
	/**
	 * SAMFileHeader built from input
	 */
	private SAMFileHeader header;
	/**
	 * RefNameMap for canonicalization, built off the header
	 */
	private RefNameMap<String> refNameMap;
	/**
	 * Determination of whether or not the input is pre-sorted
	 */
	private boolean presorted;
	
	/**
	 * Read counters
	 */
	private final CountingMap<String> counters;
	
	/**
	 * Constructor using SamReaders
	 * 
	 * @param alignedReader aligned BAM reader
	 * @param ubamReader unaligned BAM reader
	 * @param seqDictReader seq dict reader (may be null)
	 * @param pgReader PG header reader (may be null)
	 * @param pgId PG ID to set on RGs (may be null)
	 */
	public SimpleMergeBamAlignment(SamReader alignedReader,
			SamReader ubamReader, SamReader seqDictReader,
			SamReader pgReader, String pgId) {
		this(alignedReader, ubamReader, seqDictReader, pgReader, pgId,
				DEFAULT_EXTRA_OPTION);
	}

	/**
	 * Constructor using SamReaders
	 * 
	 * @param alignedReader aligned BAM reader
	 * @param ubamReader unaligned BAM reader
	 * @param seqDictReader seq dict reader (may be null)
	 * @param pgReader PG header reader (may be null)
	 * @param pgId PG ID to set on RGs (may be null)
	 * @param extraOption what to do if ubam has extra records
	 */
	public SimpleMergeBamAlignment(SamReader alignedReader,
			SamReader ubamReader, SamReader seqDictReader,
			SamReader pgReader, String pgId, ExtraOption extraOption) {
		this.alignedReader = alignedReader;
		this.ubamReader = ubamReader;
		this.seqDictReader = seqDictReader;
		this.pgReader = pgReader;
		this.pgId = pgId;
		this.extraOption = extraOption;
		this.header = null;
		this.presorted = true;
		this.counters = CountingMap.newInitialized(0, ALN_READ_CTR, WRITE_CTR);
	}

	/**
	 * Sets whether to operate in blat mode
	 */
	public void setBlatMode(boolean blatMode) {
		this.blatMode = blatMode;
	}

	/**
	 * Sets whether the non-primary flag should be ignored
	 */
	public void setIgnoreNonPrimary(boolean ignoreNonPrimary) {
		this.ignoreNonPrimary = ignoreNonPrimary;
	}

	public SAMFileHeader getHeader() {
		return header;
	}

	public void setHeader(SAMFileHeader header) {
		this.header = header;
		this.refNameMap = SeqDictRefNameMapFactory.newInstance(
				header.getSequenceDictionary());
	}

	public boolean isPresorted() {
		return presorted;
	}

	public void setPresorted(boolean presorted) {
		this.presorted = presorted;
	}

	public CountingMap<String> getCounters() {
		return counters;
	}

	/**
	 * Uses the readers and an optional sort order to build the output header
	 * and determine if output is in the same sort order as input.
	 * 
	 * The resulting header and sort order may be fetched using getters
	 * 
	 * @param sortOrder optional output sort order to use, null to keep aligned
	 * 					sort order
	 */
	public void buildHeader(SortOrder sortOrder) {
		// Use aligned file header to start out with
		SAMFileHeader header = alignedReader.getFileHeader();
		presorted = true;
		
		// Reset sort order, if a different one was specified
		if(sortOrder != null && header.getSortOrder() != sortOrder) {
			header.setSortOrder(sortOrder);
			presorted = false;
		}
		
		// Copy read groups from ubam
		header.setReadGroups(ubamReader.getFileHeader().getReadGroups());
		
		// Replace seq dict with one from seq dict file if specified
		if(seqDictReader != null) {
			header.setSequenceDictionary(
					seqDictReader.getFileHeader().getSequenceDictionary());
			try {
				seqDictReader.close();
			}
			catch(IOException e) {
				// No-op
			}
		}
		
		// Add program groups from PG file if specified
		if(pgReader != null) {
			for(SAMProgramRecord pg:
					pgReader.getFileHeader().getProgramRecords()) {
				header.addProgramRecord(pg);
			}
		}
		
		// Set PGID in read groups if specified
		if(pgId != null) {
			for(SAMReadGroupRecord rg: header.getReadGroups()) {
				rg.setAttribute("PG", pgId);
			}
		}
		
		// Set it (this also sets up the ref name map)
		this.setHeader(header);
	}
	
	/**
	 * Performs the merge of the input readers to a new output writer open on
	 * the indicated file.
	 * 
	 * @param outputFile output file
	 */
	public void merge(File outputFile) {
		SAMFileWriterFactory factory = new SAMFileWriterFactory();
		SAMFileWriter writer = factory.makeSAMOrBAMWriter(
				header, presorted, outputFile);
		try {
			merge(writer);
		}
		finally {
			writer.close();
		}
	}
	
	/**
	 * Performs the merge of the input readers to the output writer.
	 * 
	 * @param writer output writer
	 */
	public void merge(SAMFileWriter writer) {
		// All option is special
		if(extraOption == ExtraOption.ALL) {
			mergeAllUbam(writer);
			return;
		}

		// Whether or not we're in trim mode
		boolean trim = extraOption == ExtraOption.TRIM;
		
		// Loop over aligned records
		SAMRecordIterator ubamIterator = ubamReader.iterator();
		ReadId readId;
		SAMRecord ubamRecord = null;
		ReadId ubamReadId = null;
		for(SAMRecord record: alignedReader) {
			counters.increment(ALN_READ_CTR);
			readId = getReadId(record);
			
			// Get the matching ubam record
			if(ignoreNonPrimary || !record.getNotPrimaryAlignmentFlag()) {
				do {
					if(!ubamIterator.hasNext()) {
						throw new RuntimeException("Ran out of ubam records");
					}
					ubamRecord = ubamIterator.next();
					ubamReadId = getReadId(ubamRecord);
					counters.increment(UBAM_READ_CTR);
				} while(trim && !readId.equals(ubamReadId));
			}
			
			// Verify that they are the same read
			if(ubamRecord == null || !readId.equals(ubamReadId)) {
				System.err.println("Records are out of order:");
				System.err.println("Aligned record:   " + record);
				System.err.println("Unaligned record: " + ubamRecord);
				throw new RuntimeException("Records are out of order");
			}
			
			// Copy attributes from ubam record
			copyAttributes(record, ubamRecord);
			
			// Do BLAT-mode processing
			if(blatMode) doBlatModeMerge(record, ubamRecord, readId);
			
			// Write out
			write(writer, record);
		}
	}

	/**
	 * Gets a ReadId for a record.
	 * 
	 * If in blatMode, this uses ReadId.newInstanceFlex() instead of the normal
	 * constructor.
	 * 
	 * @param record the record
	 * @return the ReadId
	 */
	private ReadId getReadId(SAMRecord record) {
		if(record == null) return null;
		return blatMode ? ReadId.newInstanceFlex(record) : new ReadId(record);
	}

	/**
	 * Writes a record, applying canonicalization of ref name, and updating
	 * counter
	 * 
	 * @param writer writer to use
	 * @param record record to write
	 */
	private void write(SAMFileWriter writer, SAMRecord record) {
		// Canonicalize reference name
		String refName = record.getReferenceName();
		if(refName != null &&
				!refName.equals(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME)) {
			String canonRefName = refNameMap.get(refName);
			if(!refName.equals(canonRefName)) {
				record.setReferenceName(canonRefName);
			}
		}
		// Write
		writer.addAlignment(record);
		// Count
		counters.increment(WRITE_CTR);
	}

	/**
	 * Merge implementation when all records from the ubam should be used.
	 * @param writer
	 */
	private void mergeAllUbam(SAMFileWriter writer) {
		// Aligned iterator
		SAMRecordIterator iterator = alignedReader.iterator();
		// Current aligned record
		SAMRecord record = iterator.hasNext() ? iterator.next() : null;
		if(record != null) counters.increment(ALN_READ_CTR);
		ReadId readId = getReadId(record);

		// Loop over ubam records
		ReadId ubamReadId;
		for(SAMRecord ubamRecord: ubamReader) {
			counters.increment(UBAM_READ_CTR);
			ubamReadId = getReadId(ubamRecord);
			
			// See if it matches the aligned record
			if(record != null && readId.equals(ubamReadId)) {
				// Copy attributes from ubam record
				copyAttributes(record, ubamRecord);
				// Do BLAT-mode processing
				if(blatMode) doBlatModeMerge(record, ubamRecord, readId);
				
				// Write it
				write(writer, record);
				
				// Go to next aligned record
				record = iterator.hasNext() ? iterator.next() : null;
				if(record != null) counters.increment(ALN_READ_CTR);
				readId = getReadId(record);
			}
			// If not, write the ubam unchanged
			else {
				write(writer, ubamRecord);
			}
		}

	}

	/**
	 * Copies attributes from a ubam record to a record.
	 * 
	 * @param record recipient record
	 * @param ubamRecord ubam record as source
	 */
	private void copyAttributes(SAMRecord record, SAMRecord ubamRecord) {
		for(SAMTagAndValue tagAndValue: ubamRecord.getAttributes()) {
			if(record.getAttribute(tagAndValue.tag) == null) {
				record.setAttribute(tagAndValue.tag, tagAndValue.value);
			}
		}
	}

	/**
	 * Performs BLAT-mode merge steps.
	 * 
	 * This is usually useful for aligned SAMs from psl2sam.pl
	 * 
	 * Copies sequence and quality from a ubam record to a record, if they are
	 * missing.  Also converts H to S in Cigar and sets pairing flags.
	 * 
	 * @param record recipient record
	 * @param ubamRecord ubam record as source
	 * @param readId record's ReadId
	 */
	private void doBlatModeMerge(SAMRecord record, SAMRecord ubamRecord,
			ReadId readId) {
		// If the read mapped to the reverse strand, then we need to RC it
		// because the SAM file always presents the reference forward strand.
		boolean reverse = record.getReadNegativeStrandFlag();
		// Seq
		if(record.getReadBases().length == 0) {
			String seq = ubamRecord.getReadString();
			if(reverse) seq = SequenceHelper.reverseComplement(seq);
			record.setReadString(seq);
		}
		// Qual
		String qual = ubamRecord.getBaseQualityString();
		if(reverse) qual = new StringBuilder(qual).reverse().toString();
		record.setBaseQualityString(qual);
		// Change H to S in CIGAR
		record.setCigarString(record.getCigarString().replace('H', 'S'));
		// Set pairing flags
		if(readId.isPaired()) {
			record.setReadName(readId.getName());
			record.setFlags(record.getFlags() | readId.getPairingFlags());
		}
	}

	/**
	 * Helper for the main method
	 */
	private static SamMainHelper main = new SamMainHelper(getMainOptions());

	private static void usage(int exitcode) {
		main.usage(
				"java " + SimpleMergeBamAlignment.class.getName() +
					" [OPTION]...\n" +
				"Simple version of Picard's MergeBamAlignment.  This " +
				"version requires that the BAMs have the same set of records " +
				"in the same order (but see -e).  It uses all of the information " +
				"from the aligned BAM, and copies read groups, and unused " +
				"attributes from the unaligned BAM.\n\n" +
				"If you have tags on the unmapped reads in the aligned BAM " +
				"that you want to keep, then this version is a good " +
				"option.\n\n" +
				"The -a, -u, and -o arguments are required.  If you " +
				"specify a sequence dictionary file using -d, then the " +
				"sequence dictionary will be read from that SAM-formatted " +
				"file, otherwise, it is taken from the aligned BAM file.\n\n" +
				"If you want to specify program groups, then use -p with a " +
				"SAM-formatted file containing program group headers.  The " +
				"-P option may be used to give the program group ID that " +
				"should be set in all of the read group headers' PG " +
				"attribute.");
		System.exit(exitcode);
	}

	public static Options getMainOptions() {
		Options options = new Options();
		options.addOption("a", true, "input aligned sam/bam file");
		options.addOption("u", true, "input ubam file");
		options.addOption("o", true, "output sam/bam file");
		options.addOption("d", true, "sequence dictionary file");
		options.addOption("p", true, "program group file");
		options.addOption("P", true, "program group ID for read groups");
		options.addOption("b", false,
				"treat input as sam from psl2sam.pl run on a BLAT result");
		options.addOption("e", true, "action to perform if ubam/aligned " +
				"records mismatch (assuming extra ubam).  Options are " +
				"'ERROR' (default), 'TRIM' do not output unmatched ubam " +
				"records, and 'ALL' output all ubam records");
		MainHelper.addLongOpt(options, "ignore-non-primary-flag", null,
				"Ignore the non-primary flag; use this for mappings produced " +
				"by BSMAP where the non-primary flag is overloaded for " +
				"another purpose");
		return options;
	}

	/**
	 * Main method
	 */
	public static void main(String[] args) throws IOException {
		main.addValidationStringencyOption();
		main.addSortOrderOption();
		if(args.length == 0) usage(0);
		main.setArgs(args);
		
		// Parse command line
		SamReader alignedReader = main.openSamReader("a");
		SamReader ubamReader = main.openSamReader("u");
		File outputFile = main.getOutputFile("o");
		main.setStrict(false);
		File seqDictFile = main.getInputFile("d");
		File pgFile = main.getInputFile("p");
		String pgId = main.getArgString("P");
		ExtraOption extraOption = main.getArgEnum("e", ExtraOption.class);
		SortOrder sortOrder = main.getSortOrder();
		boolean blatMode = main.hasArg("b");
		boolean ignoreNonPrimary = main.hasArg("ignore-non-primary-flag");
	
		// Open readers
		SamReader seqDictReader = null;
		if(seqDictFile != null) {
			seqDictReader = main.openSamReader(seqDictFile);
		}
		SamReader pgReader = null;
		if(pgFile != null) {
			pgReader = main.openSamReader(pgFile);
		}
		
		// Create the SimpleMergeBamAlignment object
		if(extraOption == null) extraOption = DEFAULT_EXTRA_OPTION;
		SimpleMergeBamAlignment merger = new SimpleMergeBamAlignment(
				alignedReader, ubamReader, seqDictReader, pgReader, pgId,
				extraOption);
		merger.setBlatMode(blatMode);
		merger.setIgnoreNonPrimary(ignoreNonPrimary);
		
		// Build the header
		merger.buildHeader(sortOrder);
		
		// Do the merge
		try {
			merger.merge(outputFile);
		}
		finally {
			alignedReader.close();
			ubamReader.close();
			System.err.println("Stats:\n" + merger.getCounters());
		}
	}
}
