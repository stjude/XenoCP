package org.stjude.compbio.sam;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import org.apache.commons.cli.Options;
import org.stjude.compbio.util.DelimitedListBuilder;
import org.stjude.compbio.util.MainHelper.CommandLineException;

/**
 * Given a bam file and a sequence dictionary file, writes a new bam file for
 * the new sequence dictionary.  This will update the header and handle
 * differences in sequence name and order (for coordinate-sorted files)
 * 
 * @author mrusch
 */
public class MigrateSequenceDictionary {
	/**
	 * Input BAM reader
	 */
	private SamReader reader;
	/**
	 * New sequence dictionary
	 */
	private SAMSequenceDictionary newSeqDict;
	/**
	 * Long sequence cutoff
	 */
	private int minLongSeqLen = 0;
	/**
	 * Whether to ignore long sequence errors
	 */
	private boolean ignoreLongSeqErrors = false;

	// Index/cache data
	
	/**
	 * Old sequence dictionary
	 */
	private SAMSequenceDictionary oldSeqDict;
	/**
	 * Length-based index to the (old) seq dict
	 */
	private SortedMap<Integer, SAMSequenceRecord> oldSeqDictIndex;
	/**
	 * New header
	 */
	private SAMFileHeader newHeader;

	// Result of the matching
	
	/**
	 * Set of entries in the old sequence dictionary that are unmatched in new
	 */
	private Set<SAMSequenceRecord> unmatchedOld;
	/**
	 * Set of entries in the new sequence dictionary that are unmatched in old
	 */
	private Set<SAMSequenceRecord> unmatchedNew;
	/**
	 * Mapping from old sequence index to new sequence index
	 */
	private Map<Integer, Integer> seqIndexMap;
	
	/**
	 * Constructor using a SamReader and SAMSequenceDictionary
	 */
	public MigrateSequenceDictionary(SamReader reader, SAMSequenceDictionary newSeqDict) {
		/* No longer a requirement
		if(reader.getFileHeader().getSortOrder() != SortOrder.coordinate) {
			throw new IllegalArgumentException(
					"Input file must be coordinate-sorted (other sorting not "
					+ "yet supported)");
		}
		*/
		
		this.reader = reader;
		this.newSeqDict = newSeqDict;
		initIndexes();
	}

	/**
	 * Initializes the "index" data structures
	 */
	private void initIndexes() {
		// Old sequence dictionary
		oldSeqDict = reader.getFileHeader().getSequenceDictionary();

		// Map from length to old seq dict entry
		this.oldSeqDictIndex = new TreeMap<Integer, SAMSequenceRecord>();
		for(SAMSequenceRecord sequence: oldSeqDict.getSequences()) {
			this.oldSeqDictIndex.put(sequence.getSequenceLength(), sequence);
		}

		// Map from old seq dict index to new seq dict index
		performMatching();
		
		this.newHeader = reader.getFileHeader().clone();
		this.newHeader.setSequenceDictionary(newSeqDict);
	}
	
	/**
	 * Performs matching.
	 * 
	 * Requires:
	 * - seqIndexMap to exist and be empty
	 * - oldSeqDict to be set
	 * - oldSeqDictIndex to be set
	 * 
	 * Will create
	 * - unmatchedOld
	 * - unmatchedNew
	 * - seqIndexMap
	 */
	private void performMatching() {
		// Init outputs
		seqIndexMap = new HashMap<Integer, Integer>();
		unmatchedOld = new HashSet<SAMSequenceRecord>();
		unmatchedNew = new HashSet<SAMSequenceRecord>();

		// Initialized unmatched to all old sequences (remove as matched)
		unmatchedOld.addAll(oldSeqDict.getSequences());
		
		// Loop over the new sequence dictionary and process matches
		for(SAMSequenceRecord newSequence: newSeqDict.getSequences()) {
			// Find the match
			SAMSequenceRecord oldSequence =
					oldSeqDictIndex.get(newSequence.getSequenceLength());
			if(oldSequence == null) {
				unmatchedNew.add(newSequence);
			}
			else {
				unmatchedOld.remove(oldSequence);
				seqIndexMap.put(
						oldSequence.getSequenceIndex(),
						newSequence.getSequenceIndex());
			}
		}
		
		// Process unmatched old sequences
		for(SAMSequenceRecord oldSequence: unmatchedOld) {
			seqIndexMap.put(
					oldSequence.getSequenceIndex(),
					SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX);
		}
	}
	
	public int getMinLongSeqLen() {
		return minLongSeqLen;
	}
	public void setMinLongSeqLen(int minLongSeqLen) {
		this.minLongSeqLen = minLongSeqLen;
	}
	/**
	 * Automatically set the long sequence cutoff.
	 * 
	 * Look for the biggest gap between sizes as calculated using division
	 */
	public void autoSetMinLongSeqLen() {
		float bestGap = 0;
		int bestLength = 0;
		int prevLength = 0;
		
		for(int length: oldSeqDictIndex.keySet()) {
			if(prevLength > 0) {
				float gap = (float)length / prevLength;
				if(gap > bestGap) {
					bestGap = gap;
					bestLength = length;
				}
			}
			prevLength = length;
		}
		
		this.minLongSeqLen = bestLength;
	}
	
	public boolean isIgnoreLongSeqErrors() {
		return ignoreLongSeqErrors;
	}
	public void setIgnoreLongSeqErrors(boolean ignoreLongSeqErrors) {
		this.ignoreLongSeqErrors = ignoreLongSeqErrors;
	}

	public SAMFileHeader getNewHeader() {
		return newHeader;
	}

	/**
	 * Perform validation of the seq dict matches and throw an exception if
	 * there is a problem.
	 * 
	 * @throws InvalidStateException if there is a problem with the headers
	 */
	public void validate() throws IllegalStateException {
		// Initialize list of error messages
		List<String> errMsgs = new LinkedList<String>();
		
		for(SAMSequenceRecord sequence: unmatchedOld) {
			checkUnmatched(sequence, "OLD", errMsgs);
		}
		for(SAMSequenceRecord sequence: unmatchedNew) {
			checkUnmatched(sequence, "NEW", errMsgs);
		}
		if(!errMsgs.isEmpty()) {
			throw new IllegalStateException(
					"The following long sequences were unmatched: "
					+ DelimitedListBuilder.implode(", ", errMsgs));
		}
	}

	/**
	 * Process an unmatched sequence for the validate() method
	 * 
	 * @param sequence
	 * @param oldOrNew
	 * @param errMsgs list of error messages to append iff there is a problem
	 */
	private void checkUnmatched(
			SAMSequenceRecord sequence,
			String oldOrNew,
			List<String> errMsgs) {
		int len = sequence.getSequenceLength();
		if(len >= minLongSeqLen) {
			errMsgs.add(
					sequence.getSequenceName() +
					"(" + oldOrNew + ",len=" + len + ")");
		}
	}
	
	/**
	 * Performs the migration to the given output writer
	 * @param writer open SAM/BAM writer
	 */
	public void migrate(SAMFileWriter writer) {
		// If coordinate-sorted, need to also handle reordering
		if(reader.getFileHeader().getSortOrder() == SortOrder.coordinate) {
			// Loop over the new sequence dictionary and process matches
			for(SAMSequenceRecord newSequence: newSeqDict.getSequences()) {
				// Find the match
				SAMSequenceRecord oldSequence =
						oldSeqDictIndex.get(newSequence.getSequenceLength());
				if(oldSequence == null) {
					System.err.println(
							"No matching sequence for "
							+ newSequence.getSequenceName());
					return;
				}
				
				// Query and migrate everything from the old sequence
				System.err.println(
						"Migrating " + oldSequence.getSequenceName() +
						" to " + newSequence.getSequenceName());
				migrate(writer,
						reader.query(oldSequence.getSequenceName(), 0, 0, false));
			}
			
			// Process unmatched old sequences
			for(SAMSequenceRecord oldSequence: unmatchedOld) {
				System.err.println(
						"Migrating " + oldSequence.getSequenceName() +
						" to unmapped/no-reference");
				migrate(writer,
						reader.query(oldSequence.getSequenceName(), 0, 0, false));
			}
			
			// Process no-ref
			try {
				System.err.println("Copying unmapped/no-reference");
				migrate(writer, reader.queryUnmapped());
			}
			// If there are no no-ref reads, some versions of Picard throw an
			// exception; need to catch and ignore this one
			catch(RuntimeException e) {
				if(!e.getMessage().equals("IOException seeking to unmapped reads")) {
					throw e;
				}
			}
		}
		// Otherwise (not coordinate-sorted), just do the whole file in one shot
		else {
			migrate(writer, reader.iterator());
		}
	}

	/**
	 * Migrates all records in an iterator and closes the iterator
	 * @param writer 
	 * @param iterator
	 */
	private void migrate(SAMFileWriter writer, SAMRecordIterator iterator) {
		while(iterator.hasNext()) {
			migrate(writer, iterator.next());
		}
		iterator.close();
	}

	/**
	 * Migrates a single record
	 * @param writer 
	 * @param record
	 */
	private void migrate(SAMFileWriter writer, SAMRecord record) {
		// Get reference indexes and update the header
		Integer oldReferenceIndex = record.getReferenceIndex();
		Integer oldMateReferenceIndex = record.getMateReferenceIndex();
		record.setHeader(getNewHeader());
		
		// Set new reference index (self)
		if(oldReferenceIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
			Integer newReferenceIndex = seqIndexMap.get(oldReferenceIndex);
			
			// If the new one is no-ref, then need to update a bunch of fields
			if(newReferenceIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
				record.setReadUnmappedFlag(true);
				record.setAttribute(
						MoreSAMTag.ORIG_REF.name(),
						record.getReferenceName());
				record.setReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
				record.setAttribute(
						MoreSAMTag.ORIG_POS.name(),
						record.getAlignmentStart());
				record.setAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
				record.setAttribute(
						MoreSAMTag.ORIG_CIGAR.name(),
						record.getCigarString());
				record.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
			}
			// Otherwise, just update the reference
			else {
				record.setReferenceIndex(newReferenceIndex);
			}
		}
		
		// Set new reference index (mate)
		if(oldMateReferenceIndex != SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
			Integer newMateReferenceIndex =
					seqIndexMap.get(oldMateReferenceIndex);
			
			// If the new one is no-ref, then need to update a bunch of fields
			if(newMateReferenceIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX) {
				record.setMateUnmappedFlag(true);
				record.setMateReferenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
				record.setMateAlignmentStart(SAMRecord.NO_ALIGNMENT_START);
			}
			// Otherwise, just update the reference
			else {
				record.setMateReferenceIndex(newMateReferenceIndex);
			}
		}

		// Write it
		writer.addAlignment(record);
	}

	/**
	 * Helper for the main method
	 */
	private static SamMainHelper main = new SamMainHelper(getMainOptions());

	private static void usage(int exitcode) {
		main.usage(
				"java " + MigrateSequenceDictionary.class.getName() + " [OPTION]... BAM SEQDICT\n" +
				"Writes a new version of the BAM using the new sequence "
				+ "dictionary SEQDICT.  This will update the header, the ref "
				+ "names, and the order (if BAM is coordinate-sorted).\n\n"
				+ "Note: sequence length is used to match up old and new "
				+ "references, and 'long' and 'short' sequences are treated "
				+ "differently, with long sequences treated more strictly. "
				+ "The distinction is made using a length cutoff, which "
				+ "defaults to 0 (treats all as long) and can be calculated "
				+ "automatically by finding the largest gap in size.");
		System.exit(exitcode);
	}

	public static Options getMainOptions() {
		Options options = new Options();
		options.addOption("o", true, "output result file (required)"); // if not stdout");
		options.addOption("l", true, "min long sequence length (default: 0)");
		options.addOption("L", false, "calculate -l using biggest gap");
		options.addOption("f", false, "force; i.e. ignore long ref errors");
		return options;
	}

	public static void main(String[] args) throws CommandLineException, IllegalStateException, IOException {
		main.addValidationStringencyOption();
		if(args.length == 0) usage(0);
		
		// Parse command line
		main.setArgs(args);
		File bamFile = main.getInputFile(0);
		File seqDictFile = main.getInputFile(1);
		File outputFile = main.getOutputFile("o");
		main.setStrict(false);
		Integer minLongLen = main.getArgInt("l");
		boolean autoMinLongLen = main.hasArg("L");
		boolean force = main.hasArg("force");

		// Ingest sequence dictionary file
		SAMSequenceDictionary seqDict =
				SamHeaderHelper
					.ingestHeader(seqDictFile)
					.getSequenceDictionary();

		// Open the input file
		SamReader reader = main.openSamReader(bamFile);

		// Wrap the rest to make sure the reader gets closed
		try {
			// Create the worker and set options
			MigrateSequenceDictionary worker =
					new MigrateSequenceDictionary(reader, seqDict);
			worker.setIgnoreLongSeqErrors(force);
			if(autoMinLongLen) {
				worker.autoSetMinLongSeqLen();
				System.err.println("Automatically set -l to " +
						worker.getMinLongSeqLen());
			}
			else if(minLongLen != null) {
				worker.setMinLongSeqLen(minLongLen);
			}
	
			// Unless force was set, perform validation
			if(!force) {
				try {
					worker.validate();
				}
				catch(IllegalStateException e) {
					System.err.println("Validation failed: " + e.getMessage());
					System.err.println(
							"Use -f or a new -l if you want to continue anyway");
					System.exit(1);
				}
			}
			
			// Open the output writer
			// TODO Remove this hardcoding and move to SamMainHelper
			//SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
			SAMFileWriterFactory.setDefaultCreateMd5File(true);
			SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(
					worker.getNewHeader(), true, outputFile);
			
			
			// Do the work
			try {
				worker.migrate(writer);
			}
			finally {
				// Close writer
				writer.close();
			}
		}
		finally {
			reader.close();
		}
	}

}
