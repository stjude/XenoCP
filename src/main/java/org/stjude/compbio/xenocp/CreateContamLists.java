package org.stjude.compbio.xenocp;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;

import org.apache.commons.cli.Options;
import org.stjude.compbio.sam.MultiQuerynameSortedIterator;
import org.stjude.compbio.sam.ReadId;
import org.stjude.compbio.sam.SamMainHelper;
import org.stjude.compbio.util.CountingMap;
import org.stjude.compbio.util.io.FormatException;


/**
 * Creates lists of records that appear to be contamination.
 * 
 * A source bam file is provided as well as several contamination bam files.  If
 * a read is mapped better in any of the contamination bams, then it is added to
 * the contamination lists.
 * 
 * For more information, see usage instructions.
 * 
 */
public class CreateContamLists {
	private static final String READ_CTR = "read";
	private static final String CHECK_CTR = "checked";
	private static final String WRITE_CTR = "written";
	private static final String TIE_CTR = "ties";

	private static final String EDIT_DISTANCE_ATTR = "NM";
	
	/**
	 * Multi-iterator for the records aligned against contaminant references
	 */
	private MultiQuerynameSortedIterator contamMultiIterator;

	/**
	 * Read counters 
	 */
	private CountingMap<String> counters;
	
	public CreateContamLists(MultiQuerynameSortedIterator contamMultiIterator) {
		this.contamMultiIterator = contamMultiIterator;
		this.counters = CountingMap.newInitialized(0,
				READ_CTR, CHECK_CTR, WRITE_CTR, TIE_CTR);
	}

	public CountingMap<String> getCounters() {
		return counters;
	}

	/**
	 * Automatically processes all records from a reader, writing contaminant
	 * read names and "tie" records to the writers.
	 * 
	 * @param reader input reader
	 * @param writer contaminant read names writer
	 * @param tieWriter writer for "tie" records (can be null)
	 */
	public void autoProcess(SamReader reader, PrintStream writer,
			SAMFileWriter tieWriter) {
		for(SAMRecord record: reader) {
			counters.increment(READ_CTR);
			
			// Get ReadId
			ReadId readId = ReadId.newInstanceFlex(record);
			
			// Skip if there are no possible contaminant records
			if(!readId.equals(contamMultiIterator.peekReadId())) continue;
			
			// Get the possible contamination records
			List<SAMRecord> contamRecords = contamMultiIterator.nextList();
			
			// Check for contamination
			if(isCheckable(record) && !contamRecords.isEmpty()) {
				counters.increment(CHECK_CTR);
				int cmp = compareToContam(record, contamRecords);
				// Contamination
				if(cmp < 0) {
					write(writer, readId);
					counters.increment(WRITE_CTR);
				}
				// Tie
				else if(cmp == 0) {
					tieWriter.addAlignment(record);
					counters.increment(TIE_CTR);
				}
				// Good read: do nothing
			}
		}
		
		// Make sure we have read all contam records
		SAMRecord nextContam = contamMultiIterator.peek();
		if(nextContam != null) {
			System.err.println(
					"There were unused contaminant records.  The first " +
					"unused record is " + nextContam);
			throw new RuntimeException(
					"Unused contam records; probable ordering problem!");
		}	
	}

	/**
	 * Returns true if an input record should be checked for contamination
	 *  
	 * @param record input record to check
	 * @return true if it should be checked against contams
	 */
	private boolean isCheckable(SAMRecord record) {
		// If it has no mismatch attribute, then we can't check it
		if(record.getAttribute(EDIT_DISTANCE_ATTR) == null) return false;
	
		// Passed all checks, so return true
		return true;
	}

	/**
	 * Compares a primary record to a whole list of contam records.
	 * 
	 * Returns < 0 if a contam record is better, 0 if the best contam record is
	 * a tie, and > 0 if the main record is better.
	 * 
	 * @param record main (original) record
	 * @param contamRecords records mapping the read to contaminant genomes
	 */
	private static int compareToContam(SAMRecord record,
			List<SAMRecord> contamRecords) {
		boolean hasTie = false;
		for(SAMRecord contam: contamRecords) {
			// Skip unmapped contamRecords.  We do this so that unmapped
			// records in original are not counted as ties with unmapped
			// contam records.
			if(contam.getReadUnmappedFlag()) continue;
			// See which is better
			int cmp = cmpRecords(record, contam);
			// Contamination
			if(cmp < 0) {
				return cmp;
			}
			// Tie
			else if(cmp == 0) {
				hasTie = true;
			}
		}
		return hasTie ? 0 : 1;
	}

	/**
	 * Compares two records
	 * 
	 * @param first
	 * @param second
	 * @return 1 if first is better, 0 for tie, -1 if second is better
	 */
	private static int cmpRecords(SAMRecord first, SAMRecord second) {
		// First check mappedness
		if(first.getReadUnmappedFlag()) {
			return second.getReadUnmappedFlag() ? 0 : -1;
		}
		else if(second.getReadUnmappedFlag()) {
			return 1;
		}
		
		// Next, score it
		return score(first) - score(second);
	}

	/**
	 * Calculates the custom score for an alignment
	 * 
	 * @param record the record to score
	 * @return the score (M - num gaps - mismatches)
	 */
	private static int score(SAMRecord record) {
		// Collect stats from CIGAR
		int totalAlignedLength = 0;
		int numGaps = 0;
		int totalGapLength = 0;
		for(CigarElement ce: record.getCigar().getCigarElements()) {
			switch(ce.getOperator()) {
			case M:
				totalAlignedLength += ce.getLength();
				break;
			case I:
			case D:
				++numGaps;
				totalGapLength += ce.getLength();
				break;
			default:
				// No-op
			}
		}
		
		// Get edit distance, which is num mismatches + total gap len
		int editDistance = record.getIntegerAttribute(EDIT_DISTANCE_ATTR);
		
		// Do the final calculation
		/*
		 * Our score is num_matches - num_gaps.  However, we have to do some
		 * work to compute that, and the result is the formula in the return
		 * statement.
		 * 
		 * Proof that the formula gives us num_matches - num_gaps:
		 * 
		 * We know:
		 * total_aligned_length = num_matches + num_mismatches
		 * edit_distance = num_mismatches + total_gap_length
		 * 
		 * So:
		 *   total_aligned_length - edit_distance + total_gap_length - num_gaps
		 * = total_aligned_length - num_mismatches - num_gaps
		 * = num_matches - num_gaps
		 */
		return totalAlignedLength - editDistance + totalGapLength - numGaps; 
	}

	/**
	 * Writes a read id to a PrintStream
	 * 
	 * @param out output PrintStream
	 * @param readId read id to print
	 */
	private static void write(PrintStream out, ReadId readId) {
		out.println(readId.getName());
	}
	
	/**
	 * Helper for the main method
	 */
	private static SamMainHelper main = new SamMainHelper(getMainOptions());

	/**
	 * Prints usage information and exits
	 * @param exitcode exit code (0 for normal)
	 */
	private static void usage(int exitcode) {
		main.usage(
				"java " + CreateContamLists.class.getName() +
				" [OPTION] CONTAM_BAM1,...\n\n"
			+
			"Looks for contamination in a sam/bam file by finding a\n" +
			"better mapping in any of a list of contamination bams.\n" +
			"The main output is a list of read names that are called as\n" +
			"contamination.  You may also use -t to write records that\n" +
			"were called as 'ties' to an output sam/bam file for\n" +
			"review.\n\n" +
			"This will operate on single- or paired-end data, but the\n" +
			"input and contam bams must be 'logically' both single or\n" +
			"or paired.  A read is considered 'logically' paired if it\n" +
			"is actually paired or if it is single end but has a\n" +
			"suffix of [non-alphanumeric char][12] that indicates\n" +
			"read number in the pair.\n\n" +
			"CAVEAT: All bams must have the NM tag, containing edit\n" +
			"distance!\n\n" + 
			"CAVEAT: All contamination bams must be sorted by queryname");
		System.exit(exitcode);
	}

	/**
	 * Returns options for the main method
	 */
	public static Options getMainOptions() {
		Options options = new Options();
		options.addOption("i", true, "input file to use if not stdin");
		options.addOption("o", true, "output contam reads file if not stdout");
		options.addOption("t", true, "output sam/bam to which to write ties");
		return options;
	}

	/**
	 * @param args
	 * @throws IOException 
	 * @throws FormatException 
	 */
	public static void main(String[] args) throws IOException, FormatException {
		main.addValidationStringencyOption();
		if(args.length == 0) usage(0);
		main.setArgs(args);
	
		// Read args
		main.setStrict(false);
		SamReader reader = main.openSamReader("i");
		File outputFile = main.getOutputFile("o");
		File tieFile = main.getOutputFile("t");

		// Open contamination files
		int numContam = main.getArgs().size();
		SamReader[] contamReaders = new SamReader[numContam];
		for(int i = 0; i < numContam; i++) {
			File file = main.getInputFile(i);
			contamReaders[i] = main.openSamReader(file);
		}
		MultiQuerynameSortedIterator contamMultiIterator =
				new MultiQuerynameSortedIterator(contamReaders, true);
		
		// Create the worker
		CreateContamLists worker = new CreateContamLists(contamMultiIterator);
		
		// Open output
		PrintStream writer;
		if(outputFile == null) {
			writer = System.out;
		}
		else {
			try {
				writer = new PrintStream(
						new FileOutputStream(outputFile));
			} catch (FileNotFoundException e) {
				// File was already validated, so this should not happen
				throw new RuntimeException(e);
			}
		}

		// Open tie output
		SAMFileWriter tieWriter = null;
		if(tieFile != null) {
			SAMFileWriterFactory factory = new SAMFileWriterFactory();
			tieWriter = factory.makeSAMOrBAMWriter(
					reader.getFileHeader(), true, tieFile);
		}
		
		// Create the contam lists
		try {
			worker.autoProcess(reader, writer, tieWriter);
		}
		finally {
			reader.close();
			if(outputFile != null) writer.close();
			if(tieWriter != null) tieWriter.close();
			System.err.println(worker.getCounters());
		}
	}

}
