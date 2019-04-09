package org.stjude.compbio.mapping.standard;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.ListIterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;

import org.apache.commons.cli.Options;
import org.stjude.compbio.common.formats.fasta.IndexedFastaFetcher;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.StdRefInterval;
import org.stjude.compbio.sam.MdElts;
import org.stjude.compbio.sam.SamHeaderHelper;
import org.stjude.compbio.sam.SamMainHelper;
import org.stjude.compbio.sam.block.Block;
import org.stjude.compbio.sam.block.Blocklist;
import org.stjude.compbio.sam.block.BlocklistFactory;
import org.stjude.compbio.sam.block.DirectionalBlockIterator;
import org.stjude.compbio.sam.block.DirectionalBlockIterator.Direction;
import org.stjude.compbio.util.CountingMap;
import org.stjude.compbio.util.JavaHeap;
import org.stjude.compbio.util.MainHelper;

/**
 * Extension of Picard's CleanSam that provides extra functionality for post-
 * mapping cleanup.
 */
public class CleanMappedSam {
	private static final String READ_CTR = "records read";
	private static final String TRIMMED_CTR = "records trimmed";
	private static final String EXTENDED_CTR = "records extended";
	private static final String DIFIXED_CTR = "records with dangling D/I fixed";
	private static final String WRITE_CTR = "records written";

	/**
	 * Reader
	 */
	private final SamReader reader;
	/**
	 * Whether to check for reads that extend beyond the reference
	 */
	private boolean trimEnabled;
	/**
	 * Whether to check for soft-clip that could be converted to alignment
	 */
	private boolean extendEnabled;
	/**
	 * Whether to check for dangling I/D that should be deleted/converted to S
	 */
	private boolean fixDanglingDiEnabled;
	/**
	 * PGID to assign to records (or null to disable)
	 */
	private String pgid;
	/**
	 * Whether to remove paths from PG headers' CL attributes
	 */
	private boolean pgPathCleaningEnabled;
	/**
	 * Whether to remove lowercase cl attribute from PG headers
	 */
	private boolean pgClCleaningEnabled;
	/**
	 * Fetcher for reference sequence when extend is enabled
	 */
	private IndexedFastaFetcher fetcher;
	/**
	 * The writer
	 */
	private SAMFileWriter writer;
	
	/**
	 * BlocklistFactory for trimming and extending
	 */
	private BlocklistFactory factory;
	
	/**
	 * Read counters
	 */
	private final CountingMap<String> counters;

	/**
	 * Whether to print out read names of trimmed/extended records
	 */
	private boolean verbose;
	
	/**
	 * Blocklist for current record
	 */
	private Blocklist blocklist;

	/**
	 * Constructor using reader and trim/extend settings
	 * 
	 * @param reader the reader
	 * @param trimEnabled whether to check for alignment off reference end
	 * @param extendEnabled whether to check for soft clip that could be aligned
	 * @param pgPathCleaningEnabled whether to remove paths from PG headers' CLs
	 * @param pgClCleaningEnabled whether to remove PG headers' lowercase cl
	 */
	public CleanMappedSam(SamReader reader, boolean trimEnabled,
			boolean extendEnabled, boolean fixDanglingDiEnabled,
			boolean pgidEnabled, boolean pgPathCleaningEnabled,
			boolean pgClCleaningEnabled) {
		this.reader = reader;
		this.trimEnabled = trimEnabled;
		this.extendEnabled = extendEnabled;
		this.fixDanglingDiEnabled = fixDanglingDiEnabled;
		this.setPgid(pgidEnabled);
		this.pgPathCleaningEnabled = pgPathCleaningEnabled;
		this.pgClCleaningEnabled = pgClCleaningEnabled;
		this.fetcher = null;
		this.factory = new BlocklistFactory(true, null);
		this.counters = CountingMap.newInitialized(0,
				READ_CTR, TRIMMED_CTR, EXTENDED_CTR, WRITE_CTR);
		this.verbose = false;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Reads the PGID from the reader if PGID is enabled
	 */
	private void setPgid(boolean pgidEnabled) {
		if(!pgidEnabled) {
			this.pgid = null;
			return;
		}
		
		Collection<SAMProgramRecord> pgs =
				SamHeaderHelper.getProgramRecordTips(reader.getFileHeader());
		this.pgid = pgs.isEmpty() ? null : pgs.iterator().next().getId();
	}
	
	/**
	 * Sets reference FASTA to be used for blocklist factory (required for
	 * extending to work)
	 * @throws FileNotFoundException if the file is not found
	 */
	public void setRefFasta(File file) throws FileNotFoundException {
		this.fetcher = (file == null) ? null : new IndexedFastaFetcher(file); 
		factory.setRefFetcher(fetcher);
	}
	
	/**
	 * Sets the writer without performing any header cleaning
	 */
	public void setWriter(SAMFileWriter writer) {
		this.writer = writer;
	}
	
	/**
	 * Builds a clean header and uses it to build the output writer.
	 * 
	 * @return the new writer (so you can close it at the end)
	 */
	public SAMFileWriter setWriter(File outputFile, SortOrder sortOrder) {
		SAMFileHeader header = reader.getFileHeader();
		cleanHeader(header);
		
		boolean presorted = true;
		if(sortOrder != null && header.getSortOrder() != sortOrder) {
			header.setSortOrder(sortOrder);
			presorted = false;
		}
		setWriter(new SAMFileWriterFactory().makeSAMOrBAMWriter(
				header, presorted, outputFile));
		return writer;
	}

	/**
	 * Cleans the header, using settings
	 * @param header header to clean
	 */
	private void cleanHeader(SAMFileHeader header) {
		Pattern pathPattern = Pattern.compile("[^ \t]*/[^ \t]*");
		
		// Clean PG headers
		for(SAMProgramRecord pg: header.getProgramRecords()) {
			// Replace paths in the CL attrib
			if(pgPathCleaningEnabled) {
				String cl = pg.getAttribute("CL");
				if(cl != null) {
					Matcher matcher = pathPattern.matcher(cl);
					pg.setAttribute("CL", matcher.replaceAll("[path]"));
				}
			}
			
			// Remove the cl attrib
			if(pgClCleaningEnabled) {
				pg.setAttribute("cl", null);
			}
		}
		
	}

	public CountingMap<String> getCounters() {
		return counters;
	}

	/**
	 * Cleans the input data and writes to the writer
	 * 
	 * @param writer output writer
	 */
	public void clean(SAMFileWriter writer) {
		for(SAMRecord record: reader) {
			counters.increment(READ_CTR);
			clean(record);
			
			// Write out
			writer.addAlignment(record);
			counters.increment(WRITE_CTR);
		}
	}

	/**
	 * Cleans a single record
	 * @param record record to be cleaned (will be modified)
	 */
	public void clean(SAMRecord record) {
		blocklist = null;
		
		// Only trim/extend if record is mapped
		if(!record.getReadUnmappedFlag()) {
			// Check for trimming and do so if necessary
			if(trimEnabled) checkAndTrim(record);
			
			// Fix dangling D/I before extending
			if(fixDanglingDiEnabled) fixDanglingDi(record);
			
			// Check for soft-clip extension and do so if necessary
			if(extendEnabled) checkAndExtend(record);
		}
		
		// Add PGID if appropriate
		if(pgid != null) record.setAttribute(SAMTag.PG.name(), pgid);
	}

	private void initBlocklist(SAMRecord record) {
		if(blocklist == null) blocklist = factory.newBlocklist(record);
	}
	
	/**
	 * Checks if a record needs to be trimmed due to running off the end of the
	 * reference, and trims it if necessary
	 */
	private void checkAndTrim(SAMRecord record) {
		// Check, and return if OK
	    int seqLen = record.getHeader().getSequence(
	    		record.getReferenceIndex()).getSequenceLength();
		if(record.getAlignmentEnd() <= seqLen) return;
		
		// Not OK, time to trim
		initBlocklist(record);
		
		// Iterate from the right, removing whole blocks until we get
		// inside the reference
		ListIterator<Block> iterator = new DirectionalBlockIterator(
				blocklist, Direction.TO_LEFT_ANCHOR);
		Block block = null;
		int scLength = 0;
		while(iterator.hasNext()) {
			block = iterator.next();
			if(block.isHardTrim()) continue;
			if(block.getRefInterval().getLeft() >= seqLen) {
				scLength += block.getReadInterval().length();
				iterator.remove();
			}
			else {
				break;
			}
		}
		// Now block either spans the ref end or is inside the ref
		RefInterval refInterval = block.getRefInterval();
		if(refInterval.getRight() > seqLen) {
			int totalLen = block.length();
			int inRefLen = seqLen - refInterval.getLeftInt();
			iterator.set(Block.newTemplate(inRefLen, block.getOperator()));
			scLength += totalLen - inRefLen;
		}
		// Back up one spot and add the soft-clip
		iterator.previous();
		iterator.add(Block.newTemplate(scLength, "S"));
		
		// Update the record
		record.setCigar(blocklist.toCigar());
		if(record.getAttribute(SAMTag.MD.name()) != null) {
			record.setAttribute(SAMTag.MD.name(),
					MdElts.listToString(blocklist.toMdEltList()));
		}
		
		if(verbose) System.err.println("Trimmed " + record);
		counters.increment(TRIMMED_CTR);
	}

	/**
	 * Fixes any dangling D/I elements at the end of the alignment.
	 * 
	 * Dangling D are removed, and dangling I are converted to S.
	 */
	private void fixDanglingDi(SAMRecord record) {
		// If there's no D or I, then there's nothing to do
		String cigarString = record.getCigarString();
		if(!(cigarString.contains("D") || cigarString.contains("I"))) return;
		
		// Make sure the blocklist has been created
		initBlocklist(record);

		// We will check in both directions, one at a time
		ListIterator<Block> iterator;
		boolean fixed = false;
		
		// Left side first
		iterator = new DirectionalBlockIterator(
				blocklist, Direction.TO_RIGHT_ANCHOR);
		fixed = fixDanglingDi(iterator) || fixed;
		
		// Then right
		// Get sequence length, so we don't go outside the sequence
		iterator = new DirectionalBlockIterator(
				blocklist, Direction.TO_LEFT_ANCHOR);
		fixed = fixDanglingDi(iterator) || fixed;
		
		if(fixed) {
			// Update the record
			record.setAlignmentStart(blocklist.getAlignmentStart());
			record.setCigar(blocklist.toCigar());
			if(record.getAttribute(SAMTag.MD.name()) != null) {
				record.setAttribute(SAMTag.MD.name(),
						MdElts.listToString(blocklist.toMdEltList()));
			}
			
			if(verbose) System.err.println("Fixed dangling D/I in " + record);
			counters.increment(DIFIXED_CTR);
		}
	}

	/**
	 * Fixes dangling D/I (if found) on one end of a blocklist.
	 *  
	 * @param iterator iterator starting from the end to fix.	 * 
	 * @return whether any change was made
	 */
	private boolean fixDanglingDi(ListIterator<Block> iterator) {
		// Get the first block
		if(!iterator.hasNext()) return false;
		Block block = iterator.next();
		
		// If it's hard trim then skip to next
		while(block.isHardTrim()) {
			if(!iterator.hasNext()) return false;
			block = iterator.next();
		}
		
		// Now, if we're on S, then remember that it exists, and move to next
		Block sclip = null;
		if(block.isClippingOrTrimming()) {
			sclip = block;
			if(!iterator.hasNext()) return false;
			block = iterator.next();
		}
		
		// If we're now on an alignment block, then nothing to do
		if(block.isAligned()) return false;
		
		// Now, process all I/D (and N) until the first M
		int readBasesToSclip = 0;
		while(!block.isAligned()) {
			if(block.getForm().advancesRead) {
				readBasesToSclip += block.length();
			}
			iterator.remove();
			if(!iterator.hasNext()) return false;
			block = iterator.next();
		}
		
		// If there were read bases to sclip, then add/edit the sclip block
		if(readBasesToSclip > 0) {
			// Back up one (we just passed the M block)
			iterator.previous();
			
			// If there is an existing Sclip block, then extend it
			if(sclip != null) {
				sclip = Block.newTemplate(
						sclip.length() + readBasesToSclip, sclip.getOperator());
				iterator.previous();
				iterator.set(sclip);
			}
			// Otherwise, insert one
			else {
				sclip = Block.newTemplate(
						readBasesToSclip, CigarOperator.SOFT_CLIP);
				iterator.add(sclip);
			}
		}
		
		return true;
	}

	/**
	 * Checks if a record can be extended, and does so if possible
	 */
	private void checkAndExtend(SAMRecord record) {
		// If there's no soft-clipping, then there's nothing to do
		if(!record.getCigarString().contains("S")) return;
		
		// Make sure the blocklist has been created
		initBlocklist(record);

		// We will check in both directions, one at a time
		ListIterator<Block> iterator;
		boolean extended = false;
		
		// Left side first
		iterator = new DirectionalBlockIterator(
				blocklist, Direction.TO_RIGHT_ANCHOR);
		extended = checkAndExtend(iterator, null) || extended;
		
		// Then right
		// Get sequence length, so we don't go outside the sequence
	    int seqLen = record.getHeader().getSequence(
	    		record.getReferenceIndex()).getSequenceLength();
		iterator = new DirectionalBlockIterator(
				blocklist, Direction.TO_LEFT_ANCHOR);
		extended = checkAndExtend(iterator, seqLen) || extended;
		
		if(extended) {
			// Update the record
			record.setAlignmentStart(blocklist.getAlignmentStart());
			record.setCigar(blocklist.toCigar());
			if(record.getAttribute(SAMTag.MD.name()) != null) {
				record.setAttribute(SAMTag.MD.name(),
						MdElts.listToString(blocklist.toMdEltList()));
			}
			
			if(verbose) System.err.println("Extended " + record);
			counters.increment(EXTENDED_CTR);
		}
	}

	/**
	 * Checks to see if there is the S-M signature required for extending.
	 * 
	 * If there is, then the iterator will be positioned after the M, so you
	 * will probably want to call previous() to get the M block after this is
	 * done.  If there is not, then the iterator position is not guaranteed.
	 *  
	 * @param iterator iterator positioned at either end of read, oriented in
	 * @param rightSeqLen null if processing left end, sequence length if right
	 * @return amount to extend, or 0 to not extend
	 */
	private boolean checkAndExtend(ListIterator<Block> iterator, Integer rightSeqLen) {
		// Get the first block
		if(!iterator.hasNext()) return false;
		Block block = iterator.next();
		// If it's hard trim then skip to next
		if(block.isHardTrim()) {
			if(!iterator.hasNext()) return false;
			block = iterator.next();
		}
		// Make sure it's soft-clip
		if(block.getOperator() != CigarOperator.S) return false;
		Block s = block;
		// Get the next one, and make sure it's a match
		if(!iterator.hasNext()) return false;
		if(!iterator.next().isAligned()) return false;
		
		// Now, do directional
		return (rightSeqLen == null)
				? checkExtendLeft(iterator, s)
				: checkExtendRight(iterator, s, rightSeqLen);
	}

	/**
	 * Checks how much we can extend left, assuming the S-M signature has
	 * already been found.
	 * 
	 * @param iterator iterator to use for extending, positioned after M
	 * @param sc the soft-clip block
	 * @return amount to extend (0 if we cannot extend)
	 */
	private boolean checkExtendLeft(ListIterator<Block> iterator, Block sc) {
		// Get the match block, and use that to determine S-M boundary
		Block match = iterator.previous();
		long boundary = match.getRefInterval().getLeft();
		
		// Get true maximum extension (cannot extend off end of ref)
		int max = Math.min(sc.length(), (int)boundary);
		
		// Get read/ref bases
		byte[] readBases = Arrays.copyOfRange(
				sc.getReadBases(), sc.length() - max, sc.length());
		byte[] refBases = fetcher.fetchBytes(new StdRefInterval(
				sc.getRefInterval().getRefName(),
				boundary - max, boundary));
		
		// Iterate from right and determine how many are perfect match
		int extendLen = 0;
		for(int i = readBases.length - 1; i >= 0; i--) {
			if(readBases[i] != refBases[i]) break;
			++extendLen;
		}
		
		return extend(iterator, sc, match, extendLen);
	}

	/**
	 * Checks how much we can extend right, assuming the S-M signature has
	 * already been found.
	 * 
	 * @param iterator iterator to use for extending, positioned after M
	 * @param sc the soft-clip block
	 * @param seqLen reference sequence length (also caps extension)
	 * @return amount to extend (0 if we cannot extend)
	 */
	private boolean checkExtendRight(ListIterator<Block> iterator, Block sc,
			int seqLen) {
		// Get the match block, and use that to determine S-M boundary
		Block match = iterator.previous();
		long boundary = match.getRefInterval().getRight();
		
		// Get true maximum extension (cannot extend off end of ref)
		int max = Math.min(sc.length(), seqLen - (int)boundary);
		
		// Get read/ref bases
		byte[] readBases = Arrays.copyOfRange(
				sc.getReadBases(), 0, max);
		byte[] refBases = fetcher.fetchBytes(new StdRefInterval(
				sc.getRefInterval().getRefName(),
				boundary, boundary + max));
		
		// Iterate from left and determine how many are perfect match
		int extendLen = 0;
		for(int i = 0; i < readBases.length; i++) {
			if(readBases[i] != refBases[i]) break;
			++extendLen;
		}
		
		return extend(iterator, sc, match, extendLen);
	}

	/**
	 * Performs extension by a certain amount.
	 * 
	 * Iterator must be positioned between M and S having just read M.
	 * 
	 * @param iterator iterator
	 * @param sc SC block
	 * @param match match block
	 * @param extendLen amount to extend (if 0 false is returned)
	 * @return true if a change was made
	 */
	private boolean extend(ListIterator<Block> iterator, Block sc, Block match,
			int extendLen) {
		if(extendLen == 0) return false;
		
		// Extend the match block
		iterator.set(Block.newTemplate(
				extendLen + match.length(), match.getOperator()));
		
		// Shrink or eliminate the SC block
		iterator.previous();
		if(extendLen == sc.length()) {
			iterator.remove();
		}
		else {
			iterator.set(Block.newTemplate(sc.length() - extendLen, "S"));
		}
		
		return true;
	}

	/**
	 * Helper for the main method
	 */
	private static SamMainHelper main = new SamMainHelper(getMainOptions());

	private static void usage(int exitcode) {
		main.usage(
				"java " + CleanMappedSam.class.getName() +
					" [OPTION]...\n" +
				"Reimplementation of Picard's CleanSam with some extra " +
				"options.");
		System.exit(exitcode);
	}

	public static Options getMainOptions() {
		Options options = new Options();
		options.addOption("i", true, "input sam/bam file");
		options.addOption("o", true, "output sam/bam file");
		options.addOption("t", false, "enabled trimming for alignments that " +
				"run off the end of the reference (bwa has this problem).");
		options.addOption("e", false, "enabled extending alignments by " +
				"identifying soft-clipping that could align in-place " +
				"(requires reference FASTA; bwa has ths problem).");
		options.addOption("d", false, "enable fixing alignments that have " +
				"dangling I/D elements");
		options.addOption("p", false, "assign the PGID to the reads (only " +
				"applicable if there is a PGID in the input header).");
		options.addOption("P", false, "remove paths from PG headers");
		options.addOption("c", false, "remove lowercase cl attribute from PG " +
				"headers (STAR adds cl and CL attributes)");
		return options;
	}

	/**
	 * Main method
	 */
	public static void main(String[] args) throws IOException {
		main.addRefFastaOption();
		main.addSortOrderOption();
		main.addValidationStringencyOption();
		if(args.length == 0) usage(0);
		main.setArgs(args);
		
		// Parse command line
		File outputFile = main.getOutputFile("o");
		boolean trimEnabled = main.hasArg("t");
		boolean extendEnabled = main.hasArg("e");
		boolean fixDanglingDiEnabled = main.hasArg("d");
		boolean pgidEnabled = main.hasArg("p");
		boolean pgPathCleaningEnabled = main.hasArg("P");
		boolean pgClCleaningEnabled = main.hasArg("c");
		main.setStrict(false);
		SortOrder sortOrder = main.getSortOrder();
		File refFasta = main.getRefFastaFile();
		
		// Check extend v refFasta
		if(extendEnabled && refFasta == null) {
			throw new MainHelper.CommandLineException(
					"Cannot perform extension without a reference FASTA file");
		}
		
		JavaHeap.report("after parsing cli");
		
		// Open reader
		SamReader reader = main.openSamReader("i");
		
		JavaHeap.report("after opening SamReader");
		
		// Create the worker
		CleanMappedSam worker = new CleanMappedSam(reader, 
				trimEnabled, extendEnabled, fixDanglingDiEnabled, pgidEnabled,
				pgPathCleaningEnabled, pgClCleaningEnabled);
		worker.setVerbose(true);
		
		JavaHeap.report("after creating worker");
		
		// Set ref fasta
		worker.setRefFasta(refFasta);
		
		JavaHeap.report("after setting reference FASTA");
		
		// Create the writer
		SAMFileWriter writer = worker.setWriter(outputFile, sortOrder);
				
		JavaHeap.report("after creating writer");
		
		// Do the cleansing
		try {
			worker.clean(writer);
		}
		finally {
			JavaHeap.report("end");
			reader.close();
			writer.close();
			System.err.println("Stats:\n" + worker.getCounters());
		}
	}
}
