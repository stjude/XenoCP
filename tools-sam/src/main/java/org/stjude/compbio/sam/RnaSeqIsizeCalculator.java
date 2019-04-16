package org.stjude.compbio.sam;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;

import org.apache.commons.cli.Options;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.stjude.compbio.common.RefNameMap;
import org.stjude.compbio.common.formats.refFlat.RefFlatReader;
import org.stjude.compbio.common.formats.seqdict.SeqDictRefNameMapFactory;
import org.stjude.compbio.common.gene.GeneHitGroup;
import org.stjude.compbio.common.gene.GeneHitLibrary;
import org.stjude.compbio.common.gene.LibraryIndexer;
import org.stjude.compbio.common.gene.StdGeneHitLibrary;
import org.stjude.compbio.common.gene.TxHit;
import org.stjude.compbio.common.gene.TxHitGroup;
import org.stjude.compbio.common.gene.TxModel;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.index.RefIntervalIndex;
import org.stjude.compbio.util.MainHelper;
import org.stjude.compbio.util.MainHelper.CommandLineException;
import org.stjude.compbio.util.io.FormatException;

/**
 * Calculates mean and stddev of insert sizes for RNA-Seq data, using refFlat
 * data as a guide.
 */
public class RnaSeqIsizeCalculator {
	/**
	 * GeneHitLibrary to use to find large exons
	 */
	private GeneHitLibrary library;
	/**
	 * Index on the library
	 */
	private RefIntervalIndex<TxHit> locIndex;

	// Options
	
	/**
	 * Read count limit 
	 */
	private int limit = 100000;
	/**
	 * Percent to trim for trimmed mean/stddev
	 */
	private double trimPct = 0.1;
	/**
	 * Minimum exon size to use
	 */
	private int minExonSize = 300;
	/**
	 * Whether to print verbose messages
	 */
	private boolean verbose = false;

	// Used during processing
	private int[] betweens;
	private int[] isizes;
	private int n;

	/**
	 * Constructs a new RnaSeqIsizeCalculator using the given anno library
	 */
	public RnaSeqIsizeCalculator(GeneHitLibrary library) {
		this.library = library;
		LibraryIndexer indexer =
				new LibraryIndexer(false, LibraryIndexer.Key.EXON);
		locIndex = indexer.newTxHitIndex(library);
	}
	
	public int getLimit() {
		return limit;
	}
	public void setLimit(int limit) {
		this.limit = limit;
	}
	public double getTrimPct() {
		return trimPct;
	}
	public void setTrimPct(double trimPct) {
		this.trimPct = trimPct;
	}
	public int getMinExonSize() {
		return minExonSize;
	}
	public void setMinExonSize(int minExonSize) {
		this.minExonSize = minExonSize;
	}
	public boolean isVerbose() {
		return verbose;
	}
	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Statistics reported by getStats()
	 */
	public static enum Stat {
		EXONS_USED, FRAGMENTS, TRIMMED_FRAGMENTS, RETAINED_FRAGMENTS,
		BETWEEN_MEAN, BETWEEN_STDDEV, ISIZE_MEAN, ISIZE_STDDEV
	}
	
	/**
	 * Gets isize stats for a SAM/BAM file
	 * @param reader SamReader to use (must be indexed)
	 * @throws IllegalArgumentException if reader is not indexed
	 */
	public Map<Stat, Double> getStats(SamReader reader) {
		// Assert indexed
		if(!reader.hasIndex()) throw new IllegalArgumentException(
				"Input file must be indexed");

		// Read the ref sequence names into a RefNameMap in case a different
		// convention is used in the BAM than in the refFlat file
		SAMSequenceDictionary seqDict =
				reader.getFileHeader().getSequenceDictionary();
		RefNameMap<String> refNameMap =
				SeqDictRefNameMapFactory.newInstance(seqDict);
		
		// Gather between values
		this.betweens = new int[limit];
		this.isizes = new int[limit];
		this.n = 0;
		int exons = 0;
		for(GeneHitGroup gene: library.getAllGeneHitGroups()) {
			int used = 0;
			for(TxHitGroup txHitGroup: gene.getAllTxHitGroups()) {
				for(TxHit txHit: txHitGroup.getAllTxHits()) {
					TxModel txModel = txHit.getTxModel();
					for(RefInterval exon: txModel.getExons()) {
						if(exon.length() >= minExonSize) {
							used = processBigExon(txModel, exon,
									reader, refNameMap);
							if(used > 0) {
								++exons;
								break;
							}
						}
					}
					if(used > 0) break;
				}
				if(used > 0) break;
			}
			if(n >= betweens.length) break;
		}
		
		// Calculate trimmed mean
		int trim = (int)(trimPct * n);
		int nLessTrim = n - trim;
		SummaryStatistics stats = new SummaryStatistics();
		SummaryStatistics istats = new SummaryStatistics();
		Arrays.sort(betweens, 0, n);
		Arrays.sort(isizes, 0, n);
		for(int i = trim; i < nLessTrim; i++) {
			stats.addValue(betweens[i]);
			istats.addValue(isizes[i]);
		}
		
		// Report
		Map<Stat, Double> out = new HashMap<Stat, Double>();
		out.put(Stat.EXONS_USED, (double)exons);
		out.put(Stat.FRAGMENTS, (double)n);
		out.put(Stat.TRIMMED_FRAGMENTS, 2.0 * trim);
		out.put(Stat.RETAINED_FRAGMENTS, (double)stats.getN());
		out.put(Stat.BETWEEN_MEAN, stats.getMean());
		out.put(Stat.BETWEEN_STDDEV, stats.getStandardDeviation());
		out.put(Stat.ISIZE_MEAN, istats.getMean());
		out.put(Stat.ISIZE_STDDEV, istats.getStandardDeviation());
		return out;
	}

	/**
	 * Validates that a big exon is usable, fetches reads, and adds stats.
	 * 
	 * @param ucscRecord the record from which the exon is taken
	 * @param exon the exon interval
	 * @param samFileReader reader used to retrieve SAM records
	 * @param refNameMap reference name map for the BAM 
	 * @param locIndex location index on the UcscGeneLibrary
	 * @param n number of reads used so far
	 * @param betweens unordered array of between values
	 * @return
	 */
	private int processBigExon(TxModel model, RefInterval exon,
			SamReader samFileReader, RefNameMap<String> refNameMap) {
		// Make sure there isn't another transcript that has an intron here
		List<TxHit> otherTxHits = locIndex.query(exon);
		if(otherTxHits.size() > 1) {
			for(TxHit otherTxHit: otherTxHits) {
				TxModel otherModel = otherTxHit.getTxModel();
				if(otherModel.equals(model)) continue;
				for(RefInterval otherExon: otherModel.getExons()) {
					// If we found an exon that overlaps this one but does not
					// contain it, then this is an ambiguous exon, so do no
					// processing, and return 0 records processed
					if(otherExon.overlaps(exon) && !otherExon.contains(exon)) {
						if(verbose) {
							System.err.println(
									"Rejecting exon " + exon + " in model " +
									model + " because of other exon " + otherExon +
									" in record " + otherTxHit);
						}
						return 0;
					}
				}
			}
		}
		
		// Query for records
		SAMRecordIterator iterator = samFileReader.queryOverlapping(
				refNameMap.get(exon.getRefName()),
				exon.getLeftInBaseInt(),
				exon.getRightInt());
		// Loop over records
		int used = 0;
		while(iterator.hasNext()) {
			SAMRecord record = iterator.next();
			
			// Filter for good records
			
			// Must be mapped with mate mapped
			if(record.getReadUnmappedFlag()) continue;
			if(record.getMateUnmappedFlag()) continue;
			
			// Require mapping quality > 0 (try to get uniquely mapped reads)
			if(record.getMappingQuality() < 1) continue;
			
			// Mate must be on the same chromosome
			if(!record.getMateReferenceName().equals(
					record.getReferenceName())) continue;
			
			// Get left positions in interbase (for testing and use later)
			int left = record.getAlignmentStart() - 1;
			int mateLeft = record.getMateAlignmentStart() - 1;
			
			// Mate must be to the right, but with left end still in the exon
			if(mateLeft <= left || mateLeft >= exon.getRightInt()) continue;
			
			// Done with filtering, now add to stats
			int right = record.getAlignmentEnd();
			int betweenSize = mateLeft - right;
			betweens[n] = betweenSize;
			isizes[n] = betweenSize + 2 * record.getReadLength();
			// Print to stderr if verbose
			if(verbose) System.err.printf(
					"[%8d] %,11d %,11d %,11d %,11d ???? %s\n",
					n, betweenSize, left, right, mateLeft, record.getReadName());
			// Increment and break if necessary
			++used;
			++n;
			if(n >= betweens.length) break;
		}
		
		iterator.close();
		return used;
	}

	/**
	 * Helper for the main method
	 */
	private static SamMainHelper main = new SamMainHelper(getMainOptions());

	private static void usage(int exitcode) {
		main.usage(
				"java " + RnaSeqIsizeCalculator.class.getName() +
				" [OPTION]... BAM\n" +
				"Calculates insert size statistics using read pairs whose\n" +
				"inner edges fall inside large exons.  The exon size should\n" +
				"be greater than the distance in between reads of the\n" +
				"largest inserts you want to be counted.  If it is too\n" +
				"large, though, then the sample size will be small.");
		System.exit(exitcode);
	}

	public static Options getMainOptions() {
		Options options = new Options();
		options.addOption("n", true,
				"max number of inserts to sample (default 100000");
		options.addOption("t", true,
				"portion to trim; range: [0,1), default: 0.1");
		options.addOption("e", true, "min exon size (default 300)");
		options.addOption("v", false, "verbose");
		return options;
	}

	public static void main(String[] args) throws CommandLineException, IllegalStateException, IOException {
		main.addGeneAnnoOption();
		if(args.length == 0) usage(0);
		
		// Parse command line
		main.setArgs(args);
		File annoFile = main.getGeneAnnoFile();
		main.setStrict(false);
		int limit = main.getArgInt("n", 100000);
		double trimPct = main.getArgDouble("t", 0.1);
		int minExonSize = main.getArgInt("e", 300);
		boolean verbose = main.hasArg("v");
		
		// Open refFlat file
		System.err.println("Reading refFlat file");
		StdGeneHitLibrary library;
		try {
			library = new RefFlatReader(annoFile).readAllGeneHits();
		} catch (FormatException e) {
			throw new MainHelper.CommandLineException(e);
		}
		
		// Construct and configure worker
		RnaSeqIsizeCalculator worker = new RnaSeqIsizeCalculator(library);
		worker.setLimit(limit);
		worker.setTrimPct(trimPct);
		worker.setMinExonSize(minExonSize);
		worker.setVerbose(verbose);
		
		// Open BAM file
		File file = main.getInputFile(0);
		System.err.println("Reading from file: " + file);
		SamReader reader = main.openSamReader(file);
		if(!reader.hasIndex()) {
			reader.close();
			throw new MainHelper.CommandLineException(
					"Input file must be indexed");
		}
				

		// Get the stats
		try {
			Map<Stat, Double> stats = worker.getStats(reader);
			
			// Report
			System.out.println("Exons:         " + stats.get(Stat.EXONS_USED));
			System.out.println("Inserts:       " + stats.get(Stat.FRAGMENTS));
			System.out.println("  Trimmed:     " + stats.get(Stat.TRIMMED_FRAGMENTS));
			System.out.println("  Retained:    " + stats.get(Stat.RETAINED_FRAGMENTS));
			System.out.println("Between Mean:  " + stats.get(Stat.BETWEEN_MEAN));
			System.out.println("Between Stdev: " + stats.get(Stat.BETWEEN_STDDEV));
			System.out.println("Isize Mean:    " + stats.get(Stat.ISIZE_MEAN));
			System.out.println("Isize Stdev:   " + stats.get(Stat.ISIZE_STDDEV));
			System.out.println(
					"\nBetween Mean and Stdev are for the 'between sizes';\n"
					+ "that is, the distance between the inner ends of the\n"
					+ "reads (it may be negative).  The insert sizes are\n"
					+ "estimated by adding the between size and twice the\n"
					+ "length of the left-most read.  The right-most read is\n"
					+ "not used as an optimization/simplification.");
		}
		finally {
			reader.close();
		}
	}
}
