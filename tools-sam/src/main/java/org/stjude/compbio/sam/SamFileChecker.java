package org.stjude.compbio.sam;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import picard.PicardException;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

import org.stjude.compbio.common.formats.fasta.IndexedFastaFetcher;
import org.stjude.compbio.sam.block.Block;
import org.stjude.compbio.sam.block.Blocklist;
import org.stjude.compbio.sam.block.BlocklistFactory;
import org.stjude.compbio.util.Count;
import org.stjude.compbio.util.CountingMap;

/**
 * Reads through a SAM/BAM file and gathers stats that can help with assessing
 * the alignment.
 */
public class SamFileChecker {
	public static enum Metric {
		RECORDS, MAPPED_RECORDS, DUP_RECORDS, SC_RECORDS, SC_BASES,
		SKIP_RECORDS, SKIPS, MATCHED_BASES, ALIGNED_BASES
	}
	
	private BlocklistFactory factory = null;
	private CountingMap<Metric> counters;
	private SortedMap<Integer, long[]> blockLenToMatchCounts;
	private SortedSet<Block> worstBlocks;
	
	private SortedSet<Block> initWorstBlocks() {
		return new TreeSet<Block>(new Comparator<Block>() {
			public int compare(Block o1, Block o2) {
				Count matchCount1 = o1.getMatchCount();
				Count matchCount2 = o2.getMatchCount();
				return (matchCount1.total - matchCount1.count) -
						(matchCount2.total - matchCount2.count);
			}
		});
		
	}

	public void setFactory(BlocklistFactory factory) {
		this.factory = factory;
	}

	public void check(File file) {
		SamReaderFactory.setDefaultValidationStringency(
				ValidationStringency.SILENT);
		check(SamReaderFactory.makeDefault().open(file));
	}
	
	private long[] getMatchCounts(int blockLen) {
		if(!blockLenToMatchCounts.containsKey(blockLen)) {
			blockLenToMatchCounts.put(blockLen, new long[blockLen + 1]);
		}
		return blockLenToMatchCounts.get(blockLen);
	}

	public void check(SamReader reader) {
		this.counters = CountingMap.newInitialized(0, Metric.values());
		this.blockLenToMatchCounts = new TreeMap<Integer, long[]>();
		this.worstBlocks = initWorstBlocks();
		boolean printed0Match = false;
		boolean printed20Match = false;

		for(SAMRecord record: reader) {
			// Basic counts
			counters.increment(Metric.RECORDS);
			if(record.getReadUnmappedFlag()) {
				continue;
			}
			else {
				counters.increment(Metric.MAPPED_RECORDS);
			}
			if(record.getDuplicateReadFlag()) {
				counters.increment(Metric.DUP_RECORDS);
			}
			
			// CIGAR counts
			String cigarStr = record.getCigarString();
			Cigar cigar = record.getCigar();
			
			// SC
			if(cigarStr.contains("S")) {
				counters.increment(Metric.SC_RECORDS);
				for(CigarElement ce: cigar.getCigarElements()) {
					if(ce.getOperator() == CigarOperator.S) {
						counters.increment(Metric.SC_BASES, ce.getLength());
					}
				}
			}
			// N
			if(cigarStr.contains("N")) {
				counters.increment(Metric.SKIP_RECORDS);
				int numSkips = cigarStr.split("N").length - 1;
				counters.increment(Metric.SKIPS, numSkips);
			}
			
			// Matches and mismatches
			if(factory != null) {
				Blocklist blocklist = factory.newBlocklist(record);
				for(Block b: blocklist) {
					if(!b.isAligned()) continue;
					try {
						Count matchCount = b.getMatchCount();
						counters.increment(
								Metric.MATCHED_BASES, matchCount.count);
						counters.increment(
								Metric.ALIGNED_BASES, matchCount.total);
						int len = b.length();
						long[] matchCounts = getMatchCounts(len);
						++matchCounts[matchCount.count];
						
						// Special case prints
						if(matchCount.total >= 80 && matchCount.count < 20) {
							boolean print = false;
							if(matchCount.count == 0 && !printed0Match) {
								print = true;
								printed0Match = true;
							}
							else if(matchCount.count > 0 && !printed20Match) {
								print = true;
								printed20Match = true;
							}
							if(print) {
								System.err.println("LOOK AT " + record + " @ " +
										record.getReferenceName() + ":" +
										record.getAlignmentStart() + " (" +
										b.getCigarElement() + ")");
							}
						}
						
						// Worst blocks
						worstBlocks.add(b);
						if(worstBlocks.size() > 10) {
							worstBlocks.remove(worstBlocks.last());
						}
					} catch(PicardException e) {
						e.printStackTrace();
					}
				}
			}
		}
	}

	public Map<Integer, long[]> getBlockLenToMatchCounts() {
		return blockLenToMatchCounts;
	}

	public CountingMap<Metric> getCounters() {
		return counters;
	}
	
	
	public SortedSet<Block> getWorstBlocks() {
		return worstBlocks;
	}

	private double base = 0;
	
	public char countToAscii(long count) {
		// 0 gets a blank, end of story
		if(count == 0) return ' ';
		
		// Compute what is hopefully a reasonable base to use for log scaling.
		// We want to normalize such that most values are in the range 0-9.
		if(base == 0) {
			double middle = counters.getCount(Metric.RECORDS) /
					(blockLenToMatchCounts.size() / 2);
			base = Math.pow(middle, 0.2);
		}
		
		// Now, use the log, and convert everything >= 10 to X
		double log = Math.floor(Math.log(count) / Math.log(base));
		if(log >= 10) return 'X';
		return Double.toString(log).charAt(0);
	}
	
	
	public static void main(String[] args) throws IOException {
		System.err.println("Usage: Arg 1 = BAM, Arg 2 = FASTA\n");
		
		BlocklistFactory factory = new BlocklistFactory(false, 
				new IndexedFastaFetcher(new File(args[1])));
		
		SamFileChecker checker = new SamFileChecker();
		checker.setFactory(factory);
		
		checker.check(new File(args[0]));
		
		System.out.println("Stats:\n" + checker.getCounters());
		System.out.println();
		System.out.println("Block Len to Aligned Bases:\n");
		System.out.println("Blk  Aligned   1  Bases  2         3         4         5         6         7         8         9         C");
		System.out.println("Len  0123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012");
		for(Map.Entry<Integer, long[]> entry:
				checker.getBlockLenToMatchCounts().entrySet()) {
			System.out.printf("%3d  ", entry.getKey());
			for(long count: entry.getValue()) {
				System.out.print(checker.countToAscii(count));
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("10 Worst Blocks:\n");
		int i = 1;
		for(Block b: checker.getWorstBlocks()) {
			System.out.printf("%2d. %s\n", i++, b.toString());
		}
	}
}
