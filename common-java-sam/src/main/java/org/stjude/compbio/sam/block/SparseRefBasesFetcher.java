package org.stjude.compbio.sam.block;

import java.util.Arrays;
import java.util.List;

import org.stjude.compbio.common.interval.LongInterval;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.StdLongInterval;
import org.stjude.compbio.common.interval.index.LongIntervalIndex;
import org.stjude.compbio.common.interval.tree.Entry;

/**
 * RefBasesFetcher that uses isolated areas of known ref bases (as could be
 * inferred from an MD tag)
 */
public class SparseRefBasesFetcher implements RefBasesFetcher {
	/**
	 * Reference name
	 */
	private String refName;
	/**
	 * Unstranded (everything is forward) index of known intervals and their
	 * sequences.  Each byte array should be the same length as its interval
	 * key.
	 */
	private LongIntervalIndex<byte[]> knownSeqIntervals;
	/**
	 * Backup fetcher for when the requested interval is outside the sparse
	 * region.
	 */
	private RefBasesFetcher backup;
	
	/**
	 * Constructs a sparse fetcher with the given backup fetcher.
	 *  
	 * @param backup backup fetcher or null for no backup
	 */
	public SparseRefBasesFetcher(String refName, RefBasesFetcher backup) {
		this.refName = refName;
		this.backup = backup;
		this.knownSeqIntervals = new LongIntervalIndex<byte[]>();
	}

	/**
	 * Constructs a sparse fetcher with no backup.
	 */
	public SparseRefBasesFetcher(String refName) {
		this(refName, null);
	}

	/**
	 * Adds known bases starting at a position
	 * 
	 * @param left left position of the run of known bases
	 * @param bases the known bases
	 */
	public void addKnownBases(long left, byte[] bases) {
		knownSeqIntervals.add(
				new StdLongInterval(left, left + bases.length), bases);
	}
	
	public byte[] getRefBases(RefInterval interval) {
		// Check for ref match
		if(!interval.getRefName().equals(refName)) {
			return getRefBasesBackup(interval);
		}
		
		// Query
		List<Entry<Long, LongInterval, byte[]>> entries =
				knownSeqIntervals.queryEntries(interval);
		
		// Make sure there is only one entry, and that it contains the query
		// interval
		if(entries.size() != 1) return getRefBasesBackup(interval);
		LongInterval entryInterval = entries.get(0).getInterval();
		if(!entryInterval.contains(interval)) {
			return getRefBasesBackup(interval);
		}
		
		// Return the bases
		byte[] entryBases = entries.get(0).getData();
		long offset = entryInterval.getLeft();
		return Arrays.copyOfRange(entryBases,
				(int)(interval.getLeft() - offset),
				(int)(interval.getRight() - offset));
	}
	
	private byte[] getRefBasesBackup(RefInterval interval) {
		return backup == null ? null : backup.getRefBases(interval);
	}

}
