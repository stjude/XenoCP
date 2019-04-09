package org.stjude.compbio.common.interval.index;

import java.util.List;

import org.stjude.compbio.common.interval.LongInterval;
import org.stjude.compbio.common.interval.tree.Entry;

/**
 * Index using LongIntervals
 *
 * @param <D> data type
 */
public class LongIntervalIndex<D>
extends GenericIntervalIndex<Long, LongInterval, D> {

	/**
	 * Constructs a LongIntervalIndex with default 10-10 rebuild strategy
	 */
	public LongIntervalIndex() {
		super();
	}

	/**
	 * Constructs a LongIntervalIndex with a custom rebuild strategy
	 */
	public LongIntervalIndex(TreeRebuildStrategy treeRebuildStrategy) {
		super(treeRebuildStrategy);
	}

	public List<D> query(int point) {
		return super.query((long)point);
	}
	
	public List<Entry<Long, LongInterval, D>> queryEntries(int point) {
		return super.queryEntries((long)point);
	}

	
}
