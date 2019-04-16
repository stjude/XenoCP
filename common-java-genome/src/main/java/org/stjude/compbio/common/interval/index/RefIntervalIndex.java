package org.stjude.compbio.common.interval.index;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.stjude.compbio.common.RefNameMap;
import org.stjude.compbio.common.RefStrand;
import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.RefStranded;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.tree.Entry;

/**
 * An index on ref intervals that takes reference name, and optionally strand,
 * into account.
 *
 * @param <D> data type
 */
public class RefIntervalIndex<D>
implements IntervalIndex<Long, RefInterval, D> {
	/**
	 * Individual indexes for ref name and strand
	 */
	private RefNameMap<Map<Strand, RefStrandIntervalIndex<D>>> indexes;
	/**
	 * Tree rebuild strategy (passed on to each individual index) 
	 */
	private TreeRebuildStrategy treeRebuildStrategy;
	/**
	 * Whether strand use is enabled (if not, then RefStranded keys are all NA)
	 */
	private boolean strandSpecific;
	
	/**
	 * Constructs a new RefIntervalIndex
	 * 
	 * @param strandEnabled whether to provide strand-specific indexes
	 */
	public RefIntervalIndex(boolean strandSpecific) {
		this(new TenTenTreeRebuildStrategy(), strandSpecific);
	}

	/**
	 * Constructs a new RefIntervalIndex with a custom rebuild strategy
	 * 
	 * @param treeRebuildStrategy tree rebuild strategy
	 * @param strandEnabled whether to provide strand-specific indexes
	 */
	public RefIntervalIndex(TreeRebuildStrategy treeRebuildStrategy,
			boolean strandSpecific) {
		this.indexes = new RefNameMap<Map<Strand,
				RefStrandIntervalIndex<D>>>();
		this.treeRebuildStrategy = treeRebuildStrategy;
		this.strandSpecific = strandSpecific;
	}
	
	/**
	 * Gets an immutable set view of ref names represented here
	 */
	public Set<String> refNameSet() {
		return Collections.unmodifiableSet(indexes.keySet());
	}
	
	/**
	 * Gets the index for a RefStranded object (usually a RefInterval)
	 * 
	 * @param refStranded RefStranded object, such as a RefInterval
	 * @return corresponding LongIndex
	 */
	public RefStrandIntervalIndex<D> getIndex(RefStranded refStranded) {
		// Get or create the ref-specific map
		String refName = refStranded.getRefName();
		Map<Strand, RefStrandIntervalIndex<D>> refIndexes;
		if(indexes.containsKey(refName)) {
			refIndexes = indexes.get(refName);
		}
		else {
			refIndexes = initRefIndexMap(refName);
		}

		// Get the index
		Strand strand =
				strandSpecific ? refStranded.getStrand() : Strand.NA;
		return refIndexes.get(strand);
	}

	/**
	 * Initializes all of the strand-specific indexes for a ref name
	 * 
	 * @param refName reference name
	 * @return the new map
	 */
	private Map<Strand, RefStrandIntervalIndex<D>> initRefIndexMap(String refName) {
		// Instantiate the map
		Map<Strand, RefStrandIntervalIndex<D>> refIndexes =
				new HashMap<Strand, RefStrandIntervalIndex<D>>();
		// Add 1 or 3 strand maps, depending on strand specificity option
		Strand[] strands = strandSpecific
				? Strand.values()
				: new Strand[] { Strand.NA };
		for(Strand strand: strands) {
			refIndexes.put(strand, new RefStrandIntervalIndex<D>(
					refName, strand, newTreeRebuildStrategy()));
		}
		// Add to indexes
		indexes.put(refName, refIndexes);
		// Return
		return refIndexes;
	}

	public TreeRebuildStrategy newTreeRebuildStrategy() {
		try {
			return treeRebuildStrategy.getClass().newInstance();
		} catch(Exception e) {
			return new TenTenTreeRebuildStrategy();
		}
	}

	public List<D> query(RefInterval interval) {
		return getIndex(interval).query(interval);
	}
	
	public Set<D> queryDistinct(RefInterval interval) {
		return new HashSet<D>(query(interval));
	}
	
	/**
	 * Performs an stabbing query
	 * @param refName reference name
	 * @param pos position
	 */
	public List<D> query(RefStranded refStranded, long pos) {
		return getIndex(refStranded).query(pos);
	}

	/**
	 * Performs an unstranded stabbing query
	 * @param refName reference name
	 * @param pos position
	 */
	public List<D> query(String refName, long pos) {
		return getIndex(new RefStrand(refName, Strand.NA)).query(pos);
	}

	public void add(RefInterval interval, D data) {
		getIndex(interval).add(interval, data);
	}

	public boolean remove(RefInterval interval, D data) {
		return getIndex(interval).remove(interval, data);
	}

	public List<Entry<Long, RefInterval, D>> getDataEntries(D data) {
		List<Entry<Long, RefInterval, D>> entries =
				new LinkedList<Entry<Long, RefInterval, D>>();
		for(Map<Strand, RefStrandIntervalIndex<D>> refIndexes: indexes.values()) {
			for(RefStrandIntervalIndex<D> index: refIndexes.values()) {
				entries.addAll(index.getDataEntries(data));
			}
		}
		return entries;
	}

	public List<Entry<Long, RefInterval, D>> queryEntries(
			RefInterval interval) {
		return getIndex(interval).queryEntries(interval);
	}

	public void addEntry(Entry<Long, RefInterval, D> entry) {
		getIndex(entry.getInterval()).addEntry(entry);
	}

	public boolean removeEntry(Entry<Long, RefInterval, D> entry) {
		return getIndex(entry.getInterval()).removeEntry(entry);
	}

	public List<Entry<Long, RefInterval, D>> removeAllEntries(
			Collection<Entry<Long, RefInterval, D>> entries) {
		List<Entry<Long, RefInterval, D>> removed =
				new LinkedList<Entry<Long, RefInterval, D>>();
		for(Entry<Long, RefInterval, D> entry: entries) {
			if(removeEntry(entry)) removed.add(entry);
		}
		return removed;
	}

	public int size() {
		int size = 0;
		for(Map<Strand, RefStrandIntervalIndex<D>> refIndexes: indexes.values()) {
			for(RefStrandIntervalIndex<D> index: refIndexes.values()) {
				size += index.size();
			}
		}
		return size;
	}

	public Set<D> dataSet() {
		Set<D> out = new HashSet<D>();
		for(Map<Strand, RefStrandIntervalIndex<D>> refIndexes: indexes.values()) {
			for(RefStrandIntervalIndex<D> index: refIndexes.values()) {
				out.addAll(index.dataSet());
			}
		}
		return out;
	}
}
