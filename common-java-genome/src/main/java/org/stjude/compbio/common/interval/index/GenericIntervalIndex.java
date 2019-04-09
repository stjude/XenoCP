package org.stjude.compbio.common.interval.index;

import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.stjude.compbio.common.interval.Interval;
import org.stjude.compbio.common.interval.tree.Entry;
import org.stjude.compbio.common.interval.tree.IntervalTree;
import org.stjude.compbio.util.HashMultimap;
import org.stjude.compbio.util.Multimap;

/**
 * Generic implementation of LinearIntervalIndex.
 * 
 * In this implementation, the tree rebuild is checked on query, not add/remove,
 * which means that if you subclass this class and override any query methods
 * without calling super, then you are responsible for calling checkRebuild().
 *
 * @param <T> coordinate type
 * @param <I> interval type
 * @param <D> data type
 */
public class GenericIntervalIndex
<T extends Comparable<T>, I extends Interval<T>, D>
implements LinearIntervalIndex<T, I, D> {
	/**
	 * The basis interval tree
	 */
	protected IntervalTree<T, I, D> tree;
	/**
	 * The Multimap used to find entries for a data object
	 */
	protected Multimap<D, Entry<T, I, D>, List<Entry<T, I, D>>> dataMap;
	/**
	 * Tree rebuild strategy 
	 */
	protected TreeRebuildStrategy treeRebuildStrategy;
	
	/**
	 * Constructor using TenTenTreeRebuildStrategy
	 */
	public GenericIntervalIndex() {
		this(new TenTenTreeRebuildStrategy());
	}
	
	/**
	 * Constructor using custom rebuild strategy
	 * 
	 * @param treeRebuildStrategy the tree rebuild strategy to use
	 */
	public GenericIntervalIndex(TreeRebuildStrategy treeRebuildStrategy) {
		// Start with a degenerate tree
		this.tree = new IntervalTree<T, I, D>();
		// And an empty map
		this.dataMap = HashMultimap.newListInstance();
		// Indicated strategy
		this.treeRebuildStrategy = treeRebuildStrategy;
		this.treeRebuildStrategy.setTree(tree);
	}
	
	/**
	 * Checks the rebuild strategy and rebuilds the tree if necessary
	 */
	protected void checkRebuild() {
		if(treeRebuildStrategy.isTreeStale()) {
			tree.rebuild();
			treeRebuildStrategy.rebuilt();
		}
	}
	
	/**
	 * Converts a list of entries to a list of data objects
	 * 
	 * @param entries entries to convert
	 * @return data objects, as list of same length
	 */
	protected List<D> entriesToData(List<Entry<T, I, D>> entries) {
		List<D> out = new LinkedList<D>();
		for(Entry<T, I, D> entry: entries) {
			out.add(entry.getData());
		}
		return out;
	}

	public List<D> query(T point) {
		return entriesToData(queryEntries(point));
	}
	
	public Set<D> queryDistinct(T point) {
		return new HashSet<D>(query(point));
	}

	public List<D> query(I interval) {
		return entriesToData(queryEntries(interval));
	}

	public Set<D> queryDistinct(I interval) {
		return new HashSet<D>(query(interval));
	}

	public void add(I interval, D data) {
		addEntry(new Entry<T, I, D>(interval, data));
	}

	public boolean remove(I interval, D data) {
		return removeEntry(new Entry<T, I, D>(interval, data));
	}

	public List<Entry<T, I, D>> getDataEntries(D data) {
		List<Entry<T, I, D>> entries = dataMap.get(data);
		if(entries == null) entries = new LinkedList<Entry<T, I, D>>();
		return entries;
	}

	public List<Entry<T, I, D>> queryEntries(T point) {
		checkRebuild();
		return tree.stabClosed(point);
	}

	public List<Entry<T, I, D>> queryEntries(I interval) {
		checkRebuild();
		return tree.queryOverlaps(interval);
	}

	public void addEntry(Entry<T, I, D> entry) {
		tree.add(entry);
		dataMap.add(entry.getData(), entry);
		treeRebuildStrategy.added(1);
		scope = (scope == null) ?
				entry.getInterval() : scope.union(entry.getInterval());
	}

	public boolean removeEntry(Entry<T, I, D> entry) {
		if(tree.remove(entry)) {
			dataMap.removeCascade(entry.getData(), entry);
			treeRebuildStrategy.removed(1);
			return true;
		}
		return false;
	}

	public List<Entry<T, I, D>> removeAllEntries(
			Collection<Entry<T, I, D>> entries) {
		List<Entry<T, I, D>> removed = tree.removeAll(entries);
		treeRebuildStrategy.removed(removed.size());
		return removed;
	}

	public int size() {
		return tree.size();
	}

	public Set<D> dataSet() {
		return dataMap.keySet();
	}

	protected Interval<T> scope = null;
	
	/**
	 * Returns the interval that spans from the lowest to highest coordinate
	 * seen.
	 * 
	 * Note: in this implementation, removing the leftmost or rightmost entry
	 * will NOT change the scope.  This may change in the future.  However,
	 * it is guaranteed that there will not be any entries outside of scope.
	 * 
	 * May return null if index is empty.
	 */
	public Interval<T> getScope() {
		return scope;
	}
}
