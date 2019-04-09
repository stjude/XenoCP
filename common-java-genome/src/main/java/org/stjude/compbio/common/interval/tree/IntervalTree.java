package org.stjude.compbio.common.interval.tree;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import org.stjude.compbio.common.interval.Interval;

/**
 * An interval tree
 *
 * @param <T> coordinate type
 * @param <I> interval type
 * @param <D> data type
 */
public class IntervalTree<T extends Comparable<T>, I extends Interval<T>, D> {
	/**
	 * Tree root node
	 */
	private Node<T, I, D> root;

	/**
	 * Build a tree from a collection of entries
	 */
	public static <T extends Comparable<T>, I extends Interval<T>, D>
	IntervalTree<T, I, D> build(List<Entry<T, I, D>> entries) {
		return new IntervalTree<T, I, D>(new Node<T, I, D>(entries));
	}
	
	/**
	 * Builds an empty, degenerate tree
	 */
	public IntervalTree() {
		this.root = new Node<T, I, D>(new LinkedList<Entry<T, I, D>>());
	}
	
	/**
	 * Construct a tree using a root node.
	 * 
	 * @param root the root node
	 */
	private IntervalTree(Node<T, I, D> root) {
		this.root = root;
	}
	
	/**
	 * Rebuilds the tree
	 */
	public void rebuild() {
		List<Entry<T, I, D>> allEntries = new LinkedList<Entry<T, I, D>>();
		root.getAllEntries(allEntries);
		this.root = new Node<T, I, D>(allEntries);
	}
	
	/**
	 * Adds an entry to the tree
	 */
	public void add(Entry<T, I, D> entry) {
		root.add(entry);
	}
	
	/**
	 * Adds an entry to the tree
	 * 
	 * @param interval the interval
	 * @param data the data object to add there
	 */
	public void add(I interval, D data) {
		add(new Entry<T, I, D>(interval, data));
	}
	
	/**
	 * Adds all entries in a collection to the tree
	 */
	public void addAll(Collection<Entry<T, I, D>> entries) {
		for(Entry<T, I, D> entry: entries) {
			add(entry);
		}
	}
	
	/**
	 * Removes an entry from the tree
	 * @return true if the entry was found and removed
	 */
	public boolean remove(Entry<T, I, D> entry) {
		return root.remove(entry);
	}
	
	/**
	 * Removes an entry from the tree
	 * 
	 * @param interval the interval
	 * @param data the data object to add there
	 * @return true if the entry was found and removed
	 */
	public boolean remove(I interval, D data) {
		return remove(new Entry<T, I, D>(interval, data));
	}
	
	/**
	 * Removes all entries in a collection from the tree
	 * @return removed entries
	 */
	public List<Entry<T, I, D>> removeAll(Collection<Entry<T, I, D>> entries) {
		List<Entry<T, I, D>> removed = new LinkedList<Entry<T, I, D>>();
		for(Entry<T, I, D> entry: entries) {
			if(remove(entry)) removed.add(entry);
		}
		return removed;
	}

	/**
	 * Performs a stabbing query, treating intervals as open
	 * 
	 * @param point query point
	 * @return list of hits as entries
	 */
	public List<Entry<T, I, D>> stabOpen(T point) {
		List<Entry<T, I, D>> results = new LinkedList<Entry<T, I, D>>();
		root.stabOpen(point, results);
		return results;
	}

	/**
	 * Performs a stabbing query, treating intervals as closed
	 * 
	 * @param point query point
	 * @return list of hits as entries
	 */
	public List<Entry<T, I, D>> stabClosed(T point) {
		List<Entry<T, I, D>> results = new LinkedList<Entry<T, I, D>>();
		root.stabClosed(point, results);
		return results;
	}

	/**
	 * Performs an interval intersect query, treating all intervals as open.
	 * 
	 * @param interval query interval
	 * @return list of hits as entries
	 */
	public List<Entry<T, I, D>> queryOverlaps(Interval<T> query) {
		List<Entry<T, I, D>> results = new LinkedList<Entry<T, I, D>>();
		root.queryOverlaps(query, results);
		return results;
	}

	/**
	 * Performs an interval intersect query, treating all intervals as open.
	 * 
	 * @param interval query interval
	 * @return list of hits as entries
	 */
	public List<Entry<T, I, D>> queryOverlapsOrAbuts(Interval<T> query) {
		List<Entry<T, I, D>> results = new LinkedList<Entry<T, I, D>>();
		root.queryOverlapsOrAbuts(query, results);
		return results;
	}
	
	/**
	 * Returns the total number of entries in the tree.
	 */
	public int size() {
		return root.size();
	}
}
