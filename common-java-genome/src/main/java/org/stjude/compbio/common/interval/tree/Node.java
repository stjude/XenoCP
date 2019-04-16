package org.stjude.compbio.common.interval.tree;

import java.util.LinkedList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.stjude.compbio.common.interval.GenericInterval;
import org.stjude.compbio.common.interval.Interval;

/**
 * IntervalTree node
 *
 * @param <T> coordinate type
 * @param <I> interval type
 * @param <D> data type
 */
class Node<T extends Comparable<T>, I extends Interval<T>, D> {
	/**
	 * Node "center" point (if null, then left and right must be null)
	 */
	private T center;
	/**
	 * Entries that overlap the center
	 */
	private List<Entry<T, I, D>> entries;
	/**
	 * Left node (may be null)
	 */
	private Node<T, I, D> left;
	/**
	 * Right node (may be null)
	 */
	private Node<T, I, D> right;
	/**
	 * Total number of entries in this node and descendant nodes
	 */
	private int size;
	
	/**
	 * Method for determining center point from a list of entries
	 * 
	 * @param entries entries to use to determine center
	 * @return center point
	 */
	private static <T extends Comparable<T>, I extends Interval<T>, D> T
	determineCenter(List<Entry<T, I, D>> entries) {
		if(entries.isEmpty()) return null;
		return entries.get(entries.size() / 2).getInterval().getLeft();
	}
	
	public Node(List<Entry<T, I, D>> entries) {
		this(entries, determineCenter(entries));
	}

	private Node(List<Entry<T, I, D>> entries, T center) {
		this(center);
		this.build(entries);
	}
	
	private Node(T center) {
		this.center = center;
		this.entries = new LinkedList<Entry<T, I, D>>();
		this.left = null;
		this.right = null;
		this.size = 0;
	}
	
	/**
	 * Builder for child nodes
	 */
	private class ChildNodeBuilder {
		private List<Entry<T, I, D>> entryList =
				new LinkedList<Entry<T, I, D>>();
		private SortedSet<T> points =
				new TreeSet<T>();
		/**
		 * Adds an entry which will need to go in the child node
		 */
		public void addEntry(Entry<T, I, D> entry) {
			entryList.add(entry);
			points.add(entry.getInterval().getLeft());
			points.add(entry.getInterval().getRight());
		}
		/**
		 * Picks the center as the middle endpoint
		 */
		private T getCenter() {
			int index = 0;
			int target = points.size() / 2;
			T newCenter = null;
			for(T point: points) {
				if(++index >= target) {
					newCenter = point;
					break;
				}
			}
			return newCenter;
		}
		/**
		 * Builds and returns the child node
		 */
		public Node<T, I, D> toNode() {
			if(entryList.isEmpty()) return null;
			return new Node<T, I, D>(entryList, getCenter());
		}
	}
	
	/**
	 * Builds out the tree under this node, with a predetermined center.
	 *  
	 * @param entries the entries to add here and below
	 * @param center center to assign
	 */
	private void build(List<Entry<T, I, D>> entries) {
		// Otherwise, fill out left/right lists
		ChildNodeBuilder leftBuilder = new ChildNodeBuilder();
		ChildNodeBuilder rightBuilder = new ChildNodeBuilder();
		for(Entry<T, I, D> entry: entries) {
			I interval = entry.getInterval();
			if(isLeft(interval)) {
				leftBuilder.addEntry(entry);
			}
			else if(isRight(interval)) {
				rightBuilder.addEntry(entry);
			}
			else {
				this.entries.add(entry);
			}
		}
		
		// Recursively build child nodes
		this.left = leftBuilder.toNode();
		this.right = rightBuilder.toNode();
		
		// Set size
		this.size = entries.size();
	}

	/**
	 * Returns true if interval is wholly to the left of center
	 */
	private boolean isLeft(Interval<T> interval) {
		if(center == null) return false;
		return interval.getRight().compareTo(center) < 0;
	}

	/**
	 * Returns true if interval is wholly to the right of center
	 */
	private boolean isRight(Interval<T> interval) {
		if(center == null) return false;
		return interval.getLeft().compareTo(center) > 0;
	}

	/**
	 * Returns the child node for the interval, or null if it crosses center
	 */
	private Node<T, I, D> getChild(Interval<T> interval) {
		if(isLeft(interval)) {
			return left;
		}
		else if(isRight(interval)) {
			return right;
		}
		else {
			return null;
		}
	}
	
	/**
	 * Adds an entry to this node or a descendant.
	 */
	public void add(Entry<T, I, D> entry) {
		// Get the interval and assert it isn't backwards (backwards intervals
		// cause problems!)
		I interval = entry.getInterval();
		if(interval.getLeft().compareTo(interval.getRight()) > 0) { 
			throw new IllegalArgumentException(
					"Cannot add backwards interval to tree: " + interval);
		}
		
		// Now, add to this node or child
		Node<T, I, D> node = getChild(interval);
		if(node == null) {
			entries.add(entry);
		}
		else {
			node.add(entry);
		}
		++size;
	}
	
	/**
	 * Removes an entry from this node or a descendant.
	 * @return true if entry was found/removed
	 */
	public boolean remove(Entry<T, I, D> entry) {
		Node<T, I, D> node = getChild(entry.getInterval());
		boolean removed;
		if(node == null) {
			removed = entries.remove(entry);
		}
		else {
			removed = node.remove(entry);
		}
		if(removed) --size;
		return removed;
	}
	
	/**
	 * Performs a stabbing query, treating intervals as open
	 * 
	 * @param point query point
	 * @param results list of entries to be appended to
	 */
	public void stabOpen(T point, List<Entry<T, I, D>> results) {
		queryOverlaps(new GenericInterval<T>(point, point), results);
	}
	
	/**
	 * Performs a stabbing query, treating intervals as closed
	 * 
	 * @param point query point
	 * @param results list of entries to be appended to
	 */
	public void stabClosed(T point, List<Entry<T, I, D>> results) {
		queryOverlapsOrAbuts(new GenericInterval<T>(point, point), results);
	}
	
	/**
	 * Performs an interval intersect query, treating all intervals as open.
	 * 
	 * @param interval query interval
	 * @param results list of entries to be appended to
	 */
	public void queryOverlaps(Interval<T> query,
			List<Entry<T, I, D>> results) {
		queryOverlaps(query, results, false);
	}
	
	/**
	 * Performs an interval intersect query, treating all intervals as closed,
	 * and including intervals that overlap only by a single endpoint.
	 * 
	 * @param interval query interval
	 * @param results list of entries to be appended to
	 */
	public void queryOverlapsOrAbuts(Interval<T> query,
			List<Entry<T, I, D>> results) {
		queryOverlaps(query, results, true);
	}
	
	/**
	 * Performs an interval intersect query, with abut option.
	 * 
	 * @param interval query interval
	 * @param results list of entries to be appended to
	 * @param abuts whether or not to include intervals that abut the query
	 */
	private void queryOverlaps(Interval<T> query,
			List<Entry<T, I, D>> results, boolean abuts) {
		// Add left unless right only
		if(!isRight(query) && left != null) {
			left.queryOverlaps(query, results, abuts);
		}
		
		// Add center
		for(Entry<T, I, D> entry: entries) {
			if(abuts
					? entry.getInterval().overlapsOrAbuts(query)
					: entry.getInterval().overlaps(query)) {
				results.add(entry);
			}
		}
		
		// Add right unless left only
		if(!isLeft(query) && right != null) {
			right.queryOverlaps(query, results, abuts);
		}
	}
	
	/**
	 * Adds all entries in or under this node to the result list
	 * 
	 * @param results list of entries to be appended to
	 */
	public void getAllEntries(List<Entry<T, I, D>> results) {
		if(left != null) left.getAllEntries(results);
		results.addAll(entries);
		if(right != null) right.getAllEntries(results);
	}
	
	/**
	 * Returns the total number of entries in or under this node.
	 */
	public int size() {
		return size;
	}
}
