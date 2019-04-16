package org.stjude.compbio.util;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Aggregates a list of comparators to compare two objects, and applies them in
 * order.
 *
 * This operates by using the first comparator in the list first.  If the result
 * is 0, then the second is used, and so-on.
 * 
 * @param <T> type of object being compared
 */
public class AggregateComparator<T> implements Comparator<T> {
	/**
	 * Comparators being aggregated
	 */
	private List<Comparator<T>> comparators;

	/**
	 * Constructs an empty aggregate comparator that must have comparators
	 * added.
	 */
	public AggregateComparator() {
		this(new ArrayList<Comparator<T>>());
	}
	
	/**
	 * Constructs an aggregate comparator with the given comparators.
	 */
	public AggregateComparator(List<Comparator<T>> comparators) {
		this.comparators = comparators;
	}

	public void add(Comparator<T> comparator) {
		comparators.add(comparator);
	}
	
	public int compare(T o1, T o2) {
		int cmp = 0;
		
		for(Comparator<T> comparator: comparators) {
			cmp = comparator.compare(o1, o2);
			if(cmp != 0) return cmp;
		}
		
		return cmp;
	}

}
