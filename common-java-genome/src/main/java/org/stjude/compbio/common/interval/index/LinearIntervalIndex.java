package org.stjude.compbio.common.interval.index;

import java.util.List;
import java.util.Set;

import org.stjude.compbio.common.interval.Interval;
import org.stjude.compbio.common.interval.tree.Entry;

/**
 * A data structure holding data objects tied to intervals, with the ability to
 * efficiently query based on points or intervals.
 * 
 * Most IntervalIndex implementations will probably use an IntervalTree under
 * the hood, but IntervalIndex is designed to be more user-friendly.
 * 
 * @param <T> coordinate type
 * @param <I> interval type
 * @param <D> data type
 */
public interface LinearIntervalIndex
<T extends Comparable<T>, I extends Interval<T>, D>
extends IntervalIndex<T, I, D> {
	/**
	 * Retrieve all data objects at intervals that intersect the point.
	 * 
	 * Intervals are considered to intersect their endpoints.
	 * 
	 * If there is more than one hit entry with equal data objects, then both
	 * equal data objects are returned in the list.
	 * 
	 * @param point query point
	 * @return data points whose intervals intersect the point
	 */
	List<D> query(T point);

	/**
	 * Retrieve all distinct data objects at intervals that intersect the point.
	 * 
	 * Intervals are considered to intersect their endpoints.
	 * 
	 * If there is more than one hit entry with equal data objects, then only
	 * one data object is returned in the set
	 * 
	 * @param point query point
	 * @return set of distinct data points whose intervals intersect the point
	 */
	Set<D> queryDistinct(T point);

	
	
	// Entries versions
	
	/**
	 * Retrieve all entries at intervals that intersect the point.
	 * 
	 * Intervals are considered to intersect their endpoints.
	 * 
	 * @param point query point
	 * @return entries whose intervals intersect the point
	 */
	List<Entry<T, I, D>> queryEntries(T point);
	
	/**
	 * Returns an interval that contains all other intervals in the index
	 */
	Interval<T> getScope();
}
