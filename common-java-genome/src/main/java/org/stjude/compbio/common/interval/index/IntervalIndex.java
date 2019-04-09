package org.stjude.compbio.common.interval.index;

import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.stjude.compbio.common.interval.Interval;
import org.stjude.compbio.common.interval.tree.Entry;

/**
 * An interval search data structure, used to perform interval overlap queries
 * and return associated data objects.
 * 
 * The index is made up of entries, where each entry has an interval and a data
 * object.  You can then query the index with an interval, and it can return
 * entries or data objects with intervals that overlap the query interval.
 * 
 * This is similar to an interval tree, and many implementations will probably
 * use interval trees.  It differs in two key ways:
 * 1. It hides building logic from the user
 * 2. It does not allow stabbing queries; however, the LinearIntervalIndex adds
 *    the ability to do stabbing queries.  They are not allowed in the base
 *    interface because some intervals contain information orthogonal to the
 *    coordinates, and this information may be necessary to do lookups.
 *
 * @param <T> coordinate type
 * @param <I> interval type
 * @param <D> data type
 */
public interface IntervalIndex<T extends Comparable<T>, I extends Interval<T>, D> {

	/**
	 * Retrieve all data objects at intervals that overlap the interval.
	 * 
	 * Intervals are NOT considered to overlap intervals that merely abut them.
	 * 
	 * If there is more than one hit entry with equal data objects, then both
	 * equal data objects are returned in the list.
	 * 
	 * @param interval query interval
	 * @return data points whose intervals overlap the query interval
	 */
	List<D> query(I interval);

	/**
	 * Retrieve all distinct data objects at intervals that overlap the
	 * interval.
	 * 
	 * Intervals are NOT considered to overlap intervals that merely abut them.
	 * 
	 * If there is more than one hit entry with equal data objects, then only
	 * one data object is returned in the set
	 * 
	 * @param interval query interval
	 * @return set of distinct data points whose intervals overlap the query
	 */
	Set<D> queryDistinct(I interval);

	/**
	 * Adds an entry to the index
	 * 
	 * @param interval the interval
	 * @param data the data object at that interval
	 * @throws UnsupportedOperationException if index is unmodifiable
	 */
	void add(I interval, D data);

	/**
	 * Removes an entry from the index
	 * 
	 * @param interval the interval
	 * @param data the data object at that interval
	 * @return true if the entry was found/removed
	 * @throws UnsupportedOperationException if index is unmodifiable
	 */
	boolean remove(I interval, D data);

	/**
	 * Retrieves all entries for the given data object
	 * 
	 * This is an optional operation.
	 * 
	 * @param data the data object to query
	 * @return corresponding entries, may be empty, but never null
	 * @throws UnsupportedOperationException if not enabled
	 */
	List<Entry<T, I, D>> getDataEntries(D data);

	/**
	 * Retrieve all entries at intervals that overlap the interval.
	 * 
	 * Intervals are NOT considered to overlap intervals that merely abut them.
	 * 
	 * @param interval query interval
	 * @return entries whose intervals overlap the query interval
	 */
	List<Entry<T, I, D>> queryEntries(I interval);

	/**
	 * Adds an entry to the index
	 * 
	 * @param entry the entry to add
	 * @throws UnsupportedOperationException if index is unmodifiable
	 */
	void addEntry(Entry<T, I, D> entry);

	/**
	 * Removes an entry from the index
	 * 
	 * @param entry the entry to remove
	 * @return true if the entry was found/removed
	 * @throws UnsupportedOperationException if index is unmodifiable
	 */
	boolean removeEntry(Entry<T, I, D> entry);

	/**
	 * Removes all of the specified entries from the index
	 * 
	 * @param entries the entries to remove
	 * @return entries removed
	 * @throws UnsupportedOperationException if index is unmodifiable
	 */
	List<Entry<T, I, D>> removeAllEntries(
			Collection<Entry<T, I, D>> entries);

	/**
	 * Returns the number of entries in the index
	 */
	int size();
	
	/**
	 * Returns a collection view of all data objects in the index.
	 * 
	 * The collection may or may not be backed by the index.
	 */
	Set<D> dataSet();
}