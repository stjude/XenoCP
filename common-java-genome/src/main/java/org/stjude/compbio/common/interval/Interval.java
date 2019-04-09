package org.stjude.compbio.common.interval;

/**
 * Represents an interval between two endpoints.
 * 
 * Endpoints must not be null.  Implementations need to enforce this.
 *
 * @param <T> endpoint type
 * @param <I> type of interval returned by interval-returning methods
 */
public interface Interval<T extends Comparable<T>>
extends Comparable<Interval<T>> {
	/**
	 * Gets the left (lower) endpoint.
	 */
	T getLeft();
	/**
	 * Gets the right (upper) endpoint
	 */
	T getRight();
	
	/**
	 * Returns true if pos is an endpoint of this interval
	 * 
	 * @param point point to test
	 */
	boolean isEndpoint(T point);
	/**
	 * Returns true if pos is contained in the CLOSED interval.
	 * 
	 * The interval is treated as closed; hence if either endpoint is equal to 
	 * pos, then true is returned.
	 * 
	 * @param point point to test
	 */
	boolean containsClosed(T point);
	/**
	 * Returns true if pos is contained in the OPEN interval.
	 * 
	 * The interval is treated as open; hence if either endpoint is equal to 
	 * pos, then false is returned.
	 * 
	 * @param point point to test
	 */
	boolean containsOpen(T point);
	/**
	 * Returns true if the specified interval is contained in this interval.
	 * 
	 * @param other interval to test
	 */
	boolean contains(Interval<? extends T> other);
	/**
	 * Returns true if interval overlaps this interval.
	 * 
	 * @param other interval to test
	 */
	boolean overlaps(Interval<? extends T> other);
	/**
	 * Returns true if interval overlaps or abuts this interval.
	 * 
	 * @param other interval to test
	 */
	boolean overlapsOrAbuts(Interval<? extends T> other);

	/**
	 * Returns the intersection of this interval with the given interval.
	 * 
	 * @throws IllegalArgumentException if intervals don't overlap
	 */
	Interval<T> intersect(Interval<? extends T> other)
	throws IllegalArgumentException;
	/**
	 * Returns the union of this interval with the given interval.
	 * 
	 * This also includes everything in between, if the intervals do not overlap.
	 */
	Interval<T> union(Interval<? extends T> other);
}
