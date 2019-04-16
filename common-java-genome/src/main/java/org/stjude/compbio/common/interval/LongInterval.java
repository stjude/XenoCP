package org.stjude.compbio.common.interval;

public interface LongInterval extends Interval<Long> {
	/**
	 * Gets the left endpoint as an int
	 */
	int getLeftInt();
	
	/**
	 * Gets the right endpoint as an int
	 */
	int getRightInt();
	
	/**
	 * Returns the interval length
	 */
	long length();
	
	/**
	 * Gets the midpoint
	 */
	long getMidpoint();
	
	/**
	 * Returns length of overlapping interval.
	 * 
	 * Return value is <= 0, but otherwise undefined for non-overlapping
	 * intervals.
	 */
	long overlapLength(Interval<? extends Long> other);
	
	/**
	 * Returns a new interval which is the same as this but with endpoints
	 * shifted
	 */
	LongInterval shift(long amount);
	/**
	 * Returns a new interval which is the same as this but with the left
	 * endpoint shifted
	 * 
	 * @param amount amount to shift OUTWARD (negative to shrink)
	 * @throws IllegalArgumentException if left would become greater than right 
	 */
	LongInterval expandLeft(long amount);
	/**
	 * Returns a new interval which is the same as this but with the right
	 * endpoint shifted
	 * 
	 * @param amount amount to shift OUTWARD (negative to shrink)
	 * @throws IllegalArgumentException if left would become greater than right 
	 */
	LongInterval expandRight(long amount);

	/**
	 * Returns the intersection of this interval with the given interval.
	 * 
	 * @throws IllegalArgumentException if intervals don't overlap
	 */
	LongInterval intersect(Interval<? extends Long> other)
	throws IllegalArgumentException;
	/**
	 * Returns the union of this interval with the given interval.
	 * 
	 * This also includes everything in between, if the intervals do not overlap.
	 */
	LongInterval union(Interval<? extends Long> other);
}
