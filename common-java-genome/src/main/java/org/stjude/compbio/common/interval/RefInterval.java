package org.stjude.compbio.common.interval;

import org.stjude.compbio.common.RefStranded;

/**
 * An interval on a strand of a reference.
 * 
 * Two RefIntervals should be considered equal only if their ref names, strands,
 * and intervals are equal.
 */
public interface RefInterval extends LongInterval, RefStranded {
	/**
	 * Gets left position in in-base coordinates
	 */
	long getLeftInBase();
	/**
	 * Gets left position as int in in-base coordinates
	 */
	int getLeftInBaseInt();
	
	/**
	 * Gets start position (using orientation)
	 */
	long getStart();
	/**
	 * Gets start position (using orientation) in in-base coordinates
	 */
	long getStartInBase();
	
	/**
	 * Gets end position (using orientation)
	 */
	long getEnd();
	/**
	 * Gets end position (using orientation) in in-base coordinates
	 */
	long getEndInBase();
	
	/**
	 * Returns true if the in-base pos is contained in the interval.
	 * 
	 * @param point point to test, in in-base coordinates
	 */
	boolean containsInBase(long point);

	/**
	 * Gets the interval between the start and the given point
	 */
	RefInterval getStartInterval(Long point);
	/**
	 * Gets the interval between the given point and the end
	 */
	RefInterval getEndInterval(Long point);
	
	/**
	 * Returns a new interval which is the same as this but with endpoints
	 * shifted
	 */
	RefInterval shift(long amount);
	/**
	 * Returns a new interval which is the same as this but with the start
	 * endpoint shifted
	 * 
	 * @param amount amount to shift OUTWARD (negative to shrink)
	 * @throws IllegalArgumentException if left would become greater than right 
	 */
	RefInterval expandStart(long amount);
	/**
	 * Returns a new interval which is the same as this but with the end
	 * endpoint shifted
	 * 
	 * @param amount amount to shift OUTWARD (negative to shrink)
	 * @throws IllegalArgumentException if left would become greater than right 
	 */
	RefInterval expandEnd(long amount);
	
	/**
	 * Returns the intersection of this interval with the given interval.
	 * 
	 * @throws IllegalArgumentException if intervals don't overlap
	 */
	RefInterval intersect(Interval<? extends Long> other)
	throws IllegalArgumentException;
	/**
	 * Returns the union of this interval with the given interval.
	 * 
	 * This also includes everything in between, if the intervals do not overlap.
	 */
	RefInterval union(Interval<? extends Long> other);
}
