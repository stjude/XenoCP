package org.stjude.compbio.common.interval;

/**
 * A RefInterval that additionally has a name
 */
public interface NamedRefInterval extends RefInterval {
	/**
	 * Returns the interval's name
	 */
	String getName();

	/**
	 * Gets the interval between the start and the given point
	 */
	NamedRefInterval getStartInterval(Long point);
	/**
	 * Gets the interval between the given point and the end
	 */
	NamedRefInterval getEndInterval(Long point);
	
	/**
	 * Returns a new interval which is the same as this but with endpoints
	 * shifted
	 */
	NamedRefInterval shift(long amount);
	/**
	 * Returns a new interval which is the same as this but with the start
	 * endpoint shifted
	 * 
	 * @param amount amount to shift OUTWARD (negative to shrink)
	 * @throws IllegalArgumentException if left would become greater than right 
	 */
	NamedRefInterval expandStart(long amount);
	/**
	 * Returns a new interval which is the same as this but with the end
	 * endpoint shifted
	 * 
	 * @param amount amount to shift OUTWARD (negative to shrink)
	 * @throws IllegalArgumentException if left would become greater than right 
	 */
	NamedRefInterval expandEnd(long amount);
	
	/**
	 * Returns the intersection of this interval with the given interval.
	 * 
	 * @throws IllegalArgumentException if intervals don't overlap
	 */
	NamedRefInterval intersect(Interval<? extends Long> other)
	throws IllegalArgumentException;
	/**
	 * Returns the union of this interval with the given interval.
	 * 
	 * This also includes everything in between, if the intervals do not overlap.
	 */
	NamedRefInterval union(Interval<? extends Long> other);
}
