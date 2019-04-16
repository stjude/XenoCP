package org.stjude.compbio.common.interval;

/**
 * Generic implementation of interval
 * 
 * Protected setters are provided to make this a good base class for other
 * implementations.
 * 
 * @param <T> endpoint type
 */
public class GenericInterval<T extends Comparable<T>>
extends AbstractInterval<T, GenericInterval<T>> {
	/**
	 * Construct a generic interval using endpoints
	 * 
	 * @param left left (lower) endpoint
	 * @param right right (upper) endpoint
	 * @throws NullPointerException if either endpoint is null
	 */
	public GenericInterval(T left, T right) {
		super(left, right);
	}
	
	/**
	 * Copy Constructor
	 */
	public GenericInterval(Interval<? extends T> src) {
		this(src.getLeft(), src.getRight());
	}

	@Override
	protected GenericInterval<T> copy() {
		return new GenericInterval<T>(this);
	}
}
