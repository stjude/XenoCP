package org.stjude.compbio.common.interval;

/**
 * Abstract base implementation of LongInterval
 * 
 * @param <T> endpoint type
 * @param <I> interval type returned by interval-returning methods
 */
public abstract class AbstractLongInterval<I extends AbstractLongInterval<I>>
extends AbstractInterval<Long, I> 
implements LongInterval {
	/**
	 * Constructs a StdLongInterval using endpoints
	 * 
	 * @param left left endpoint
	 * @param right right endpoint
	 * @throws NullPointerException if either is null
	 */
	protected AbstractLongInterval(Long left, Long right) {
		super(left, right);
	}

	/**
	 * Copy constructor
	 */
	protected AbstractLongInterval(Interval<? extends Long> src) {
		super(src);
	}

	public int getLeftInt() {
		Long endpoint = getLeft();
		if(endpoint > Integer.MAX_VALUE) {
			throw new IllegalStateException(
					"Left endpoint is too big for an int: " + endpoint);
		}
		return endpoint.intValue();
	}

	public int getRightInt() {
		Long endpoint = getRight();
		if(endpoint > Integer.MAX_VALUE) {
			throw new IllegalStateException(
					"Right endpoint is too big for an int: " + endpoint);
		}
		return endpoint.intValue();
	}

	public long length() {
		return getRight() - getLeft();
	}

	public long getMidpoint() {
		return (getLeft() + getRight()) / 2;
	}

	public long overlapLength(Interval<? extends Long> other) {
		return Math.min(getRight(), other.getRight()) -
				Math.max(getLeft(), other.getLeft());
	}

	public I shift(long amount) {
		I out = copy();
		out.setShift(amount);
		return out;
	}
	public I expandLeft(long amount) {
		long newLeft = getLeft() - amount;
		if(newLeft > getRight()) throw new IllegalArgumentException(
				"Cannot expand left by " + amount);
		
		I out = copy();
		out.setLeft(newLeft);
		return out;
	}
	public I expandRight(long amount) {
		long newRight = getRight() + amount;
		if(getLeft() > newRight) throw new IllegalArgumentException(
				"Cannot expand right by " + amount);
		
		I out = copy();
		out.setRight(newRight);
		return out;
	}
	
	/**
	 * Sets left and right by shifting by an amount.
	 */
	protected void setShift(long amount) {
		setLeft(getLeft() + amount);
		setRight(getRight() + amount);
	}
}
