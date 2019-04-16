package org.stjude.compbio.common.interval;

/**
 * Abstract base implementation of interval
 * 
 * @param <T> endpoint type
 * @param <I> interval type returned by interval-returning methods
 */
public abstract class AbstractInterval<
	T extends Comparable<T>, I extends AbstractInterval<T, I>>
implements Interval<T> {
	private T left;
	private T right;
	
	/**
	 * Construct a generic interval using endpoints
	 * 
	 * @param left left (lower) endpoint
	 * @param right right (upper) endpoint
	 * @throws NullPointerException if either endpoint is null
	 */
	protected AbstractInterval(T left, T right) {
		if(left == null || right == null) {
			throw new NullPointerException(
					"Interval endpoints must be non-null; got " + left + 
					", " + right);
		}
		this.left = left;
		this.right = right;
	}
	
	/**
	 * Copy Constructor
	 */
	protected AbstractInterval(Interval<? extends T> src) {
		this.left = src.getLeft();
		this.right = src.getRight();
	}
	
	/**
	 * Abstract copy method, which must be implemented to provide support for
	 * union/intersect.
	 */
	protected abstract I copy();

	/**
	 * Sets left endpoint.
	 * 
	 * @param left left (lower) endpoint
	 * @throws NullPointerException if left is null
	 */
	protected void setLeft(T left) {
		if(left == null) throw new NullPointerException();
		this.left = left;
	}

	/**
	 * Sets right endpoint.
	 * 
	 * @param right right (upper) endpoint
	 * @throws NullPointerException if right is null
	 */
	protected void setRight(T right) {
		if(right == null) throw new NullPointerException();
		this.right = right;
	}

	public T getLeft() {
		return left;
	}

	public T getRight() {
		return right;
	}

	public boolean isEndpoint(T point) {
		return point.equals(left) || point.equals(right);
	}

	public boolean containsClosed(T point) {
		return left.compareTo(point) <= 0 &&
				point.compareTo(right) <= 0;
	}

	public boolean containsOpen(T point) {
		return left.compareTo(point) < 0 &&
				point.compareTo(right) < 0;
	}

	public boolean contains(Interval<? extends T> other) {
		return left.compareTo(other.getLeft()) <= 0 &&
				right.compareTo(other.getRight()) >= 0;
	}

	public boolean overlaps(Interval<? extends T> other) {
		return left.compareTo(other.getRight()) < 0 &&
				right.compareTo(other.getLeft()) > 0;
	}

	public boolean overlapsOrAbuts(Interval<? extends T> other) {
		return left.compareTo(other.getRight()) <= 0 &&
				right.compareTo(other.getLeft()) >= 0;
	}

	public I intersect(Interval<? extends T> other)
			throws IllegalArgumentException {
		I out = copy();
		out.setIntersect(other);
		return out;
	}
	
	/**
	 * Sets left/right by intersecting with the other interval
	 */
	protected void setIntersect(Interval<? extends T> other) {
		// Max left, min right
		T newLeft = (left.compareTo(other.getLeft()) >= 0) ?
				left : other.getLeft();
		T newRight = (right.compareTo(other.getRight()) <= 0) ?
				right : other.getRight();
		if(newLeft.compareTo(newRight) > 0) {
			throw new IllegalArgumentException(
					"Interval for intersection, " + other +
					", does not overlap this interval: " + this);
		}
		this.left = newLeft;
		this.right = newRight;
	}

	public I union(Interval<? extends T> other) {
		I out = copy();
		out.setUnion(other);
		return out;
	}
	
	/**
	 * Sets left/right by performing a union operation with the other interval
	 */
	protected void setUnion(Interval<? extends T> other) {
		// Min left, max right
		this.left = (left.compareTo(other.getLeft()) <= 0) ?
				left : other.getLeft();
		this.right = (right.compareTo(other.getRight()) >= 0) ?
				right : other.getRight();
	}

	public int compareTo(Interval<T> o) {
		int cmp = left.compareTo(o.getLeft());
		if(cmp != 0) return cmp;
		return right.compareTo(o.getRight());
	}

	@Override
	public int hashCode() {
		return left.hashCode() + 31 * right.hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if(this == obj) return true;
		if(obj == null) return false;
		if(!(obj instanceof Interval)) return false;
		
		Interval<?> other = (Interval<?>) obj;
		if(left == null) {
			if(other.getLeft() != null) return false;
		}
		else if(!left.equals(other.getLeft())) {
			return false;
		}
		if(right == null) {
			if(other.getRight() != null) return false;
		}
		else if(!right.equals(other.getRight())) {
			return false;
		}
		
		return true;
	}

	@Override
	public String toString() {
		return left + "-" + right;
	}
}
