package org.stjude.compbio.common.interval;

public abstract class AbstractRefInterval<I extends AbstractRefInterval<I>>
extends AbstractLongInterval<I>
implements RefInterval {
	protected AbstractRefInterval(long left, long right) {
		super(left, right);
	}

	protected AbstractRefInterval(RefInterval src) {
		this(src.getLeft(), src.getRight());
	}
	
	public long getLeftInBase() {
		return getLeft() + 1;
	}

	public int getLeftInBaseInt() {
		return getLeftInt() + 1;
	}

	public long getStart() {
		return getStrand().isReverse() ? getRight() : getLeft();
	}

	public long getStartInBase() {
		return getStrand().isReverse() ? getRight() : getLeftInBase();
	}

	public long getEnd() {
		return getStrand().isReverse() ? getLeft() : getRight();
	}

	public long getEndInBase() {
		return getStrand().isReverse() ? getLeftInBase() : getRight();
	}

	public boolean containsInBase(long point) {
		return getLeft() < point && point <= getRight();
	}

	public I getStartInterval(Long point) {
		I out = copy();
		out.setEnd(point);
		return out;
	}

	protected void setEnd(Long point) {
		if(getStrand().isReverse()) setLeft(point);
		else setRight(point);
	}

	public I getEndInterval(Long point) {
		I out = copy();
		out.setStart(point);
		return out;
	}

	protected void setStart(Long point) {
		if(getStrand().isReverse()) setRight(point);
		else setLeft(point);
	}

	public I expandStart(long amount) {
		if(getStrand().isReverse()) return expandRight(amount);
		else return expandLeft(amount);
	}

	public I expandEnd(long amount) {
		if(getStrand().isReverse()) return expandLeft(amount);
		else return expandRight(amount);
	}

	@Override
	public int hashCode() {
		int orient = (getStrand() == null) ? 3 : getStrand().getOrient();
		return 37 * super.hashCode() + orient * getRefName().hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if(obj == this) return true;
		if(obj == null) return false;
		if(!(obj instanceof RefInterval)) return false;
		
		// Test on the basis of endpoints
		if(!super.equals(obj)) return false;

		RefInterval other = (RefInterval) obj;
		// Test on the basis of ref name
		if(getRefName() == null) {
			if(other.getRefName() != null) return false;
		}
		else if(!getRefName().equals(other.getRefName())) {
			return false;
		}
		// Test on the basis of strand
		return getStrand() == other.getStrand();
	}

	public int compareTo(RefInterval o) {
		int cmp = getRefName().compareTo(o.getRefName());
		if(cmp != 0) return cmp;
		
		cmp = getStrand().compareTo(o.getStrand());
		if(cmp != 0) return cmp;
		
		return super.compareTo(o);
	}

	@Override
	public String toString() {
		return getRefName() + ":" + super.toString() + "(" +
				getStrand() + ")";
	}
	
	
}
