package org.stjude.compbio.common.interval;

public class StdLongInterval extends AbstractLongInterval<StdLongInterval> {
	/**
	 * Constructs a StdLongInterval using endpoints
	 * 
	 * @param left left endpoint
	 * @param right right endpoint
	 * @throws NullPointerException if either is null
	 */
	public StdLongInterval(Long left, Long right) {
		super(left, right);
	}

	/**
	 * Copy constructor
	 */
	public StdLongInterval(Interval<? extends Long> src) {
		super(src);
	}

	public StdLongInterval(int left, int right) {
		this((long)left, (long)right);
	}

	@Override
	protected StdLongInterval copy() {
		return new StdLongInterval(this);
	}
	
	/**
	 * Constructs a StdLongInterval which represents a single in-base position
	 * 
	 * @param inBasePos in-base position
	 * @return interval of length 1
	 */
	public StdLongInterval inBaseBaseInterval(long inBasePos) {
		return new StdLongInterval(inBasePos - 1, inBasePos);
	}
}
