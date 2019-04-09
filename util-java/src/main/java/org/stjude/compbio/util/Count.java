package org.stjude.compbio.util;

/**
 * Simple immutable struct-like class to hold the result of a count of
 * some subclass of items out of a set, with a method to provide the
 * ratio.
 */
public class Count {
	/**
	 * Number of items counted
	 */
	public final int count;
	/**
	 * Total number of items considered
	 */
	public final int total;
	/**
	 * The ratio of count / total as a double
	 */
	public double ratio() {
		return divide(count, total);
	}
	public Count(int count, int total) {
		this.count = count;
		this.total = total;
	}
	/**
	 * Convenience function to divide two ints and return a double
	 * @param num numerator
	 * @param denom denominator
	 * @return quotient as a double
	 */
	public static double divide(int num, int denom) {
		return ((double)num) / denom;
	}
}