package org.stjude.compbio.common;

/**
 * Strand/orientation relative to a reference
 */
public enum Strand {
	NA(0, "N/A", "N/A", false),
	FORWARD(1, "+", "fwd", false),
	REVERSE(-1, "-", "rev", true);
	
	/**
	 * Minimum strand when ordering
	 */
	public static final Strand MIN_STRAND = NA;
	/**
	 * Maximum strand when ordering
	 */
	public static final Strand MAX_STRAND = REVERSE;
	
	private final int orient;
	private final String str;
	private final String threeCharStr;
	private final boolean reverse;
	
	private Strand(
			int orient, String str, String threeCharStr, boolean reverse) {
		this.orient = orient;
		this.str = str;
		this.threeCharStr = threeCharStr;
		this.reverse = reverse;
	}

	/**
	 * Get integer orientation representation of the strand (-1, 0, or 1)
	 */
	public int getOrient() {
		return orient;
	}

	/**
	 * Returns true for the reverse strand.
	 * 
	 * @return true for reverse, false for forward or N/A
	 */
	public boolean isReverse() {
		return reverse;
	}

	/**
	 * Returns the strand that is the reverse of this strand
	 */
	public Strand flip() {
		return valueOfOrient(-1 * orient);
	}
	
	/**
	 * Gets RefStranded for an integer orientation
	 * @param orient integer orientation
	 * @return FORWARD for positive, REVERSE for negative, null for 0
	 */
	public static Strand valueOfOrient(int orient) {
		if(orient > 0) return FORWARD;
		if(orient < 0) return REVERSE;
		return NA;
	}
	
	/**
	 * Gets RefStranded for +/- representation
	 *
	 * @param str +/- strand representation
	 * @return FORWARD or REVERSE
	 * @throws IllegalArgumentException if chr is not + or -
	 */
	public static Strand valueOfStr(String str) {
		for(Strand strand: values()) {
			if(strand.str.equals(str)) return strand;
		}
		throw new IllegalArgumentException(
				"Invalid strand representation: " + str);
	}
	
	/**
	 * Converts to a character ('+', '-', or null)
	 */
	public Character toCharacter() {
		return (str.length() == 1) ? str.charAt(0) : null;
	}

	/**
	 * Gets 3-character string representation of fwd, rev, or N/A
	 */
	public String toThreeCharString() {
		return threeCharStr;
	}

	/**
	 * Gets string representation of +, -, or N/A
	 */
	@Override
	public String toString() {
		return str;
	}
}
