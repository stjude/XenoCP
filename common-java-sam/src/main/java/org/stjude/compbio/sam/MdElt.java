package org.stjude.compbio.sam;

import java.util.Arrays;
import java.util.regex.Pattern;

/**
 * Data structure representing one element in an MD tag
 */
public class MdElt {
	/**
	 * Element type: match, mismatch, or deletion
	 */
	public static enum Type {
		MATCH("^[0-9]*"),
		MISMATCH("^[A-Za-z]*"),
		DELETION("^\\^[A-Za-z]*");
		
		public final Pattern pattern;
		
		private Type(String patternStr) {
			this.pattern = Pattern.compile(patternStr);
		}
	}
	
	/**
	 * The type: match, mismatch, or deletion
	 */
	public final Type type;
	/**
	 * Number of reference bases covered
	 */
	public final int length;
	/**
	 * The reference bases for MISMATCH or DELETION, or null for MATCH.
	 */
	public final byte[] bases;
	
	public MdElt(Type type, int length, byte[] bases) {
		this.type = type;
		this.length = length;
		this.bases = bases;
	}
	
	/**
	 * Returns a new sub-MdElt starting at position start in this element
	 * @param start starting position (length to skip)
	 * @return new MdElt
	 */
	public MdElt subMdElt(int start) {
		return new MdElt(
				type,
				length - start,
				bases == null ? null :
					Arrays.copyOfRange(bases, start, bases.length));
	}

	/**
	 * Returns the canonical String representation of the MdElt, as would appear
	 * in an MD string in a SAMRecord.
	 */
	@Override
	public String toString() {
		switch(type) {
		case MATCH:
			return Integer.toString(length);
		case MISMATCH:
			return new String(bases);
		case DELETION:
			return "^" + new String(bases);
		default:
			throw new RuntimeException("Impossible type"); 
		}
	}
}
