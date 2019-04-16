package org.stjude.compbio.common.seq;

import java.util.HashMap;
import java.util.Map;

/**
 * Simple helper methods for working with sequences
 */
public class SequenceHelper {
	/**
	 * Map from each base to its complement
	 */
	private static Map<Byte, Byte> baseComplements = new HashMap<Byte, Byte>();
	static {
		String[] ucPairs = { "AT", "CG", "NN" };
		for(String ucPair: ucPairs) {
			for(String pair: new String[] { ucPair, ucPair.toLowerCase() }) {
				byte[] bases = pair.getBytes();
				baseComplements.put(bases[0], bases[1]);
				baseComplements.put(bases[1], bases[0]);
			}
		}
	}
	
	/**
	 * Default TranslTable
	 */
	private static TranslTable defaultTranslTable = null;
	static {
		try {
			defaultTranslTable = TranslTable.getStandardTranslTable();
		} catch(Exception e) {
			// Don't do anything--this may not be a problem
		}
	}
	
	/**
	 * Returns the reverse complement of a sequence.
	 * 
	 * Supports ACGTN, and preserves case.
	 * 
	 * @param seq original sequence
	 * @return reverse complement
	 * @throws IllegalArgumentException if the sequence cannot be RC-ed
	 */
	public static String reverseComplement(String seq)
	throws IllegalArgumentException {
		return new String(reverseComplement(seq.getBytes()));
	}

	/**
	 * Returns the reverse complement of a sequence.
	 * 
	 * Supports ACGTN, and preserves case.
	 * 
	 * @param seq original sequence
	 * @return reverse complement
	 * @throws IllegalArgumentException if the sequence cannot be RC-ed
	 */
	public static byte[] reverseComplement(byte[] bases)
	throws IllegalArgumentException {
		byte[] out = bases.clone();
		try {
			reverseComplementInPlace(out);
		} catch(IllegalArgumentException e) {
			throw new IllegalArgumentException(
					"Could not reverse complement sequence {" +
							new String(bases) + "}", e);
		}
		return out;
	}

	/**
	 * Reverse complements a sequence in-place.
	 * 
	 * Supports ACGTN, and preserves case.
	 * 
	 * @param bases sequence, which is overwritten with reverse complement
	 * @throws IllegalArgumentException if the sequence cannot be RC-ed
	 */
	public static void reverseComplementInPlace(byte[] bases) 
	throws IllegalArgumentException {
		// Reverse-complement from both ends
		int i = 0;
		int j = bases.length - 1;
		try {
			while(i < j) {
				byte newAtJ = complementBase(bases[i]);
				bases[i] = complementBase(bases[j]);
				bases[j] = newAtJ;
				++i;
				--j;
			}
			// Reverse-complement the middle base, if there is one
			if(i == j) bases[i] = complementBase(bases[i]);
		} catch(IllegalArgumentException e) {
			throw new IllegalArgumentException(
					"Could not complement base(s) at positions " + i + 
					" and/or " + j, e);
		}
	}
	
	/**
	 * Complements a base, throwing a meaningful exception if it cannot be done.
	 * 
	 * @param base the base to complement
	 * @return the complement
	 * @throws IllegalArgumentException if the base cannot be complemented
	 */
	public static byte complementBase(byte base)
	throws IllegalArgumentException {
		if(baseComplements.containsKey(base)) {
			return baseComplements.get(base);
		}
		else {
			throw new IllegalArgumentException(
					"Cannot complement base {" + base + "}");
		}
	}
	
	/**
	 * Gets the default TranslTable used for translation
	 */
	public static TranslTable getDefaultTranslTable() {
		return defaultTranslTable;
	}

	/**
	 * Sets the default TranslTable used for translation
	 * @param translTable new default TranslTable
	 */
	public static void setDefaultTranslTable(TranslTable translTable) {
		SequenceHelper.defaultTranslTable = translTable;
	}
	
	/**
	 * Returns a SeqTranslator using the default transl table
	 */
	public static SeqTranslator newSeqTranslator() {
		return new SeqTranslator(defaultTranslTable);
	}
	
	/**
	 * Returns a SeqTranslator using the default transl table, with options
	 * @param strictCodon whether to require valid codons (if false uses X)
	 * @param strictFrame whether to require in-frame input
	 * @param append whether to append extra nuc bases at end of aa seq
	 */
	public static SeqTranslator newSeqTranslator(
			boolean strictCodon, boolean strictFrame, boolean append) {
		return new SeqTranslator(defaultTranslTable,
				strictCodon, strictFrame, append);
	}
}
