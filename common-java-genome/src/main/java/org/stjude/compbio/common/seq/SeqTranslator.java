package org.stjude.compbio.common.seq;

/**
 * Translate nucleic acid sequence into amino acid sequence according to a
 * TranslTable (translation table).
 * 
 * There are various options that you can set
 */
public class SeqTranslator {
	/**
	 * Translation table
	 */
	private final TranslTable tt;
	/**
	 * Whether to throw an exception if a codon does not exist in the transl
	 * table.
	 */
	private boolean strictCodon;
	/**
	 * Whether to throw an exception if the input sequence length is not a
	 * multiple of 3.
	 */
	private boolean strictFrame;
	/**
	 * Whether to append extra nucleotides in lowercase at the end of the
	 * protein sequence.
	 */
	private boolean append;
	/**
	 * Extra nucleotides from the last run
	 */
	private String extra = null;

	/**
	 * Construct a lenient sequence translator that does not append.
	 * @param tt translation table
	 */
	public SeqTranslator(TranslTable tt) {
		this(tt, false, false, false);
	}
	
	/**
	 * Construct a sequence translator with all options
	 * @param tt translation table
	 * @param strictCodon whether to require valid codons (if false uses X)
	 * @param strictFrame whether to require in-frame input
	 * @param append whether to append extra nuc bases at end of aa seq
	 */
	public SeqTranslator(TranslTable tt,
			boolean strictCodon, boolean strictFrame, boolean append) {
		if(tt == null) {
			throw new NullPointerException("Translation table is null");
		}
		this.tt = tt;
		this.strictCodon = strictCodon;
		this.strictFrame = strictFrame;
		this.append = append;
	}
	
	/**
	 * Translates nucleic acid sequence to amino acid.
	 * 
	 * If the input sequence length is not a multiple of three then an exception
	 * is thrown in strict mode, and in lenient mode the extra bases are saved
	 * to the "extra" data member, and are appended if append is turned on.
	 * 
	 * @param seq sequence to translate
	 * @return translated sequence
	 */
	public String translate(String seq) {
		// Assert in-frame if in strict mode
		if(strictFrame && seq.length() % 3 != 0) {
			throw new IllegalArgumentException(
					"Sequence length is not a multiple of 3: " + seq);
		}
		// Build output
		StringBuilder out = new StringBuilder();
		String codon = null;
		// Iterate over whole codons
		int i;
		for(i = 3; i <= seq.length(); i += 3) {
			codon = seq.substring(i - 3, i);
			char aa = 'X';
			try {
				aa = tt.translate(codon);
			}
			catch(NullPointerException e) {
				if(strictCodon) {
					throw new IllegalArgumentException(
						"Untranslatable codon {" + codon + "} in " + seq, e);
				}
			}
			out.append(aa);
		}
		// Remember extra and append if enabled
		extra = seq.substring(i - 3, seq.length());
		if(append) out.append(extra.toLowerCase());
		return out.toString();
	}
	
}
