package org.stjude.compbio.common.formats.fasta;

/**
 * Represents a FASTA record.
 * 
 * The entire sequence is held in memory, so this is not necessarily suited for
 * large sequences, such as entire chromosomes.
 */
public class FastaRecord {
	/**
	 * Defline, including leading ">"
	 */
	private String defline;
	/**
	 * Sequence
	 */
	private String seq;
	
	/**
	 * Constructs a FASTA record from defline and sequence
	 * 
	 * @param defline the entire defline, including leading '>'
	 * @param seq sequence, without any line breaks
	 */
	public FastaRecord(String defline, String seq) {
		this.defline = defline;
		this.seq = seq;
	}

	/**
	 * @return the entire defline, including leading '>'
	 */
	public String getDefline() {
		return defline;
	}

	/**
	 * @return the sequence without any line breaks
	 */
	public String getSeq() {
		return seq;
	}
	
	/**
	 * Returns the first word in the defline, without leading '>'
	 */
	public String getName() {
		return defline.substring(1).split("\\s", 2)[0];
	}
}
