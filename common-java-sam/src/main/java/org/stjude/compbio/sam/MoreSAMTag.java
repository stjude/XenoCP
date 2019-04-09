package org.stjude.compbio.sam;

/**
 * SAM tags that are not represented in Picard's SAMTag enum
 * 
 * Since the set of standard SAM tags changes, it is recommended to use the
 * static constants, rather than the true enum values.
 */
public enum MoreSAMTag {
  YR, OP, OC, ZC, ZM, XU, XS, XX;
  
  /**
   * Reference name in original mapping, before transformation
   */
  public static final MoreSAMTag ORIG_REF = YR;
  /**
   * Position in original mapping, before transformation
   */
  public static final MoreSAMTag ORIG_POS = OP;
  /**
   * CIGAR in original mapping, before transformation
   */
  public static final MoreSAMTag ORIG_CIGAR = OC;
  /**
   * Mate's CIGAR
   */
  public static final MoreSAMTag MATE_CIGAR = ZC;
  /**
   * Mate's MD
   */
  public static final MoreSAMTag MATE_MD = ZM;
  /**
   * Declared to be interesting by Bambino's SAMExtractUnmapped
   */
  public static final MoreSAMTag SEU_INTERESTING = XU;
  /**
   * Transcribed strand: + or -
   */
  public static final MoreSAMTag TX_STRAND = XS;
  /**
   * Extra sequences, as a comma separated list of label=SEQ/QUAL,...
   */
  public static final MoreSAMTag EXTRA_SEQS = XX;
}
