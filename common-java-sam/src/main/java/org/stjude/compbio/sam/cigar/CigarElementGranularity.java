package org.stjude.compbio.sam.cigar;

/**
 * The level of granularity to use for aligned regions.
 */
public enum CigarElementGranularity {
	/**
	 * Coarse: each continuously aligned interval gets one block, regardless
	 * of base matches and mismatches.
	 */
	COARSE,
	/**
	 * Fine: each continuously matching or mismatching interval gets one
	 * block.
	 */
	FINE,
	/**
	 * Use the level of granularity present in the source
	 */
	SOURCE
}