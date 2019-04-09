package org.stjude.compbio.sam.block;

import htsjdk.samtools.CigarOperator;

/**
 * Whether the block covers read bases, reference bases, both, or neither
 */
public enum BlockForm { 
	READ_ONLY(true, false),
	REF_ONLY(false, true),
	BOTH(true, true),
	NEITHER(false, false),
	REF_SKIP(false, false, true),
	;

	public final boolean coversRead;
	public final boolean advancesRead;
	public final boolean coversRef;
	public final boolean advancesRef;
	private BlockForm(boolean coversRead, boolean coversRef) {
		this(coversRead, coversRef, coversRef);
	}
	private BlockForm(boolean coversRead, boolean coversRef,
			boolean advancesRef) {
		this.coversRead = coversRead;
		this.advancesRead = coversRead;
		this.coversRef = coversRef;
		this.advancesRef = advancesRef;
	}
	
	/**
	 * Gets the Covers instance based on booleans for read/ref coverage
	 * 
	 * @param coversRead whether or not read bases are covered
	 * @param coversRef whetehr or not ref bases are covered
	 * @return Covers value
	 */
	public final static BlockForm valueOf(
			boolean coversRead, boolean coversRef) {
		if(coversRead) {
			return coversRef ? BOTH : READ_ONLY;
		}
		else {
			return coversRef ? REF_ONLY : NEITHER;
		}
	}
	
	/**
	 * Gets the Covers instance based on a CigarOperator
	 * 
	 * @param op CigarOperator whose coverage to describe
	 * @return Covers value
	 */
	public final static BlockForm valueOf(CigarOperator op) {
		// Handle N as a special case
		if(op == CigarOperator.N) return BlockForm.REF_SKIP;
		// Then, base it on op covers info
		return valueOf(op.consumesReadBases(), op.consumesReferenceBases());
	}
}