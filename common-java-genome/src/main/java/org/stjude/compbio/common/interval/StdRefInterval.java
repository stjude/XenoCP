package org.stjude.compbio.common.interval;

import org.stjude.compbio.common.RefStrand;
import org.stjude.compbio.common.Strand;

public class StdRefInterval
extends AbstractRefInterval<StdRefInterval> {
	/**
	 * Reference name/strand
	 */
	private RefStrand refStrand;

	/**
	 * Constructs a StdRefInterval using location information
	 * 
	 * @param refName reference name
	 * @param strand strand
	 * @param left left (lower) endpoint
	 * @param right right (upper) endpoint
	 */
	public StdRefInterval(String refName, Strand strand, long left,
			long right) {
		super(left, right);
		this.refStrand = new RefStrand(refName, strand);
	}

	/**
	 * Copy constructor
	 */
	public StdRefInterval(RefInterval src) {
		this(src.getRefName(), src.getStrand(),
				src.getLeft(), src.getRight());
	}

	/**
	 * Constructs an unstranded (strand == NA) ref interval
	 *  
	 * @param refName reference name
	 * @param left left (lower) endpoint
	 * @param right right (upper) endpoint
	 */
	public StdRefInterval(String refName, long left, long right) {
		this(refName, Strand.NA, left, right);
	}

	@Override
	protected StdRefInterval copy() {
		return new StdRefInterval(this);
	}

	public String getRefName() {
		return refStrand.getRefName();
	}

	public Strand getStrand() {
		return refStrand.getStrand();
	}

}
