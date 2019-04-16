package org.stjude.compbio.common.interval;

import org.stjude.compbio.common.RefStrand;
import org.stjude.compbio.common.Strand;

public class StdNamedRefInterval
extends AbstractRefInterval<StdNamedRefInterval> 
implements NamedRefInterval {
	/**
	 * Reference name/strand
	 */
	private RefStrand refStrand;
	/**
	 * Interval name
	 */
	private String name;
	
	/**
	 * Constructs a StdRefInterval using location information
	 * 
	 * @param refName reference name
	 * @param strand strand
	 * @param left left (lower) endpoint
	 * @param right right (upper) endpoint
	 * @param name interval name
	 */
	public StdNamedRefInterval(String refName, Strand strand, long left,
			long right, String name) {
		super(left, right);
		this.refStrand = new RefStrand(refName, strand);
		this.name = name;
	}

	/**
	 * Quasi-copy constructor (uses RefInterval + name, rather than 
	 * NamedRefInterval)
	 */
	public StdNamedRefInterval(RefInterval src, String name) {
		this(src.getRefName(), src.getStrand(),
				src.getLeft(), src.getRight(), name);
	}

	/**
	 * Constructs an unstranded (strand == NA) ref interval
	 *  
	 * @param refName reference name
	 * @param left left (lower) endpoint
	 * @param right right (upper) endpoint
	 * @param name interval name
	 */
	public StdNamedRefInterval(String refName, long left, long right,
			String name) {
		this(refName, Strand.NA, left, right, name);
	}

	@Override
	protected StdNamedRefInterval copy() {
		return new StdNamedRefInterval(this, this.name);
	}

	public String getRefName() {
		return refStrand.getRefName();
	}

	public Strand getStrand() {
		return refStrand.getStrand();
	}

	public String getName() {
		return this.name;
	}


	@Override
	public String toString() {
		return name + "@" + super.toString();
	}
}
