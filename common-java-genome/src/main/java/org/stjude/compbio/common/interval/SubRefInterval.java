package org.stjude.compbio.common.interval;

import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.RefStranded;

/**
 * RefInterval which is part of another RefStranded object.
 * 
 * It derives it's ref name and strand from the parent object.
 */
public class SubRefInterval
extends AbstractRefInterval<SubRefInterval> {
	/**
	 * Parent RefStranded
	 */
	private RefStranded parent;
	/**
	 * Ordinal, or null if not applicable
	 */
	private Integer ordinal;

	/**
	 * Constructs a SubRefInterval using parent, ordinal, and reference
	 * coordinates
	 * 
	 * @param parent containing RefStranded object
	 * @param ordinal ordinal with respect to the parent, or null if N/A
	 * @param left left (lower) endpoint
	 * @param right right (upper) endpoint
	 */
	public SubRefInterval(RefStranded parent, Integer ordinal, long left, long right) {
		super(left, right);
		this.parent = parent;
		this.ordinal = ordinal;
	}

	/**
	 * Constructs a SubRefInterval without ordinal
	 * 
	 * @param parent containing RefStranded object
	 * @param left left (lower) endpoint
	 * @param right right (upper) endpoint
	 */
	public SubRefInterval(RefStranded parent, long left, long right) {
		this(parent, null, left, right);
	}

	/**
	 * Copy constructor
	 */
	public SubRefInterval(SubRefInterval src) {
		this(src.getParent(), src.getOrdinal(), src.getLeft(), src.getRight());
	}

	public SubRefInterval(RefStranded parent, Integer ordinal,
			Interval<? extends Long> src) {
		this(parent, ordinal, src.getLeft(), src.getRight());
	}

	@Override
	protected SubRefInterval copy() {
		return new SubRefInterval(this);
	}

	/**
	 * Gets the parent RefStranded object
	 */
	public RefStranded getParent() {
		return parent;
	}

	/**
	 * Gets the ordinal with respect to the parent
	 * 
	 * @return ordinal, or null if not applicable
	 */
	public Integer getOrdinal() {
		return ordinal;
	}

	public String getRefName() {
		return parent.getRefName();
	}

	public Strand getStrand() {
		return parent.getStrand();
	}

}
