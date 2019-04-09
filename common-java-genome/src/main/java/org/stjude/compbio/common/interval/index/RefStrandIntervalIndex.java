package org.stjude.compbio.common.interval.index;

import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.RefStranded;
import org.stjude.compbio.common.interval.RefInterval;

/**
 * An interval index for a single strand of a reference
 * 
 * @param <D> data type
 */
public class RefStrandIntervalIndex<D>
extends GenericIntervalIndex<Long, RefInterval, D>
implements RefStranded {
	private String refName;
	private Strand refStrand;
	
	public RefStrandIntervalIndex(String refName, Strand refStrand) {
		super();
		this.refName = refName;
		this.refStrand = refStrand;
	}
	
	public RefStrandIntervalIndex(String refName, Strand refStrand,
			TreeRebuildStrategy treeRebuildStrategy) {
		super(treeRebuildStrategy);
		this.refName = refName;
		this.refStrand = refStrand;
	}

	public String getRefName() {
		return refName;
	}
	public Strand getStrand() {
		return refStrand;
	}
		
}
