package org.stjude.compbio.common.gene;

import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.RefStranded;
import org.stjude.compbio.common.interval.RefInterval;

/**
 * A single annotation of a named transcript somewhere in the genome.
 * 
 * The same named transcript could be annotated in other locations, which would
 * require more TxHit objects.
 */
public class TxHit implements RefStranded {
	/**
	 * The group to which this object belongs
	 */
	private TxHitGroup txHitGroup;
	/**
	 * The TxModel
	 */
	private TxModel txModel;

	/**
	 * Creates a new TxHit instance, and adds the new instance to the
	 * TxHitGroup.
	 * 
	 * @param txHitGroup the tx group to which this tx belongs
	 * @param txModel tx model
	 */
	static TxHit newInstance(TxHitGroup txHitGroup, TxModel txModel) {
		TxHit txHit = new TxHit(txHitGroup, txModel);
		txHitGroup.addTxHit(txHit);
		return txHit;
	}

	/**
	 * Constructor using all fields
	 * 
	 * @param txHitGroup the group to which this object belongs
	 * @param txModel the model
	 */
	private TxHit(TxHitGroup txHitGroup, TxModel txModel) {
		this.txHitGroup = txHitGroup;
		this.txModel = txModel;
	}
	
	public TxHitGroup getTxHitGroup() {
		return txHitGroup;
	}

	public TxModel getTxModel() {
		return txModel;
	}

	public String getTxName() {
		return txHitGroup.getName();
	}
	
	public String getGeneName() {
		return txHitGroup.getGeneName();
	}

	public String getRefName() {
		return txModel.getRefName();
	}

	public Strand getStrand() {
		return txModel.getStrand();
	}

	public RefInterval getTxInterval() {
		return txModel.getTxInterval();
	}
	
	/**
	 * Returns true if the TxHits' models are equal
	 */
	public boolean hasEqualTxModel(TxHit other) {
		return txModel.equals(other.txModel);
	}
}
