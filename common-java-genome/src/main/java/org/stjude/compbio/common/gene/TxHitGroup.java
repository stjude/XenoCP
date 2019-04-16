package org.stjude.compbio.common.gene;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class TxHitGroup {
	/**
	 * The gene group to which this tx group belongs
	 */
	private GeneHitGroup geneHitGroup;
	/**
	 * Name, which should be unique among group instances
	 */
	private String name;
	/**
	 * TxHits, stored as a map from their models to the TxHit instances, to
	 * allow for lookups by model.
	 */
	private Map<TxModel, TxHit> txHitMap;
	
	/**
	 * Creates a new TxHitGroup instance, and adds the new instance to the
	 * GeneHitGroup.
	 * 
	 * @param geneHitGroup the gene group to which this tx group belongs
	 * @param name tx name
	 */
	static TxHitGroup newInstance(GeneHitGroup geneHitGroup, String name) {
		TxHitGroup txHitGroup = new TxHitGroup(geneHitGroup, name);
		geneHitGroup.addTxHitGroup(txHitGroup);
		return txHitGroup;
	}
	
	/**
	 * Constructor using GeneHitGroup and tx name
	 * 
	 * @param geneHitGroup the gene group to which this tx group belongs
	 * @param name tx name
	 */
	private TxHitGroup(GeneHitGroup geneHitGroup, String name) {
		super();
		this.geneHitGroup = geneHitGroup;
		this.name = name;
		this.txHitMap = new HashMap<TxModel, TxHit>();
	}
	
	/**
	 * Adds a TxHit to the group.
	 * 
	 * This should only be called by TxHit.newInstance()
	 * 
	 * @throws IllegalArgumentException if the hit is for a different group
	 */
	void addTxHit(TxHit txHit) {
		if(txHit.getTxHitGroup() != this) {
			throw new IllegalArgumentException("Cannot add txHit " + txHit +
					" to this hit group, as it belongs to another.");
		}
		txHitMap.put(txHit.getTxModel(), txHit);
	}

	/**
	 * Gets the gene-level hit group to which this tx-level hit group belongs.
	 */
	public GeneHitGroup getGeneHitGroup() {
		return geneHitGroup;
	}

	/**
	 * Gets the tx name
	 */
	public String getName() {
		return name;
	}
	
	/**
	 * Gets an unmodifiable collection of all TxHits
	 */
	public Collection<TxHit> getAllTxHits() {
		return Collections.unmodifiableCollection(txHitMap.values());
	}

	/**
	 * Gets the gene name
	 */
	public String getGeneName() {
		return geneHitGroup.getName();
	}
	
	/**
	 * Gets an unmodifiable map from TxModel to TxHit
	 */
	public Map<TxModel, TxHit> getTxHitMap() {
		return Collections.unmodifiableMap(txHitMap);
	}

}
