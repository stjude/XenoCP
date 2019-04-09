package org.stjude.compbio.common.gene;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

public class GeneHitGroup {
	/**
	 * Name, which should be unique among group instances
	 */
	private String name;
	/**
	 * TxHitGroups, stored as a map from their names to the instances, to
	 * allow for lookups by name.
	 */
	private Map<String, TxHitGroup> txHitGroupMap;
	
	/**
	 * Constructor using gene name
	 * 
	 * @param name gene name
	 */
	GeneHitGroup(String name) {
		this.name = name;
		this.txHitGroupMap = new HashMap<String, TxHitGroup>();
	}

	/**
	 * Adds a TxHitGroup to the gene group
	 * 
	 * This should only be called by TxHitGroup.newInstance()
	 * 
	 * @throws IllegalArgumentException if the tx group is for a different gene 
	 * 									group
	 */
	void addTxHitGroup(TxHitGroup txHitGroup) {
		if(txHitGroup.getGeneHitGroup() != this) {
			throw new IllegalArgumentException("Cannot add tx hit group " + 
					txHitGroup +
					" to this gene hit group, as it belongs to another.");
		}
		txHitGroupMap.put(txHitGroup.getName(), txHitGroup);
	}

	/**
	 * Gets the gene name
	 */
	public String getName() {
		return name;
	}
	
	/**
	 * Gets an unmodifiable collection of all TxHitGroups
	 */
	public Collection<TxHitGroup> getAllTxHitGroups() {
		return Collections.unmodifiableCollection(txHitGroupMap.values());
	}

	/**
	 * Gets an unmodifiable map from TxHitGroup name to TxHitGroup
	 */
	public Map<String, TxHitGroup> getTxHitGroupMap() {
		return Collections.unmodifiableMap(txHitGroupMap);
	}

	/**
	 * Gets the TxHitGroup with the given name
	 * 
	 * @param txName tx name
	 * @return the TxHitGroup, or null if there is none
	 */
	public TxHitGroup getTxHitGroup(String txName) {
		return txHitGroupMap.get(txName);
	}
	
}
