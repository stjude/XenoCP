package org.stjude.compbio.common.gene;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.stjude.compbio.util.HashMultimap;
import org.stjude.compbio.util.Multimap;

/**
 * Library of gene/tx hit information
 */
public class StdGeneHitLibrary implements GeneHitLibrary {
	/**
	 * GeneHitGroups, stored as a map from their names to the instances, to
	 * allow for lookups by name.
	 */
	private Map<String, GeneHitGroup> geneHitGroupMap;
	/**
	 * Multimap from tx name to TxHitGroups (in case there is a tx name shared
	 * across multiple genes--which should not happen, but can)
	 */
	private Multimap<String, TxHitGroup, List<TxHitGroup>> txHitGroupMap;
	
	/**
	 * Constructs an empty GeneHitLibrary
	 */
	public StdGeneHitLibrary() {
		this.geneHitGroupMap = new HashMap<String, GeneHitGroup>();
		this.txHitGroupMap = HashMultimap.newListInstance();
	}
	
	public Set<String> getAllGeneNames() {
		return geneHitGroupMap.keySet();
	}

	public Collection<GeneHitGroup> getAllGeneHitGroups() {
		return geneHitGroupMap.values();
	}

	public GeneHitGroup getGeneHitGroup(String name) {
		return geneHitGroupMap.get(name);
	}

	public List<TxHitGroup> getTxHitGroups(String name) {
		return txHitGroupMap.getLazyInit(name);
	}

	public TxHitGroup getTxHitGroup(String name) throws IllegalStateException {
		List<TxHitGroup> txHitGroups = txHitGroupMap.get(name);
		if(txHitGroups == null || txHitGroups.isEmpty()) {
			return null;
		}
		else if(txHitGroups.size() > 1) {
			throw new IllegalStateException("Found " + txHitGroups.size() +
					" tx hit groups with name " + name);
		}
		else {
			return txHitGroups.get(0);
		}
	}

	/**
	 * Adds a new, empty GeneHitGroup with the given name to the library
	 * 
	 * @param geneName the gene name
	 * @return the new group
	 */
	private GeneHitGroup addNewGeneHitGroup(String geneName) {
		GeneHitGroup newGrp = new GeneHitGroup(geneName);
		geneHitGroupMap.put(geneName, newGrp);
		return newGrp;
	}
	
	/**
	 * Gets an existing gene hit group or creates a new empty one
	 * 
	 * @param geneName gene name
	 * @return new or existing group
	 */
	private GeneHitGroup getOrAddGeneHitGroup(String geneName) {
		if(geneHitGroupMap.containsKey(geneName)) {
			return geneHitGroupMap.get(geneName);
		}
		else {
			return addNewGeneHitGroup(geneName);
		}
	}

	/**
	 * Adds a new, empty TxHitGroup with the given name to the library
	 * 
	 * @param geneGrp the gene hit group
	 * @param txName the tx name
	 * @return the new group
	 */
	private TxHitGroup addNewTxHitGroup(GeneHitGroup geneGrp, String txName) {
		TxHitGroup newGrp = TxHitGroup.newInstance(geneGrp, txName);
		txHitGroupMap.add(txName, newGrp);
		return newGrp;
	}
	
	/**
	 * Gets an existing tx hit group or creates a new empty one.
	 * 
	 * @param geneName the gene name
	 * @param txName the tx name
	 * @return new or existing group
	 */
	private TxHitGroup getOrAddTxHitGroup(String geneName, String txName) {
		GeneHitGroup geneGrp = getOrAddGeneHitGroup(geneName);
		if(geneGrp.getTxHitGroupMap().containsKey(txName)) {
			return geneGrp.getTxHitGroup(txName);
		}
		else {
			return addNewTxHitGroup(geneGrp, txName);
		}
	}
	
	/**
	 * Adds a new TxHit to the library
	 * 
	 * @param txGrp the tx hit group
	 * @param txModel the defining model
	 * @return the new TxHit
	 */
	private TxHit addNewTxHit(TxHitGroup txGrp, TxModel txModel) {
		return TxHit.newInstance(txGrp, txModel);
	}

	/**
	 * Adds a new TxHit.
	 * 
	 * @param geneName the gene name
	 * @param txName the tx name
	 * @param txModel the tx model
	 * @return true if one was added, false if it already existed
	 */
	public boolean addTxHit(String geneName, String txName, TxModel txModel) {
		TxHitGroup txGrp = getOrAddTxHitGroup(geneName, txName);
		if(txGrp.getTxHitMap().containsKey(txModel)) {
			return false;
		}
		else {
			addNewTxHit(txGrp, txModel);
			return true;
		}
	}
}
