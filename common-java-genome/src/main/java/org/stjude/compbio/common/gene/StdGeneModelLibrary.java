package org.stjude.compbio.common.gene;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Standard library of gene models, which allows modification.
 */
public class StdGeneModelLibrary implements GeneModelLibrary {
	/**
	 * GeneModels, stored as a map from their names to the instances, to allow
	 * for lookups by name.
	 */
	private Map<String, GeneModel> geneModelMap;
	/**
	 * Map from spliceform name to Spliceforms
	 */
	private Map<String, Spliceform> spliceformMap;

	/**
	 * Constructs an empty StdGeneModelLibrary
	 */
	public StdGeneModelLibrary() {
		this.geneModelMap = new HashMap<String, GeneModel>();
		this.spliceformMap = new HashMap<String, Spliceform>();
	}

	public Collection<GeneModel> getGeneModels() {
		return Collections.unmodifiableCollection(geneModelMap.values());
	}

	public GeneModel getGeneModel(String name) {
		return geneModelMap.get(name);
	}

	public Spliceform getSpliceform(String name) {
		return spliceformMap.get(name);
	}

	public TxModel getTxModel(String spliceformName) {
		Spliceform spliceform = spliceformMap.get(spliceformName);
		return (spliceform == null) ? null : spliceform.getTxModel();
	}

	/**
	 * Adds a new, empty GeneModel with the given name to the library
	 * 
	 * @param geneName the gene name
	 * @return the new model
	 */
	private GeneModel addNewGeneModel(String geneName) {
		GeneModel geneModel = new GeneModel(geneName);
		geneModelMap.put(geneName, geneModel);
		return geneModel;
	}

	/**
	 * Gets an existing gene model or creates a new empty one
	 * 
	 * @param geneName gene name
	 * @return new or existing GeneModel
	 */
	private GeneModel getOrAddGeneModel(String geneName) {
		if(geneModelMap.containsKey(geneName)) {
			return geneModelMap.get(geneName);
		} else {
			return addNewGeneModel(geneName);
		}
	}

	/**
	 * Adds a new Spliceform.
	 * 
	 * @param geneName the gene name
	 * @param txName the tx name
	 * @param txModel the tx model
	 * @return true if one was added, false if it already existed
	 * @throws IllegalArgumentException if there is a conflict
	 */
	public boolean addSpliceform(
			String geneName, String txName, TxModel txModel) {
		Spliceform spliceform;
		
		// See if the tx name is known, and if so, just check for conflicts
		if(spliceformMap.containsKey(txName)) {
			spliceform = spliceformMap.get(txName);
			if(!spliceform.getGeneName().equals(geneName)) {
				throw new IllegalArgumentException("Transcript " + txName + 
						" already exists in gene " + spliceform.getGeneName() +
						".  Cannot add in gene " + geneName);
			}
			else if(!spliceform.getTxModel().equals(txModel)) {
				throw new IllegalArgumentException("Transcript " + txName + 
						" already exists in gene " + geneName +
						" with a different tx model.");
			}
			return false;
		}
		
		// Otherwise, add it
		GeneModel geneModel = getOrAddGeneModel(geneName);
		spliceform = Spliceform.newInstance(geneModel, txName, txModel);
		spliceformMap.put(txName, spliceform);
		return true;
	}
}
