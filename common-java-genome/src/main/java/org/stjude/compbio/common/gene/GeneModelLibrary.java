package org.stjude.compbio.common.gene;

import java.util.Collection;

/**
 * A library of gene models
 */
public interface GeneModelLibrary {
	/**
	 * Gets all gene models
	 */
	Collection<GeneModel> getGeneModels();
	
	/**
	 * Gets the named GeneModel
	 */
	GeneModel getGeneModel(String name);
	
	/**
	 * Gets the named Spliceform
	 */
	Spliceform getSpliceform(String name);
	
	/**
	 * Gets the TxModel for the named Spliceform
	 */
	TxModel getTxModel(String spliceformName);
}
