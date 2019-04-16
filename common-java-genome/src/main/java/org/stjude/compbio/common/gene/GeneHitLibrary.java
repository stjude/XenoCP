package org.stjude.compbio.common.gene;

import java.util.Collection;
import java.util.List;
import java.util.Set;

public interface GeneHitLibrary {
	/**
	 * Gets a set of all gene names
	 */
	Set<String> getAllGeneNames();
	
	/**
	 * Gets all GeneHitGroups in the library
	 * 
	 * @return collection of GeneHitGroups
	 */
	Collection<GeneHitGroup> getAllGeneHitGroups();

	/**
	 * Gets the gene for a gene name
	 * 
	 * @return GeneHitGroup, or null if none found
	 */
	GeneHitGroup getGeneHitGroup(String name);

	/**
	 * Gets all TxHitGroups for a name.
	 * 
	 * Typically, there should only be one TxHitGroup for a name, but that is
	 * not enforced.
	 * 
	 * @return list of TxHitGroups, empty if there are none
	 */
	List<TxHitGroup> getTxHitGroups(String name);

	/**
	 * Gets the single TxHitGroup for a name.
	 * 
	 * Typically, there should only be one TxHitGroup for a name, but that is
	 * not enforced.  This method will throw an exception if there is more than
	 * one.
	 * 
	 * @return TxHitGroup, or null if there are none
	 * @throws IllegalStateException if there are multiple TxHitGroups with the
	 * 								 same name.
	 */
	TxHitGroup getTxHitGroup(String name) throws IllegalStateException;

}