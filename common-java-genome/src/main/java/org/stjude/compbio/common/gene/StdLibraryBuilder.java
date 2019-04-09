package org.stjude.compbio.common.gene;

/**
 * Helper class for building StdGeneModelLibrary and StdGeneHitLibrary objects.
 * 
 * This class exists because the adding semantics for both are similar, even
 * though the underlying assumptions and data model differ in significant ways.
 * 
 * To use the class, specify the StdGeneModelLibrary and/or StdGeneHitLibrary
 * that you want to build up, and then call the add method for entries you want
 * to add.  If strict is on, then add calls may throw IllegalArgumentExceptions,
 * but if it's off, then the corresponding library will be disabled. 
 */
public class StdLibraryBuilder {
	private boolean strict;
	private StdGeneModelLibrary geneModelLibrary;
	private StdGeneHitLibrary geneHitLibrary;
	
	/**
	 * Constructs a library builder with the given strictness.
	 */
	public StdLibraryBuilder(boolean strict) {
		this.strict = strict;
	}

	/**
	 * Sets strictness
	 */
	public void setStrict(boolean strict) {
		this.strict = strict;
	}

	/**
	 * Starts a new StdGeneModelLibrary
	 */
	public void startGeneModelLibrary() {
		this.geneModelLibrary = new StdGeneModelLibrary();
	}

	/**
	 * Sets the StdGeneModelLibrary to built up
	 */
	public void setGeneModelLibrary(StdGeneModelLibrary stdGeneModelLibrary) {
		this.geneModelLibrary = stdGeneModelLibrary;
	}

	/**
	 * Gets the StdGeneModelLibrary (will be null if disabled or a problem was
	 * encountered in non-strict mode)
	 * 
	 * @return GeneModelLibarary instance, or null
	 */
	public StdGeneModelLibrary getGeneModelLibrary() {
		return geneModelLibrary;
	}
	
	/**
	 * Starts a new StdGeneHitLibrary
	 */
	public void startGeneHitLibrary() {
		this.geneHitLibrary = new StdGeneHitLibrary();
	}

	/**
	 * Sets the StdGeneHitLibrary to build up
	 */
	public void setGeneHitLibrary(StdGeneHitLibrary StdGeneHitLibrary) {
		this.geneHitLibrary = StdGeneHitLibrary;
	}

	/**
	 * Gets the StdGeneHitLibrary (will be null if disabled or a problem was
	 * encountered in non-strict mode)
	 * 
	 * @return StdGeneHitLibrary instance, or null
	 */
	public StdGeneHitLibrary getGeneHitLibrary() {
		return geneHitLibrary;
	}
	
	/**
	 * Adds to any and all underlying libraries
	 * 
	 * If there is a problem and strict is false, then the problematic library
	 * is disabled.
	 * 
	 * @param geneName the gene name
	 * @param txName the tx name
	 * @param txModel the tx model
	 * @return boolean if it was added to any underlying library
	 * @throws IllegalArgumentException if strict and it is thrown from
	 * 		underlying library
	 */
	public boolean addTxModel(String geneName, String txName, TxModel txModel) {
		boolean anyChange = false;
		if(geneModelLibrary != null) {
			try {
				anyChange = geneModelLibrary.addSpliceform(
						geneName, txName, txModel) || anyChange;
			} catch(IllegalArgumentException e) {
				if(strict) throw e;
				else geneModelLibrary = null;
			}
		}
		if(geneHitLibrary != null) {
			try {
				anyChange = geneHitLibrary.addTxHit(
						geneName, txName, txModel) || anyChange;
			} catch(IllegalArgumentException e) {
				if(strict) throw e;
				else geneHitLibrary = null;
			}
		}
		return anyChange;
	}
	
	/**
	 * Returns true if neither library is enabled
	 */
	public boolean isDegenerate() {
		return geneModelLibrary == null && geneHitLibrary == null;
	}
}
