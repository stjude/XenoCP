package org.stjude.compbio.common.gene;

import org.stjude.compbio.common.interval.RefInterval;

/**
 * A named wrapper around a TxModel that belongs to a gene model
 */
public class Spliceform {
	/**
	 * The gene to which this tx belongs
	 */
	private GeneModel geneModel;
	/**
	 * The tx name
	 */
	private String name;
	/**
	 * The model with all of the positional information
	 */
	private TxModel txModel;
	
	/**
	 * Creates a new Spliceform instance, and adds the new instance to the
	 * GeneModel.
	 * 
	 * @param geneModel the gene model to which this spliceform belongs
	 * @param name spliceform/tx name
	 * @param txModel tx model
	 * @return the new instance
	 */
	static Spliceform newInstance(GeneModel geneModel, String name,
			TxModel txModel) {
		Spliceform spliceform = new Spliceform(geneModel, name, txModel);
		geneModel.addSpliceform(spliceform);
		return spliceform;
	}
	
	private Spliceform(GeneModel geneModel, String name, TxModel txModel) {
		super();
		this.geneModel = geneModel;
		this.name = name;
		this.txModel = txModel;
	}

	/**
	 * Gets the gene model to which this spliceform belongs
	 */
	public GeneModel getGeneModel() {
		return geneModel;
	}

	/**
	 * Gets the gene name
	 */
	public String getGeneName() {
		return geneModel.getName();
	}
	
	/**
	 * Gets the spliceform/tx name
	 */
	public String getName() {
		return name;
	}

	/**
	 * Gets the TxModel that gives the locations.
	 */
	public TxModel getTxModel() {
		return txModel;
	}
	
	/**
	 * Gets the model's tx interval
	 */
	public RefInterval getTxInterval() {
		return txModel.getTxInterval();
	}

	@Override
	public String toString() {
		return "Spliceform [name=" + name + "]";
	}
}
