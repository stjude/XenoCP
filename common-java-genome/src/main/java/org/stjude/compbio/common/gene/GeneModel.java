package org.stjude.compbio.common.gene;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.stjude.compbio.common.interval.NamedRefInterval;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.StdNamedRefInterval;

/**
 * A gene model, which has a name and location and contains multiple spliceforms
 */
public class GeneModel {
	/**
	 * Gene name
	 */
	private String name;
	/**
	 * Spliceforms, stored as a map from their names to allow for lookups.
	 */
	private Map<String, Spliceform> spliceformMap;
	/**
	 * Gene interval, which is the union of all spliceform intervals
	 */
	private NamedRefInterval geneInterval;
	
	GeneModel(String name) {
		this.name = name;
		this.spliceformMap = new HashMap<String, Spliceform>();
		this.geneInterval = null;
	}

	/**
	 * Adds a spliceform to the model.
	 * 
	 * This should only be called by Spliceform.newInstance()
	 * 
	 * @param spliceform the new spliceform
	 * @throws IllegalArgumentException if the spliceform is for a different 
	 * 									gene model
	 */
	void addSpliceform(Spliceform spliceform) {
		if(spliceform.getGeneModel() != this) {
			throw new IllegalArgumentException("Cannot add spliceform " + 
					spliceform +
					" to this gene model, as it belongs to another.");
		}
		spliceformMap.put(spliceform.getName(), spliceform);
		
		// Update gene interval
		RefInterval txInterval = spliceform.getTxInterval();
		if(geneInterval == null) geneInterval = new StdNamedRefInterval(txInterval, name);
		else geneInterval = geneInterval.union(txInterval);
	}
	
	/**
	 * Gets the gene name
	 */
	public String getName() {
		return name;
	}
	
	/**
	 * Gets an interval covering all spliceforms
	 */
	public NamedRefInterval getInterval() {
		return geneInterval;
	}

	/**
	 * Gets an unmodifiable map mapping names to spliceforms
	 */
	public Map<String, Spliceform> getSpliceformMap() {
		return Collections.unmodifiableMap(spliceformMap);
	}
	
	/**
	 * Gets an unmodifiable collection of all spliceforms
	 */
	public Collection<Spliceform> getSpliceforms() {
		return Collections.unmodifiableCollection(spliceformMap.values());
	}
	
	/**
	 * Returns a named spliceform.
	 * 
	 * @param name spliceform name
	 * @return the Spliceform, or null if not found
	 */
	public Spliceform getSpliceform(String name) {
		return spliceformMap.get(name);
	}
	
	/**
	 * Returns the named TxModel 
	 * 
	 * @param name spliceform name
	 * @return the Spliceform's TxModel, or null if not found
	 */
	public TxModel getTxModel(String name) {
		Spliceform spliceform = spliceformMap.get(name);
		return (spliceform == null) ? null : spliceform.getTxModel();
	}
}
