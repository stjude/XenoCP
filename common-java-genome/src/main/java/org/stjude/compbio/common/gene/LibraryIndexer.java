package org.stjude.compbio.common.gene;

import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.index.RefIntervalIndex;

/**
 * Class that creates indexes from libraries.
 */
public class LibraryIndexer {
	public static enum Key {
		EXON, CDS, TX
	}
	
	/**
	 * Whether or not to create stranded indexes
	 */
	private boolean stranded;
	/**
	 * Key to use
	 */
	private Key key;
	
	/**
	 * Constructs an indexer with the indicated strandedness and key
	 * 
	 * @param stranded whether or not to produce stranded indexes
	 * @param key which interval should be used in the interval trees
	 */
	public LibraryIndexer(boolean stranded, Key key) {
		this.stranded = stranded;
		this.key = key;
	}

	/**
	 * Indexes a GeneHitLibrary at the TxHit level
	 * 
	 * @param library the library to index
	 * @return a RefIntervalIndex of TxHits
	 */
	public RefIntervalIndex<TxHit> newTxHitIndex(GeneHitLibrary library) {
		RefIntervalIndex<TxHit> index = new RefIntervalIndex<TxHit>(stranded);
		for(GeneHitGroup geneGrp: library.getAllGeneHitGroups()) {
			for(TxHitGroup txGrp: geneGrp.getAllTxHitGroups()) {
				for(TxHit txHit: txGrp.getAllTxHits()) {
					addToIndex(index, txHit.getTxModel(), txHit);
				}
			}
		}
		return index;
	}

	/**
	 * Indexes a GeneModelLibrary at the Spliceform level
	 * 
	 * @param library the library to index
	 * @return a RefIntervalIndex of Spliceforms
	 */
	public RefIntervalIndex<Spliceform> newSpliceformIndex(
			GeneModelLibrary library) {
		RefIntervalIndex<Spliceform> index =
				new RefIntervalIndex<Spliceform>(stranded);
		for(GeneModel geneModel: library.getGeneModels()) {
			for(Spliceform spliceform: geneModel.getSpliceforms()) {
				addToIndex(index, spliceform.getTxModel(), spliceform);
			}
		}
		return index;
	}

	/**
	 * Adds to the index based on the key
	 * 
	 * @param index index to add to
	 * @param txModel the tx model that is being added
	 * @param data the data object associated with it
	 */
	private <D> void addToIndex(
			RefIntervalIndex<D> index, TxModel txModel, D data) {
		switch(key) {
		case TX:
			index.add(txModel.getTxInterval(), data);
			break;
		case CDS:
			index.add(txModel.getCdsInterval(), data);
			break;
		case EXON:
			for(RefInterval exon: txModel.getExons()) {
				index.add(exon, data);
			}
		}
	}
}
