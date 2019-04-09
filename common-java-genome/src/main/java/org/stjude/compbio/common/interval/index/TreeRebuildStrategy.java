package org.stjude.compbio.common.interval.index;

import org.stjude.compbio.common.interval.tree.IntervalTree;

/**
 * Strategy for determining when the tree inside an index should be rebuilt.
 * 
 * To use this, you must initially call setTree(), and call it again whenever
 * a new tree instance is created.  Call rebuilt() whenever the tree is
 * rebuilt.
 * 
 * Also, call added() on every addition and removed() on ever removal.  When
 * you call removed(), only specify the number of entries that were actually
 * removed.
 * 
 * After any operation that modifies the tree and any calls to added/removed,
 * and before performing any querying, call isTreeStale().  If it reports true,
 * then rebuild the tree.
 */
public interface TreeRebuildStrategy  {
	/**
	 * Provides the tree to the strategy, so that it may make direct inquiries
	 */
	void setTree(IntervalTree<?, ?, ?> tree);
	
	/**
	 * Tells the strategy the tree has been rebuilt.
	 */
	void rebuilt();
	
	/**
	 * Tells the strategy that one or more entries have been added
	 * 
	 * @param num number of entries added
	 */
	void added(int num);
	
	/**
	 * Tells the strategy that one or more entries have been removed
	 * 
	 * @param num number of entries removed
	 */
	void removed(int num);
	
	/**
	 * Reports whether or not the tree should be rebuilt.
	 */
	boolean isTreeStale();
}
