package org.stjude.compbio.common.interval.index;

import org.stjude.compbio.common.interval.tree.IntervalTree;

/**
 * Indicates a tree should be rebuilt if the number of modifications is more
 * than both 10 and 1/10 of the tree size.
 */
public class TenTenTreeRebuildStrategy implements TreeRebuildStrategy {
	private IntervalTree<?, ?, ?> tree = null;
	private int mods = 0;
	
	public void setTree(IntervalTree<?, ?, ?> tree) {
		this.tree = tree;
	}
	
	public void rebuilt() {
		mods = 0;
	}

	public void added(int num) {
		mods += num;
	}

	public void removed(int num) {
		mods += num;
	}

	public boolean isTreeStale() {
		if(tree == null) {
			throw new IllegalStateException("No tree has been set");
		}
		return mods >= 10 && 10 * mods >= tree.size();
	}

}
