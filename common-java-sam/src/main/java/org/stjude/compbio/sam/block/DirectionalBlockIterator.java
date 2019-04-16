package org.stjude.compbio.sam.block;

import java.util.ListIterator;
import java.util.NoSuchElementException;

/**
 * Directional ListIterator implementation over Blocklists.
 * 
 * There are 4 "directions" of iteration, based on two anchor sides, and two
 * directions relative to the anchor.
 * 
 * The anchor can be left or right.  This affects how list modifications
 * operate.  If the anchor is left then the left end of the blocklist never
 * changes.  If the anchor is right, then the right end of the blocklist never
 * changes.
 * 
 * If you are not making modifications, then you may use the aliases LEFT and
 * RIGHT to indicate the two left-anchored directions.
 * 
 * This class cannot necessarily handle changes made to the underlying
 * Blocklist, other than through the class itself. 
 */
public class DirectionalBlockIterator implements ListIterator<Block> {
	/**
	 * Directions
	 * 
	 * @see DirectionalBlockIterator
	 */
	public static enum Direction {
		FROM_LEFT_ANCHOR(true, true),
		TO_LEFT_ANCHOR(true, false),
		FROM_RIGHT_ANCHOR(false, false),
		TO_RIGHT_ANCHOR(false, true);
		
		public static final Direction LEFT = TO_LEFT_ANCHOR;
		public static final Direction RIGHT = FROM_LEFT_ANCHOR;
		
		public static final Direction[] FROM_OUTSIDE = {
			TO_RIGHT_ANCHOR, TO_LEFT_ANCHOR
		};
		
		private final boolean anchoredLeft;
		private final boolean fromLeft;
		private Direction(boolean anchoredLeft, boolean fromLeft) {
			this.anchoredLeft = anchoredLeft;
			this.fromLeft = fromLeft;
		}
		
		private Block first(Blocklist blocklist) {
			return fromLeft ? blocklist.getHead() : blocklist.getTail();
		}
		private Block prev(Block block) {
			return fromLeft ? block.leftBlock : block.rightBlock;
		}
		private Block next(Block block) {
			return fromLeft ? block.rightBlock : block.leftBlock;
		}
	}

	/**
	 * The blocklist
	 */
	private final Blocklist blocklist;
	/**
	 * Iteration direction
	 */
	private final Direction direction;
	/**
	 * Next block
	 */
	private Block next;
	/**
	 * Prev block
	 */
	private Block prev;
	/**
	 * Removable block
	 */
	private Block removable;
	
	/**
	 * Constructs a directional block iterator going in the given direction.
	 */
	public DirectionalBlockIterator(Blocklist blocklist, Direction direction) {
		this.blocklist = blocklist;
		this.direction = direction;
		
		this.next = direction.first(blocklist);
		this.prev = null;
		this.removable = null;
	}

	public boolean hasNext() {
		return next != null;
	}

	public Block next() {
		if(!hasNext()) throw new NoSuchElementException();
		removable = next;
		prev = removable;
		next = direction.next(removable);
		return removable;
	}

	public boolean hasPrevious() {
		return prev != null;
	}

	public Block previous() {
		if(!hasPrevious()) throw new NoSuchElementException();
		removable = prev;
		next = removable;
		prev = direction.prev(removable);
		return removable;
	}

	public int nextIndex() {
		if(next == null) {
			return direction.fromLeft ? blocklist.size() : -1;
		}
		else {
			return indexOf(next); 
		}
	}

	public int previousIndex() {
		if(prev == null) {
			return direction.fromLeft ? -1 : blocklist.size();
		}
		else {
			return indexOf(prev); 
		}
	}

	private int indexOf(Block block) {
		int index = 0;
		while(block.leftBlock != null) {
			block = block.leftBlock;
			++index;
		}
		return index;
	}

	public void remove() {
		// Check
		if(removable == null) throw new IllegalStateException();
		
		// Update next/prev before removable is destroyed
		next = direction.next(removable);
		prev = direction.prev(removable);
		
		// Do the removal
		if(!direction.anchoredLeft) {
			blocklist.shift(removable.getReadLen(), removable.getRefLen());
		}
		blocklist.remove(removable);
		
		// Clear removable
		removable = null;
	}

	public void set(Block blockTemplate) {
		// Check
		if(removable == null) throw new IllegalStateException();
		
		// Determine which cursor object needs to be updated at the end
		boolean replacingNext = (next == removable);
		
		// Do the remove first
		remove();
		
		// Add in the list, and reset removable to the added block
		removable = addInBlocklist(blockTemplate);

		// Reset the appropriate cursor object
		if(replacingNext) {
			next = removable;
		}
		else {
			prev = removable;
		}
	}

	public void add(Block blockTemplate) {
		// Add in the list, and set prev equal to the new block
		prev = addInBlocklist(blockTemplate);
		
		// Clear removable
		removable = null;
	}

	/**
	 * Performs an add at the cursor location in the blocklist, but does not
	 * update the iterator.
	 * 
	 * @param blockTemplate template of new block to add
	 * @return actual new block added
	 */
	private Block addInBlocklist(Block blockTemplate) {
		// If we are anchored right, then shift the list to the left
		if(!direction.anchoredLeft) {
			blocklist.shift(
					-1 * blockTemplate.getReadLen(),
					-1 * blockTemplate.getRefLen());
		}
		
		// Determine the right cursor block, for insertion
		Block right = direction.fromLeft ? next : prev;
		
		// Add/insert and reset previous
		if(right == null) {
			return blocklist.add(
					blockTemplate.length(), blockTemplate.getOperator());
		}
		else {
			return blocklist.insert(
					right, blockTemplate.length(), blockTemplate.getOperator());
		}
	}

}
