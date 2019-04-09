package org.stjude.compbio.sam.block;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.stjude.compbio.sam.MdElt;
import org.stjude.compbio.sam.MdElts;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * List of blocks that represents an alignment.
 * 
 * A Blocklist may represent the entire alignment of a read or a pair of reads
 * or some other list of blocks.
 */
public class Blocklist implements Iterable<Block> {
	/**
	 * Leftmost block in the list (null if list is empty)
	 */
	private Block head;
	/**
	 * Rightmost block in the list (null if list is empty)
	 */
	private Block tail;
	/**
	 * Reference name
	 */
	private String refName;
	/**
	 * Left end of head block on read
	 */
	private long readLeft;
	/**
	 * Left end of head block on ref
	 */
	private long refLeft;
	/**
	 * All read bases
	 */
	private byte[] readBases;
	/**
	 * All read quals
	 */
	private byte[] readQuals;
	/**
	 * Object used to fetch reference bases
	 */
	private RefBasesFetcher refFetcher;
	
	/**
	 * A SAMRecord associated with the blocklist.
	 * 
	 * The record is used for tracking purposes, to put a read name in the
	 * toString() output.  That is the only use of the record.
	 */
	private SAMRecord record;

	/**
	 * Constructs a record-less, referenceless blocklist
	 * 
	 * Usually this is not the constructor you want to use.  Any operations
	 * involving sequence will throw exceptions when using a blocklist from
	 * this constructor!
	 * 
	 * @param refLeft 
	 * @param readLeft
	 */
	Blocklist(long refLeft, long readLeft) {
		this.record = null;
		this.refFetcher = null;
		
		this.head = null;
		this.tail = null;
		this.refName = null;
		this.refLeft = refLeft;
		this.readLeft = readLeft;
		this.readBases = new byte[0];
		this.readQuals = new byte[0];
	}

	/**
	 * Constructs an empty blocklist using a record and a RefBasesFetcher
	 * 
	 * Note that CIGAR and MD are not parsed, and no blocks are added.
	 * 
	 * @param record originating record
	 * @param refFetcher ref bases fetcher
	 */
	Blocklist(SAMRecord record, RefBasesFetcher refFetcher) {
		this.record = record;
		this.refFetcher = refFetcher;
		
		this.head = null;
		this.tail = null;
		this.refName = record.getReferenceName();
		this.refLeft = record.getAlignmentStart() - 1;
		this.readLeft = 0;
		this.readBases = record.getReadBases();
		this.readQuals = record.getBaseQualities();
	}

	/**
	 * Gets the SAMRecord that was previously set on this blocklist.
	 * 
	 * Remember that the SAMRecord is updatable and is not used by the Blocklist
	 * in any way.  This is purely informational and for convenience.
	 */
	public SAMRecord getRecord() {
		return record;
	}

	/**
	 * Sets the SAMRecord.
	 * 
	 * Remember that the SAMRecord is for reference only.
	 */
	public void setRecord(SAMRecord record) {
		this.record = record;
	}

	Block getHead() {
		return head;
	}

	Block getTail() {
		return tail;
	}

	/**
	 * Gets a block by index
	 * 
	 * @param index 0-based index in the list
	 */
	public Block get(int index) {
		int i = 0;
		for(Block block = head; block != null; block = block.rightBlock) {
			if(i == index) return block;
			++i;
		}
		throw new IndexOutOfBoundsException("Index " + index +
				" out of bounds");
	}
	
	public String getRefName() {
		return refName;
	}
	
	/**
	 * Sets reference name and left end (use this to "translate" from one
	 * coordinate system to another).
	 * 
	 * Use carefully
	 */
	public void translate(String refName, long refLeft) {
		this.refName = refName;
		this.refLeft = refLeft;
	}

	public long getReadLeft() {
		return readLeft;
	}

	public long getRefLeft() {
		return refLeft;
	} 

	public long getReadRight() {
		return tail.getReadRight();
	}

	public long getRefRight() {
		return tail.getRefRight();
	}

	byte[] getReadBases() {
		return readBases;
	}

	byte[] getReadQuals() {
		return readQuals;
	}

	RefBasesFetcher getRefFetcher() {
		return refFetcher;
	}

	/**
	 * Turns on sticky ref bases using the current non-sticky ref bases, or
	 * turns them off.
	 * 
	 * Be careful when using this.  Sticky bases are only appropriate when
	 * you are manipulating a blocklist in such a way that the resulting
	 * reference sequence is known in advance but is not available as a
	 * FASTA file.  This is a rare situation.
	 * 
	 * @param boolean whether to turn on or off
	 */
	public void setStickyRefBases(boolean stickyOn) {
		for(Block block: this) block.setStickyRefBases(stickyOn);
	}

	/**
	 * Returns the size of the underlying list
	 * 
	 * This is the same as calling getBlocks().size()
	 */
	public int size() {
		int size = 0;
		for(Block block = head; block != null; block = block.getRightBlock()) {
			++size;
		}
		return size;
	}
	
	/**
	 * Returns true if there are no blocks in the list
	 * 
	 * This is the same as calling getBlocks().isEmpty()
	 */
	public boolean isEmpty() {
		return head == null;
	}
	
	/**
	 * Adds a block to the far right side without any shifting.
	 * 
	 * @param length length of new block
	 * @param operator CigarOperator of new block
	 * @return newly added block
	 */
	public Block add(int length, CigarOperator operator) {
		Block block = new Block(this, tail, null, length, operator);
		if(tail != null) tail.rightBlock = block;
		if(head == null) head = block;
		tail = block;
		return block;
	}
	
	/**
	 * Adds a block before the indicated block, and shifts everything to the
	 * right.
	 * 
	 * @param block block to insert before
	 * @param length length of new block
	 * @param operator CigarOperator of new block
	 * @return the new block
	 */
	public Block insert(Block block, int length, CigarOperator operator) {
		Block newBlock = new Block(
				this, block.leftBlock, block, length, operator);
		if(block == head) {
			head = newBlock;
		}
		else {
			block.leftBlock.rightBlock = newBlock;
		}
		block.leftBlock = newBlock;
		return newBlock;
	}
	
	/**
	 * Removes a block, and shifts blocks on the right over to the left to fill
	 * 
	 * @param block block to remove
	 * @return block that was to the right of the removed block, or null
	 */
	public Block remove(Block block) {
		// Make sure it's in this list
		if(block.getBlocklist() != this) {
			throw new IllegalArgumentException(
					"Block " + block + " is not in this blocklist");
		}
		
		Block out = block.rightBlock;
		
		if(block == head && block == tail) {
			head = null;
			tail = null;
		}
		else if(block == head) {
			head = block.rightBlock;
			head.leftBlock = null;
		}
		else if(block == tail) {
			tail = block.leftBlock;
			tail.rightBlock = null;
		}
		else {
			block.leftBlock.rightBlock = block.rightBlock;
			block.rightBlock.leftBlock = block.leftBlock;
		}
		
		block.removed();
		return out;
	}

	/**
	 * Returns an iterator over the blocks, from left to right.
	 */
	public Iterator<Block> iterator() {
		return new Iterator<Block>() {
			private Block prev = null;
			private Block next = head;
			
			public boolean hasNext() {
				return next != null;
			}

			public Block next() {
				prev = next;
				next = next.getRightBlock();
				return prev;
			}

			public void remove() {
				Blocklist.this.remove(prev);
				prev = null;
			}
		};
	}

	/**
	 * Returns an iterator over the blocks, from right to left.
	 */
	public Iterator<Block> descendingIterator() {
		return new Iterator<Block>() {
			private Block prev = null;
			private Block next = tail;
			
			public boolean hasNext() {
				return next != null;
			}

			public Block next() {
				prev = next;
				next = next.getLeftBlock();
				return prev;
			}

			public void remove() {
				Blocklist.this.remove(prev);
				prev = null;
			}
		};
	}
	
	/**
	 * Converts the Blocklist to a Cigar
	 */
	public Cigar toCigar() {
		List<CigarElement> elts = new LinkedList<CigarElement>();
		for(Block block: this) elts.add(block.getCigarElement());
		return new Cigar(elts);
	}
	
	/**
	 * Converts the Blocklist to an MD element list
	 */
	public List<MdElt> toMdEltList() {
		MdElts.ListBuilder builder = MdElts.newListBuilder();
		for(Block block: this) {
			// Get the form
			BlockForm form = block.getForm();
			
			// Only build MD from ref-covering blocks
			if(!form.coversRef) continue;
			
			// Get refBases, as we'll need them regardless
			byte[] refBases = block.getRefBases();
			
			// Deletion is the easy case
			if(!form.coversRead) {
				builder.addDeleted(refBases);
				continue;
			}
			
			// Now, handle match/mismatch
			byte[] readBases = block.getReadBases();
			for(int i = 0; i < block.length(); i++) {
				if(readBases[i] == refBases[i]) {
					builder.addMatch();
				}
				else {
					builder.addMismatch(refBases[i]);
				}
			}
		}
		return builder.getList();
	}
	
	/**
	 * Gets the in-base reference left end of the first alignment block
	 * @return in-base left end, or 0 if there are no aligned blocks
	 */
	public int getAlignmentStart() {
		for(Block block: this) {
			if(block.isAligned()) {
				return block.getRefInterval().getLeftInBaseInt();
			}
		}
		return 0;
	}

	/**
	 * Shifts the entire blocklist by the given distance
	 *  
	 * @param readDistance amount to shift on read (negative is left)
	 * @param refDistance amount to shift on ref (negative is left)
	 */
	public void shift(long readDistance, long refDistance) {
		readLeft += readDistance;
		refLeft += refDistance;
	}

}
