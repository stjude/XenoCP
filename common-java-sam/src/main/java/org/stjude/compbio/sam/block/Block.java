package org.stjude.compbio.sam.block;

import java.util.Arrays;

import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.interval.LongInterval;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.StdLongInterval;
import org.stjude.compbio.common.interval.StdRefInterval;
import org.stjude.compbio.util.Count;

import picard.PicardException;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

/**
 * Representation of a single "block" of an alignment between a SAMRecord and a
 * reference.
 * 
 * The block is tied to a region of the read and a region of the reference.
 * These regions will either be the same length (alignment) or one will be 0-
 * length (gap).
 * 
 * This can be built up from a record's CIGAR and either an MD tag or reference
 * sequence.  The building logic is not included in here.
 * 
 * Blocks may exist at various levels of granularity, but they must be granular
 * enough that a single CigarElement may be built from a single Block.  This
 * does allow for cases where there are multiple CigarElements in the source
 * that become a single block.  For example, if the source uses = and X in
 * CIGAR, multiple = or X may be joined to form one block for the whole aligned
 * region.
 * 
 * Blocks only exist and make sense in the context of a Blocklist, although one
 * could have a single-block list.  The Block class itself contains pointers to
 * the next and previous blocks, making it essentially an element in a doubly-
 * linked list.
 * 
 * @see Blocklist
 */
public class Block {
	/**
	 * The containing Blocklist
	 */
	Blocklist blocklist;
	/**
	 * The Block to the left
	 */
	Block leftBlock;
	/**
	 * The Block to the right
	 */
	Block rightBlock;
	
	/**
	 * The block length, which is the length of the block on the read or ref,
	 * whichever is non-zero.
	 */
	private int length;
	/**
	 * CigarElement length, which will differ from length iff form is NEITHER
	 */
	private int cigarLength;
	/**
	 * The CigarOperator
	 */
	private CigarOperator operator;
	/**
	 * The form (whether read and/or ref are covered), implied by the operator
	 */
	private BlockForm form;
	
	/**
	 * Ref bases that have been "stuck" to this block, null if ref bases are not
	 * sticky.
	 */
	private byte[] stickyRefBases;
	
	/**
	 * Constructs a new template block, which is only useful for adding or
	 * setting when using DirectionalBlockIterator.
	 * 
	 * @param length cigar length
	 * @param operator cigar operator
	 * @see DirectionalBlockIterator
	 */
	public static Block newTemplate(int length, CigarOperator operator) {
		return new Block(null, null, null, length, operator);
	}
	
	/**
	 * Constructs a new template block, which is only useful for adding or
	 * setting when using DirectionalBlockIterator.
	 * 
	 * @param length cigar length
	 * @param operatorName cigar operator enum name (MNDSHIPX and "EQ"), or "="
	 * @see DirectionalBlockIterator
	 */
	public static Block newTemplate(int length, String operatorName) {
		if(operatorName.equals("=")) operatorName = "EQ";
		return newTemplate(length, CigarOperator.valueOf(operatorName));
	}
	
	/**
	 * Constructs a block using blocklist, neighbors, length and op
	 * 
	 * @param blocklist
	 * @param leftBlock
	 * @param rightBlock
	 * @param length
	 * @param operator
	 */
	Block(Blocklist blocklist, Block leftBlock, Block rightBlock,
			int cigarLength, CigarOperator operator) {
		this.blocklist = blocklist;
		this.leftBlock = leftBlock;
		this.rightBlock = rightBlock;
		this.cigarLength = cigarLength;
		this.operator = operator;
		this.form = BlockForm.valueOf(operator);
		this.length = (form == BlockForm.NEITHER) ? 0 : cigarLength;
		this.stickyRefBases = null;
	}

	/**
	 * Constructs a block using blocklist, neighbors, and cigar elt
	 * 
	 * @param blocklist
	 * @param leftBlock
	 * @param rightBlock
	 * @param cigarElt
	 */
	Block(Blocklist blocklist, Block leftBlock, Block rightBlock,
			CigarElement cigarElt) {
		this(blocklist, leftBlock, rightBlock, cigarElt.getLength(),
				cigarElt.getOperator());
	}
	
	/**
	 * Indicates this block has been removed from its list
	 */
	void removed() {
		blocklist = null;
		leftBlock = null;
		rightBlock = null;
	}

	public Blocklist getBlocklist() {
		return blocklist;
	}

	public Block getLeftBlock() {
		return leftBlock;
	}

	public Block getRightBlock() {
		return rightBlock;
	}

	public CigarOperator getOperator() {
		return operator;
	}

	public int length() {
		return length;
	}

	public CigarElement getCigarElement() {
		return new CigarElement(cigarLength, operator);
	}

	public BlockForm getForm() {
		return form;
	}

	public boolean isAligned() {
		return form == BlockForm.BOTH;
	}
	
	/**
	 * Returns true for soft-clipping and hard trimming
	 */
	public boolean isClippingOrTrimming() {
		return operator == CigarOperator.S || operator == CigarOperator.H;
	}

	public boolean isHardTrim() {
		return operator == CigarOperator.H;
	}

	long getReadLeft() {
		return (leftBlock == null) ?
				blocklist.getReadLeft() : leftBlock.getReadRight();
	}
	
	long getReadLen() {
		return (form.advancesRead ? length : 0);
	}
	
	long getReadRight() {
		return getReadLeft() + getReadLen();
	}
	
	public LongInterval getReadInterval() {
		long left = getReadLeft();
		return new StdLongInterval(left, left + getReadLen());
	}

	long getRefLeft() {
		return (leftBlock == null) ?
				blocklist.getRefLeft() : leftBlock.getRefRight();
	}
	
	long getRefLen() {
		return (form.advancesRef ? length : 0);
	}
	
	long getRefRight() {
		return getRefLeft() + getRefLen();
	}
	
	public RefInterval getRefInterval() {
		long left = getRefLeft();
		return new StdRefInterval(blocklist.getRefName(), Strand.NA,
				left, left + getRefLen());
	}

	/**
	 * Returns true if there is reference information available.
	 * 
	 * If this returns false, then getRefBases() will return null, and
	 * getMatchCount() will not work.
	 */
	public boolean hasRefBases() {
		return stickyRefBases != null || blocklist.getRefFetcher() != null;
	}

	/**
	 * Returns the reference bases if available.
	 * 
	 * If there is any reference sequence information available, then this will
	 * be a non-null array, though it will be empty if the block does not 
	 * advance the reference.  If there is no reference sequence info (no seq,
	 * no MD), then null will be returned.  If there is reference sequence info,
	 * but the block advances but does not cover the reference, then null is
	 * returned.
	 * 
	 * @return reference base array, or null if no ref seq info is available
	 */
	public byte[] getRefBases() {
		// If there are sticky bases, then return them
		if(stickyRefBases != null) return stickyRefBases.clone();
		
		// Get the fetcher, which tells us if there might be sequence available
		RefBasesFetcher refFetcher = blocklist.getRefFetcher();
		if(refFetcher == null) return null;
		
		// If this does not cover any ref bases, then return empty array or
		// null, depending on whether it advances ref (if it advances ref, but
		// does not cover it, then there is some sequence there but we don't
		// have it, so that's why we return null).
		if(!form.coversRef) return form.advancesRef ? null : new byte[0];

		// Otherwise, delegate to the fetcher
		try {
			return refFetcher.getRefBases(getRefInterval());
		} catch(PicardException e) {
			throw new PicardException("Could not fetch ref bases for block " +
					this, e);
		}
	}

	/**
	 * Returns true if the ref bases are "sticky", i.e. they move with the
	 * block and do not change.
	 */
	public boolean hasStickyRefBases() {
		return stickyRefBases != null;
	}

	/**
	 * Sets the sticky ref bases to use, or null to turn off sticky bases
	 * 
	 * Be careful when using this.  Sticky bases are only appropriate when
	 * you are manipulating a blocklist in such a way that the resulting
	 * reference sequence is known in advance but is not available as a
	 * FASTA file.  This is a rare situation.
	 * 
	 * @param stickyRefBases new sticky ref bases
	 */
	public void setStickyRefBases(byte[] stickyRefBases) {
		this.stickyRefBases = stickyRefBases;
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
		this.stickyRefBases = null;
		if(stickyOn) this.stickyRefBases = getRefBases();
	}

	/**
	 * Returns the read bases.
	 * 
	 * @return read bases as byte array, which may be empty but never null
	 */
	public byte[] getReadBases() {
		return copyOfReadArray(blocklist.getReadBases());
	}

	/**
	 * Returns the read quals.
	 * 
	 * @return read quals as byte array, which may be empty but never null
	 */
	public byte[] getReadQuals() {
		return copyOfReadArray(blocklist.getReadQuals());
	}

	/**
	 * Pulls out the portion of a read array (bases or quals) that this block
	 * covers.
	 * 
	 * @param allReadBytes all of the read bases or quals
	 * @return the portion covered by this block
	 */
	private byte[] copyOfReadArray(byte[] allReadBytes) {
		byte[] readBytes;
		if(allReadBytes == null) {
			readBytes = new byte[0];
		}
		else {
			LongInterval readInterval = getReadInterval();
			readBytes = Arrays.copyOfRange(allReadBytes,
					readInterval.getLeftInt(), readInterval.getRightInt());
		}
		return readBytes;
	}

	/**
	 * Gets the count of matching bases out of the to total number.
	 * 
	 * @return number of matching bases as a Count object
	 * @throws IllegalStateException if this is not an aligned block, or if ref
	 * 				sequence is unavailable
	 */
	public Count getMatchCount() {
		// Check for alignment
		if(!isAligned()) {
			throw new IllegalStateException(
					"Cannot count matches in unaligned block.");
		}
		// If 0 length, then it's easy
		if(length == 0) return new Count(0, 0);
		// Get ref bases, and make sure it worked
		byte[] refBases = getRefBases();
		if(refBases == null) {
			throw new IllegalStateException(
					"Could not get reference bases");
		}
		
		// Now, do the counting
		byte[] readBases = getReadBases();
		int matches = 0;
		for(int i = 0; i < refBases.length; ++i) {
			if(refBases[i] == readBases[i]) ++matches;
		}
		return new Count(matches, length);
	}
	
	/**
	 * Splits this block into two.
	 * 
	 * @param leftLength length of the resulting left block
	 * @return the two new blocks
	 */
	public Block[] split(int leftLength) {
		// Replace this block with the two split blocks
		Block[] out = new Block[2];
		out[0] = blocklist.insert(this, leftLength, operator);
		out[1] = blocklist.insert(this, length - leftLength, operator);
		blocklist.remove(this);		
		
		// Handle sticky bases
		if(stickyRefBases != null && stickyRefBases.length == length) {
			out[0].setStickyRefBases(Arrays.copyOfRange(
					stickyRefBases, 0, leftLength));
			out[1].setStickyRefBases(Arrays.copyOfRange(
					stickyRefBases, leftLength, length));
		}
		
		return out;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if(blocklist == null) {
			sb.append("template:");
		}
		else {
			SAMRecord record = blocklist.getRecord();
			if(record != null) {
				sb.append(record.getReadName());
				if(record.getReadPairedFlag()) {
					sb.append("/");
					if(record.getFirstOfPairFlag()) sb.append("1");
					else if(record.getSecondOfPairFlag()) sb.append("2");
				}
				sb.append(":");
			}
		}
		sb.append(length);
		sb.append(operator);
		if(blocklist != null) {
			sb.append("@");
			sb.append(getRefInterval());
		}
		return sb.toString();
	}
}
