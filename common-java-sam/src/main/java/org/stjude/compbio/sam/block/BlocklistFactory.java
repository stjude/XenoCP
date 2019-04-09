package org.stjude.compbio.sam.block;

import java.util.Arrays;
import java.util.ListIterator;

import org.stjude.compbio.common.formats.fasta.IndexedFastaFetcher;
import org.stjude.compbio.sam.MdElt;
import org.stjude.compbio.sam.MdElts;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

/**
 * Factory class for creating Blocklists
 */
public class BlocklistFactory {
	/**
	 * Whether to use MD information if available
	 */
	private boolean mdParsingEnabled;
	/**
	 * IndexedFastaFetcher to use to get ref sequence, or null to not use
	 */
	private IndexedFastaFetcher refFetcher;
	
	/**
	 * Current md-based ref bases fetcher.
	 * 
	 * This is stored here because after construction it needs to be filled.
	 */
	private SparseRefBasesFetcher curMdFetcher;
	
	/**
	 * Constructs a BlocklistFactory with the given settings
	 * 
	 * @param mdParsingEnabled whether or not to use MD from the record
	 * @param refFetcher fetcher to use for fetching ref seq; null to disable
	 */
	public BlocklistFactory(boolean mdParsingEnabled,
			IndexedFastaFetcher refFetcher) {
		this.mdParsingEnabled = mdParsingEnabled;
		this.refFetcher = refFetcher;
	}

	public boolean isMdParsingEnabled() {
		return mdParsingEnabled;
	}

	public void setMdParsingEnabled(boolean mdParsingEnabled) {
		this.mdParsingEnabled = mdParsingEnabled;
	}

	public IndexedFastaFetcher getRefFetcher() {
		return refFetcher;
	}

	public void setRefFetcher(IndexedFastaFetcher refFetcher) {
		this.refFetcher = refFetcher;
	}

	/**
	 * Builds a RefBasesFetcher based on settings and the SAMRecord
	 * 
	 * @param record SAMRecord, may be null
	 * @return fetcher, may be null if no fetching is enabled
	 */
	private RefBasesFetcher getRefBasesFetcher(SAMRecord record) {
		// Make FASTA-based fetcher if appropriate
		IndexedFastaRefBasesFetcher fastaFetcher = null;
		if(refFetcher != null) {
			fastaFetcher = new IndexedFastaRefBasesFetcher(refFetcher);
		}
		
		// Make MD-based fetcher if appropriate
		curMdFetcher = null;
		if(mdParsingEnabled && record != null &&
				record.getAttribute(SAMTag.MD.name()) != null) {
			curMdFetcher = new SparseRefBasesFetcher(
					record.getReferenceName(), fastaFetcher);
		}
		
		// Return the proper fetcher
		return curMdFetcher == null ? fastaFetcher : curMdFetcher;
	}
	
	/**
	 * Constructs an empty Blocklist with a read/ref position.
	 * 
	 * @param refLeft the interbase ref left position of the first added block
	 * @param readLeft the interbase read left position of the first added block
	 */
	public Blocklist newEmptyBlocklist(int refLeft, int readLeft) {
		return new Blocklist(refLeft, readLeft);
	}
	
	/**
	 * Constructs a reference-less Blocklist by parsing Cigar
	 * 
	 * You should generally parse a record instead.
	 * 
	 * @param cigar the cigar to parse
	 * @return the Blocklist (ref start == 0)
	 */
	public Blocklist newReferencelessBlocklist(Cigar cigar) {
		Blocklist blocklist = new Blocklist(0, 0);
		for(CigarElement ce: cigar.getCigarElements()) {
			blocklist.add(ce.getLength(), ce.getOperator());
		}
		return blocklist;
	}

	/**
	 * Constructs a new Blocklist by parsing a record.
	 * 
	 * Most of the work here is in handling the MD tag if that's enabled.
	 * 
	 * @param record the record to parse
	 * @return the Blocklist (empty if record has no CIGAR)
	 */
	public Blocklist newBlocklist(SAMRecord record) {
		// Construct fetcher
		RefBasesFetcher fetcher = getRefBasesFetcher(record);
		Blocklist blocklist = new Blocklist(record, fetcher);
		
		// getRefBasesFetcher() will have set curMdFetcher, if non-null, then
		// we are using MD
		boolean useMd = curMdFetcher != null;
		
		// Check for null-CIGAR case
		Cigar cigar = record.getCigar();
		if(cigar == null || cigar.getCigarElements().isEmpty()) {
			return blocklist;
		}
		
		// Get list iterator over md elements, as well as string builder to hold
		// all ref bases not yet added to sparse fetcher
		ListIterator<MdElt> mdListIter = null;
		StringBuilder refBases = null;
		long refBasesLeft = blocklist.getRefLeft();
		if(useMd) {
			mdListIter = MdElts.parse(record).listIterator();
			refBases = new StringBuilder();
		} 
		
		// Iterate over the CigarElements
		for(CigarElement cigarElt: cigar.getCigarElements()) {
			// Add the block to the list
			Block block = blocklist.add(
					cigarElt.getLength(), cigarElt.getOperator());
			
			// If MD is enabled, then add to MD bases if appropriate
			if(useMd) {
				if(block.getForm().coversRef) {
					addToKnownBasesBuilder(refBases, mdListIter, block);
				}
				else if(cigarElt.getOperator() == CigarOperator.N) { 
					// Add known bases
					addKnownBasesToMdFetcher(refBases, refBasesLeft);
					// Reset
					refBases.delete(0, refBases.length());
					refBasesLeft = block.getRefRight();
				}
			}
		}
		
		// Add final set of known bases
		if(useMd) addKnownBasesToMdFetcher(refBases, refBasesLeft);
		
		// Return the new Blocklist
		return blocklist;
	}

	private void addToKnownBasesBuilder(StringBuilder refBases,
			ListIterator<MdElt> mdListIter, Block block) {
		int blockLen = block.length();
		byte[] blockReadBases = block.getReadBases();
		int from = 0;
		while(from < blockLen) {
			// Get the next MdElt
			MdElt mdElt = mdListIter.next();
			// Easy case: if empty, we don't have to do anything
			if(mdElt.length == 0) continue;
			
			// Determine the block read bases index at the end of the block
			int to = from + mdElt.length;
			
			// If it extends past the block, then we'll only use the portion
			// that overlaps the block, and the rest of the MdElt is put back
			// into the list.
			if(to > blockLen) {
				to = blockLen;
				mdListIter.add(mdElt.subMdElt(to - from));
				mdListIter.previous();
			}
			// Otherwise, use the whole thing
			else {
				to = from + mdElt.length;
			}
			
			// Add sequence
			switch(mdElt.type) {
			case MATCH:
				refBases.append(new String(Arrays.copyOfRange(
						blockReadBases, from, to)));
				break;
				
			case MISMATCH:
			case DELETION:
				refBases.append(new String(mdElt.bases));
			}
			
			// Slide from up to to
			from = to;
		}
	}

	private void addKnownBasesToMdFetcher(StringBuilder refBases,
			long refBasesLeft) {
		if(refBases.length() == 0) return;
		curMdFetcher.addKnownBases(refBasesLeft,
				refBases.toString().getBytes());
	}

	/**
	 * Convenience method to set the MD tag on a record based on reference
	 * sequence.
	 * 
	 * @param record the record whose MD to set
	 * @return the blocklist that was created in the process
	 * @throws IllegalStateException if there is no seq db, as it is required
	 */
	public Blocklist fixMd(SAMRecord record) throws IllegalStateException {
		// Assert refFetcher existence
		if(refFetcher == null) throw new IllegalStateException(
				"You cannot fix MD without a sequence file");
		
		// Run with useMd temporarily disabled
		boolean oldUseMd = mdParsingEnabled;
		this.mdParsingEnabled = false;
		Blocklist blocklist;
		try {
			blocklist = newBlocklist(record);
			String md = MdElts.listToString(blocklist.toMdEltList());
			record.setAttribute(SAMTag.MD.name(), md);
		}
		// Restore old useMd
		finally {
			this.mdParsingEnabled = oldUseMd;
		}
		return blocklist;
	}
}
