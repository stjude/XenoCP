package org.stjude.compbio.sam;

import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

/**
 * Methods for working with Cigar
 */
public class CigarTools {

	/**
	 * Cleans up a Cigar and returns a new one.
	 * 
	 * Right now there is one cleaning operation: neighboring Is and Ds are
	 * combined together into an M (mismatch) element and at most 1 I or D
	 * element.
	 */
	public static Cigar clean(Cigar orig) {
		// Create the new empty Cigar
		List<CigarElement> cleanElts = new LinkedList<CigarElement>();
		
		// Total length M, I, D that has been read but not yet written.  These
		// will always represent the last few elements read, and the M amount
		// will need to precede the I/D amount.
		int mmLen = 0, insLen = 0, delLen = 0;
		// Pending I/D elements only (the preceding M is not represented)
		LinkedList<CigarElement> pending = new LinkedList<CigarElement>();
		for(CigarElement cigarElt: orig.getCigarElements()) {
			// Handle I/D coalescence first
			switch(cigarElt.getOperator()) {
			case I:
				insLen += cigarElt.getLength();
				pending.add(cigarElt);
				break;
				
			case D:
				delLen += cigarElt.getLength();
				pending.add(cigarElt);
				break;
				
			default:
				// If we have anything other than an I or a D and we have a 
				// pending run of Is/Ds, then we need to do a collapse
				if(!pending.isEmpty()) {
					// If there is preceding M, then push it onto beginning of
					// list
					if(mmLen > 0) {
						pending.addFirst(
								new CigarElement(mmLen, CigarOperator.M));
					}
					// Collapse the list
					List<CigarElement> collapsed =
						collapse(pending, insLen, delLen);
					// The collapsed list will be of the form "M?[ID]?"; that
					// is, it will have 0-1 M elements followed by an I or a D
					// or nothing.  It will not be empty, though.
					// If it has an M only, then do not add the M yet, but
					// rather, push it into pending
					if(collapsed.size() == 1 &&
							collapsed.get(0).getOperator() == CigarOperator.M) {
						mmLen = collapsed.get(0).getLength();
					}
					// Otherwise, add the collapsed elements, and there is no
					// pending M
					else {
						cleanElts.addAll(collapsed);
						mmLen = 0;
					}
					// Clear pending I/D list and I/D counts
					pending.clear();
					insLen = 0;
					delLen = 0;
				}
			}
			
			// Now, handle the current non-I/D element
			switch(cigarElt.getOperator()) {
			case I:
			case D:
				// Nothing to do, we handled it already
				break;
			case M:
				// For M, do not add it, just accumulate the total M length
				// that is pending
				mmLen += cigarElt.getLength();
				break;
			default:
				// For anything other than I/D/M, write any pending M element,
				// and then write the element itself
				if(mmLen > 0) {
					cleanElts.add(
							new CigarElement(mmLen, CigarOperator.M));
					mmLen = 0;
				}
				cleanElts.add(cigarElt);
			}
		}
		
		// Do final collapse/addition if necessary
		CigarElement pendingM = null;
		if(mmLen > 0) pendingM = new CigarElement(mmLen, CigarOperator.M);
		// If pending I/D then collapse
		if(!pending.isEmpty()) {
			if(pendingM != null) pending.addFirst(pendingM);
			cleanElts.addAll(collapse(pending, insLen, delLen));
		}
		// Otherwise, if pending M then just add it
		else if(pendingM != null) {
			cleanElts.add(pendingM);
		}
		
		return new Cigar(cleanElts);
	}

	/**
	 * Collapses multiple neighboring INS and DEL elements.
	 * 
	 * @param pending the list of pending elements
	 * @param insLen total I length
	 * @param delLen total D length
	 * @param mmLen true if the M element should come first
	 * @return list of collapsed elements
	 */
	private static List<CigarElement> collapse(
			List<CigarElement> pending, int insLen, int delLen) {
		// Easy case: no collapse if only 1 element
		if(pending.size() < 2) return pending;
		
		// Do the collapse
		List<CigarElement> collapsed = new LinkedList<CigarElement>();
		
		// For every base of both ins and del, convert to M
		int mmLen = Math.min(insLen, delLen);
		insLen -= mmLen;
		delLen -= mmLen;
		
		// If there is pending M, then add that to the M length to add
		if(pending.get(0).getOperator() == CigarOperator.M) {
			mmLen += pending.get(0).getLength();
		}
		
		// Make and add the M element if there is one
		if(mmLen > 0) {
			collapsed.add(new CigarElement(mmLen, CigarOperator.M));
		}
		
		// Add the gap (I or D) element of there is leftover I or D
		CigarElement gapElt = null;
		if(insLen > 0) gapElt = new CigarElement(insLen, CigarOperator.I);
		else if(delLen > 0) gapElt = new CigarElement(delLen, CigarOperator.D);
		// Do the add if there is one to add
		if(gapElt != null) collapsed.add(gapElt);
		
		return collapsed;
	}

}
