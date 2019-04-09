package org.stjude.compbio.sam.cigar;

import java.util.Iterator;
import java.util.NoSuchElementException;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

/**
 * Iterates over CigarElements, combining all adjacent alignment (M, =, X)
 * elements into single M elements.
 * 
 */
public class CoarseCigarElementIterator implements Iterator<CigarElement> {
	/**
	 * Source iterator
	 */
	private Iterator<CigarElement> iterator;
	/**
	 * Buffer holding next element (have to read ahead one element)
	 */
	private CigarElement buffer;
	
	public CoarseCigarElementIterator(Iterator<CigarElement> iterator) {
		this.iterator = iterator;
		fillBuffer();
	}

	private void fillBuffer() {
		buffer = iterator.hasNext() ? iterator.next() : null;
	}

	public boolean hasNext() {
		return buffer != null;
	}

	public CigarElement next() {
		if(!hasNext()) throw new NoSuchElementException();

		// The output will be the buffered element or a modified version of it
		CigarElement out = buffer;
		
		// Only join alignment elements.  If it's not an alignment element, then
		// just return it
		if(!isAligned(buffer)) {
			fillBuffer();
			return out;
		}
		
		// Join as long as we get alignment elements
		fillBuffer();
		while(buffer != null && isAligned(buffer)) {
			out = new CigarElement(
					out.getLength() + buffer.getLength(),
					CigarOperator.M);
			fillBuffer();
		}
		
		// If we did not join and the original was an = or X, then convert to M
		if(out.getOperator() != CigarOperator.M) {
			out = new CigarElement(out.getLength(), CigarOperator.M);
		}
		
		return out;
	}
	
	/**
	 * Returns true if CigarElement is any of the aligned types
	 */
	private static boolean isAligned(CigarElement elt) {
		CigarOperator op = elt.getOperator();
		return op.consumesReadBases() && op.consumesReferenceBases();
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}

}
