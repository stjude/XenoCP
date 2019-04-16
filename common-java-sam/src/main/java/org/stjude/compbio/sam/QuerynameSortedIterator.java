package org.stjude.compbio.sam;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;

/**
 * SAMRecordIterator for queryname-sorted data with peeking and other helper
 * functionality
 *
 */
public class QuerynameSortedIterator implements SAMRecordIterator {
	/**
	 * The basis iterator
	 */
	private SAMRecordIterator iterator;
	/**
	 * Buffered next record
	 */
	private SAMRecord next;
	
	/**
	 * Constructs from an underlying SAMRecordIterator
	 * 
	 * @param iterator iterator instance to use; will be modified 
	 */
	public QuerynameSortedIterator(SAMRecordIterator iterator) {
		this.iterator = iterator;
		iterator.assertSorted(SortOrder.queryname);
		fillBuffer();
	}
	
	/**
	 * Constructs from a SamReader
	 * 
	 * @param reader
	 */
	public QuerynameSortedIterator(SamReader reader) {
		this(reader.iterator());
	}
	
	/**
	 * Reads the next record from the iterator into the buffer
	 */
	private void fillBuffer() {
		next = iterator.hasNext() ? iterator.next() : null;
	}

	/**
	 * Returns an array with the next record and, if it has the same name, the
	 * record after next.
	 * 
	 * The returned array is always 2 elements, though one may be null.  If the
	 * next read is unpaired, then it is returned in element 0 and element 1 is
	 * null.  If the next read is paired, then it is returned in element 0 if
	 * it is the first read in the pair or element 1 if it is the second; its
	 * mate, if present, is returned in the other array element.
	 * 
	 * If the iterator is exhausted, then both elements are null
	 * 
	 * @return the next one or two records
	 */
	public SAMRecord[] nextPair() {
		SAMRecord a = next();
		SAMRecord b = null;
		// Check for paired
		if(a != null && a.getReadPairedFlag()) {
			// Check for mate
			if(next != null && next.getReadName().equals(a.getReadName())) {
				b = next();
			}
			// Return in the correct order
			if(a.getFirstOfPairFlag()) {
				return new SAMRecord[] { a, b };
			}
			else {
				return new SAMRecord[] { b, a };
			}
		}
		// Unpaired: put in first position; second position is null
		else {
			return new SAMRecord[] { a, null };
		}
	}
	
	/**
	 * Returns an array with the next record and, if it has the same name, the
	 * record after next with the next record in the first array position.
	 * 
	 * This is just like nextPair except that the first record found is placed
	 * in the first array element, and the second record found (if there is a
	 * second record) is placed in the second.  Hence, the ordering has no
	 * significance.  This means that unless the iterator is exhausted, the
	 * first element will always be non-null.
	 * 
	 * @return the next one or two records, with index 0 non-null if possible
	 */
	public SAMRecord[] nextPairUnordered() {
		SAMRecord a = next();
		SAMRecord b = null;
		// Check for paired
		if(a != null && a.getReadPairedFlag()) {
			// Check for mate
			if(next != null && next.getReadName().equals(a.getReadName())) {
				b = next();
			}
			// Return in file order
			return new SAMRecord[] { a, b };
		}
		// Unpaired: put in first position; second position is null
		else {
			return new SAMRecord[] { a, null };
		}
	}
	
	/**
	 * Gets the name of the pending record
	 * 
	 * @return name, or null if no record is pending
	 */
	public String peekName() {
		return next == null ? null : next.getReadName();
	}
	
	// Boring delegate stuff from here on out
	
	
	
	public void close() {
		iterator.close();
	}

	public boolean hasNext() {
		return next != null;
	}

	public SAMRecord next() {
		SAMRecord out = next;
		fillBuffer();
		return out;
	}

	public void remove() {
		throw new UnsupportedOperationException(
				"Due to buffering, cannot remove from QuerynameSortedIterator");
	}

	public SAMRecordIterator assertSorted(SortOrder sortOrder) {
		assert sortOrder == SortOrder.queryname;
		return this;
	}

}
