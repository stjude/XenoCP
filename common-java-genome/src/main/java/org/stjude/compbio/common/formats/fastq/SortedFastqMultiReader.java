package org.stjude.compbio.common.formats.fastq;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.fastq.FastqRecord;

/**
 * Class with functionality for reading multiple FASTQ files and returning
 * records with matching names.
 * 
 * This is useful for reading paired files where some sequences are missing
 * from one file or the other.  By default, only sequences that are present
 * in all files are returned, but you can specify that one or more inputs are
 * optional.  In that case, nulls may be returned in the output for those if
 * there is no match.
 * 
 * The files must have their records sorted by sequence name.
 * 
 */
public class SortedFastqMultiReader implements Iterable<FastqRecord[]> {
	/**
	 * The underlying readers' iterators
	 */
	private final Iterator<FastqRecord>[] iterators;
	/**
	 * The iterator indexes that are optional
	 */
	private final Set<Integer> optIndexes;
	/**
	 * Buffered next result
	 */
	private FastqRecord[] next;
	/**
	 * Pending results from optional iterators
	 */
	private FastqRecord[] pending;

	/**
	 * Constructs a SortedFastqMultiReader from underlying iterators
	 * 
	 * @param iterators iterators
	 * @param optIndexes indexes of optional iterators (or null for all req'd)
	 */
	public SortedFastqMultiReader(Iterator<FastqRecord>[] iterators,
			Collection<Integer> optIndexes) {
		this.iterators = iterators;
		this.optIndexes = (optIndexes == null) ? new HashSet<Integer>() :
			new HashSet<Integer>(optIndexes);
		if(this.optIndexes.size() >= iterators.length) {
			throw new IllegalArgumentException("Too many optional indexes");
		}
		this.next = new FastqRecord[iterators.length];
		this.pending = new FastqRecord[iterators.length];
		loadBuffer();
	}

	/**
	 * Constructs a SortedFastqMultiReader from readers
	 * 
	 * @param readers readers
	 * @param optIndexes indexes of optional readers (or null for all req'd)
	 */
	public SortedFastqMultiReader(
			List<? extends Iterable<FastqRecord>> readers,
			Collection<Integer> optIndexes) {
		this(toIterators(readers), optIndexes);
	}

	/**
	 * Constructs a SortedFastqMultiReader from readers with all iterators
	 * required.
	 * 
	 * @param readers readers
	 */
	public SortedFastqMultiReader(
			List<? extends Iterable<FastqRecord>> readers) {
		this(toIterators(readers), null);
	}

	private static Iterator<FastqRecord>[] toIterators(
			List<? extends Iterable<FastqRecord>> readers) {
		@SuppressWarnings("unchecked")
		Iterator<FastqRecord>[] iterators =	new Iterator[readers.size()];
		int i = 0;
		for(Iterable<FastqRecord> reader: readers) {
			iterators[i++] = reader.iterator();
		}
		return iterators;
	}

	private void loadBuffer() {
		// Nothing to do if already exhausted
		if(next == null) return;
		
		// Loop in a circular fashion until all at same name
		String maxName = null;
		int maxSetAt = -1;
		for(int i = 0; i != maxSetAt; i = (i + 1) % iterators.length) {
			// Skip optional on first pass
			if(optIndexes.contains(i)) continue;
			// Read until max is hit or a new one can be set
			while(true) {
				// Check for exhaustion
				if(!iterators[i].hasNext()) {
					next = null;
					return;
				}
				
				// Get the next record and process the name
				next[i] = iterators[i].next();
				String name = FastqHelper.getName(next[i]);
				name = trimSlashReadNumber(name);
				if(maxName == null || name.compareTo(maxName) > 0) {
					maxName = name;
					maxSetAt = i;
					break;
				}
				else if(maxName.equals(name)) {
					break;
				}
			}
		}
		
		// Fill in optional
		for(int i: optIndexes) {
			// Read until max is hit
			while(true) {
				// Get next record from pending or iterator
				if(pending[i] != null) {
					next[i] = pending[i];
					pending[i] = null;
				}
				else if(iterators[i].hasNext()) {
					next[i] = iterators[i].next();
				}
				else {
					next[i] = null;
					break;
				}
				
				// Process name
				String name = next[i].getReadHeader().split(" ")[0];
				if(name.compareTo(maxName) > 0) {
					pending[i] = next[i];
					next[i] = null;
					break;
				}
				else if(maxName.equals(name)) {
					break;
				}
			}
		}
	}

	public Iterator<FastqRecord[]> iterator() {
		return new Iterator<FastqRecord[]>() {

			public boolean hasNext() {
				return next != null;
			}

			public FastqRecord[] next() {
				FastqRecord[] out = trimAllSlashReadNumbers(next);
				loadBuffer();
				return out;
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
			
		};
	}

    /**
     * Trim all /1 or /2 for FastqRecords
     * 
     * @param FastqRecord[] out
	 * @return Trimmed FastqRecord[]
     */
	private static FastqRecord[] trimAllSlashReadNumbers(FastqRecord[] out) {
		FastqRecord[] processedFastqRecord = new FastqRecord[out.length];
		for(int i = 0; i < out.length; i++){
			String name = out[i].getReadHeader();
			name = trimSlashReadNumber(name);
			processedFastqRecord[i] = new FastqRecord(name, out[i].getReadString(), out[i].getBaseQualityHeader(), out[i].getBaseQualityString());
		}
		return processedFastqRecord;
	}

    /**
     * Trim /1 or /2 for a given read name.
     *
     * @param String name A read name
     * @return Trimmed read name
     */	
	private static String trimSlashReadNumber(String name) {
		if(name.charAt(name.length()-2) == '/'){
			if(name.charAt(name.length()-1) == '1' || name.charAt(name.length()-1) == '2'){
				name = name.substring(0, name.length()-2);
			}
		}
		return name;		
	}

}
