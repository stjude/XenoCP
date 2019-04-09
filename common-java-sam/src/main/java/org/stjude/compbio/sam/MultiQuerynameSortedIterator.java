package org.stjude.compbio.sam;

import java.util.ArrayList;
import java.util.List;
import java.util.PriorityQueue;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;

/**
 * Returns records from multiple SAMRecordIterators/SamReaders that are all
 * sorted in queryname order.
 * 
 * You may fetch records one at a time or you may get all records for a
 * read at a time.
 * 
 */
public class MultiQuerynameSortedIterator implements SAMRecordIterator {
	/**
	 * Underlying iterators
	 */
	private QuerynameSortedIterator[] iterators;
	/**
	 * Whether to parse SE read names to cast to PE read ID
	 */
	private boolean parseSeToPeEnabled;
	/**
	 * Priority queue holding one record from each iterator
	 */
	private PriorityQueue<Record> queue;
	/**
	 * Buffered record from off the head of the queue
	 */
	private Record nextRecord;
	/**
	 * ReadId of nextRecord
	 */
	private ReadId nextReadId;

	/**
	 * Constructor using the underlying iterators
	 * 
	 * @param iterators QuerynameSortedIterators
	 */
	public MultiQuerynameSortedIterator(QuerynameSortedIterator[] iterators) {
		this(iterators, false);
	}

	/**
	 * Constructor using the underlying iterators, optionally with parsing of
	 * SE read names to cast to PE read IDs enabled.
	 * 
	 * @param iterators QuerynameSortedIterators
	 * @param parseSeToPeEnabled whether to enable SE-to-PE casting
	 */
	public MultiQuerynameSortedIterator(QuerynameSortedIterator[] iterators,
			boolean parseSeToPeEnabled) {
		this.iterators = iterators;
		this.parseSeToPeEnabled = parseSeToPeEnabled;
		initQueue();
	}
	
	/**
	 * Constructor using SamReaders
	 * 
	 * @param readers SamReaders
	 */
	public MultiQuerynameSortedIterator(SamReader[] readers) {
		this(readers, false);
	}

	/**
	 * Constructor using SamReaders, optionally with parsing of
	 * SE read names to cast to PE read IDs enabled.
	 * 
	 * @param readers SamReaders
	 * @param parseSeToPeEnabled whether to enable SE-to-PE casting
	 */
	public MultiQuerynameSortedIterator(SamReader[] readers,
			boolean parseSeToPeEnabled) {
		this.iterators = new QuerynameSortedIterator[readers.length];
		for(int i = 0; i < readers.length; i++) {
			this.iterators[i] = new QuerynameSortedIterator(readers[i]);
		}
		this.parseSeToPeEnabled = parseSeToPeEnabled;
		initQueue();
	}
	
	/**
	 * Constructs queue and reads in two records from each iterator; also fills
	 * the one-record buffer
	 */
	private void initQueue() {
		queue = new PriorityQueue<Record>();
		for(int i = 0; i < iterators.length; i++) {
			readInput(i);
			readInput(i);
		}
		poll();
	}
	
	/**
	 * Reads one record from iterator i into the queue.
	 * 
	 * @param i iterator index
	 * @return boolean true if a record could be read, false otherwise
	 */
	private boolean readInput(int i) {
		if(iterators[i].hasNext()) {
			queue.add(new Record(iterators[i].next(), i));
			return true;
		}
		else {
			return false;
		}
	}
	
	/**
	 * Pulls head of the queue into the buffer, and pulls one record from the
	 * previous head record's iterator into the queue; returns the previously
	 * buffered record.
	 * 
	 * @return previously buffered record
	 */
	private Record poll() {
		Record out = nextRecord;
		nextRecord = queue.poll();
		if(nextRecord == null) {
			nextReadId = null;
		}
		else {
			nextReadId = parseSeToPeEnabled
					? ReadId.newInstanceFlex(nextRecord.getSamRecord())
					: new ReadId(nextRecord.getSamRecord());
			readInput(nextRecord.getIndex());
		}
		return out;
	}

	/**
	 * A record, which holds a SAMRecord and records the source iterator's index
	 */
	public static class Record implements Comparable<Record> {
		private SAMRecord samRecord;
		private int index;
		public Record(SAMRecord record, int index) {
			this.samRecord = record;
			this.index = index;
		}
		public SAMRecord getSamRecord() {
			return samRecord;
		}
		public int getIndex() {
			return index;
		}
		public int compareTo(Record other) {
			// First compare the read names
			int cmp = samRecord.getReadName().compareTo(
					other.samRecord.getReadName());
			if(cmp != 0) return cmp;
			
			// Then, unpaired comes before read1 before read2
			if(samRecord.getReadPairedFlag() !=
					other.samRecord.getReadPairedFlag()) {
				return samRecord.getReadPairedFlag() ? 1 : -1;
			}
			// For paired, 1 comes before neither comes before 2
			else if(samRecord.getReadPairedFlag()) {
				cmp = pairingToInt(samRecord) - pairingToInt(other.samRecord);
				if(cmp != 0) return cmp;
			}
			
			// Next compare the indexes
			return index - other.index;
		}
		
		private int pairingToInt(SAMRecord record) {
			if(record.getFirstOfPairFlag()) return 1;
			else if(record.getSecondOfPairFlag()) return 3;
			else return 2;
		}
	}

	// SAMRecordIterator methods
	
	public void close() {
		for(QuerynameSortedIterator iterator: iterators) iterator.close();
	}

	public boolean hasNext() {
		return nextRecord != null;
	}

	public SAMRecord next() {
		Record record = poll();
		return record == null ? null : record.getSamRecord();
	}

	public void remove() {
		throw new UnsupportedOperationException("Remove not allowed");
	}

	public SAMRecordIterator assertSorted(SortOrder sortOrder) {
		assert sortOrder == SortOrder.queryname;
		return this;
	}
	
	// Specialty methods
	
	/**
	 * Returns the next record without removing it
	 */
	public SAMRecord peek() {
		return nextRecord == null ? null : nextRecord.getSamRecord();
	}
	
	/**
	 * Returns the ReadId of the next record without removing it
	 */
	public ReadId peekReadId() {
		return nextReadId;
	}
	
	/**
	 * Returns an array of all of the records for the same read as the next
	 * record.
	 * 
	 * The array will have the same number of elements as iterators, and the
	 * positions in the array correspond to the source iterators.  If there
	 * was no record with the current read id in one or more interators, then
	 * those array elements will be null.
	 * 
	 * @return array of SAMRecords or null if iterators exhausted
	 */
	public SAMRecord[] nextArray() {
		// Get the read id we are going to be reading
		ReadId readId = nextReadId;
		if(readId == null) return null;
		
		// Add the first record
		SAMRecord[] out = new SAMRecord[iterators.length];
		Record record = poll();
		out[record.getIndex()] = record.getSamRecord();
		
		// Continue adding records until the read id changes or we are out of
		// records to add
		while(readId.equals(nextReadId)) {
			record = poll();
			out[record.getIndex()] = record.getSamRecord();
		}
		
		return out;
	}
	
	/**
	 * Returns a list of all of the records for the same read as the next
	 * record.
	 * 
	 * @return list of SAMRecords or null if iterators exhausted
	 */
	public List<SAMRecord> nextList() {
		// Get the read id we are going to be reading
		ReadId readId = nextReadId;
		if(readId == null) return null;
		
		// Add the first record
		List<SAMRecord> out = new ArrayList<SAMRecord>(iterators.length);
		out.add(poll().getSamRecord());
		
		// Continue adding records until the read id changes or we are out of
		// records to add
		while(readId.equals(nextReadId)) {
			out.add(poll().getSamRecord());
		}
		
		return out;
	}
}
