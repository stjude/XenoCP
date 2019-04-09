package org.stjude.compbio.sam;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.stjude.compbio.util.CountingMap;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.ProgressLoggerInterface;

/**
 * SAMFileWriter wrapper implementation that "normalizes" paired-end reads by
 * buffering writes for each fragment until there is a read 1 and read 2 to
 * write.
 * 
 * By contract, it will always write 0 or 2 reads for each fragment, and it will
 * always write a read 1 and a read 2 for each pair that it writes.  It can
 * handle as input:
 * - "middle" reads (with both first/1, and last/2 flags set, if written, it
 *   will have a bit unset so that it is either first or last only)
 * - non-primary reads (treated the same as primary)
 * - more than 2 reads per fragment (only 2 are written)
 * - missing read 1 or 2 (not written at all)
 * - missing paired flag but 1 or 2 is set (paired flag gets set)
 * 
 * WARNING: this implementation can use a lot of memory.
 */
public class PairNormalizingSAMFileWriter implements SAMFileWriter {
	public static final String WRITE_CTR = "records written";
	public static final String SKIP_READ_CTR = "reads skipped in written fragments";
	public static final String SKIP_FRAG_READ_CTR = "reads in skipped fragments";

	/**
	 * ProgressLogger
	 */
	private ProgressLoggerInterface progressLogger;
	
	/**
	 * Counters
	 */
	private CountingMap<String> counters;
	/**
	 * Underlying writer
	 */
	private SAMFileWriter writer;
	/**
	 * Buffer
	 */
	private Map<String, BufferEntry> buffer;
	/**
	 * Total records sent to addAlignment
	 */
	private int total;
	
	public class BufferEntry {
		private List<SAMRecord> records;
		private int flags;
		private boolean used = false;
		public BufferEntry(SAMRecord record) {
			records = new ArrayList<SAMRecord>(2);
			add(record);
		}
		public void add(SAMRecord record) {
			if(used) return;
			records.add(record);
			flags = flags | record.getFlags();
		}
		public boolean isCompletedBy(SAMRecord record) {
			return ((flags | record.getFlags()) & 0xC0) == 0xC0;
		}
		public void markUsed() {
			used = true;
			records = null;
		}
	}
	
	/**
	 * Constructs a pair normalizing writer using an underlying writer
	 */
	public PairNormalizingSAMFileWriter(SAMFileWriter writer) {
		this.counters = CountingMap.newInitialized(0,
				WRITE_CTR, SKIP_READ_CTR, SKIP_FRAG_READ_CTR);
		this.writer = writer;
		this.buffer = new HashMap<String, BufferEntry>();
		this.total = 0;
	}

	public CountingMap<String> getCounters() {
		return counters;
	}

	public void addAlignment(SAMRecord record) {
		++total;
		BufferEntry entry = buffer.get(record.getReadName());
		if(entry == null) {
			buffer.put(record.getReadName(), new BufferEntry(record));
		}
		else if(entry.used) {
			counters.increment(SKIP_READ_CTR);
		}
		else if(entry.isCompletedBy(record)) {
			SAMRecord[] pair = { entry.records.get(0), record };
			int first = (pair[0].getFirstOfPairFlag() &&
					pair[1].getSecondOfPairFlag()) ? 0 : 1; 
			pair[first].setFlags(pair[first].getFlags() & 0xFF3E | 0x41);
			pair[1-first].setFlags(pair[1-first].getFlags() & 0xFF3E | 0x81);
			writer.addAlignment(pair[0]);
			writer.addAlignment(pair[1]);
			progressLogger.record(pair);
			counters.increment(WRITE_CTR, 2);
			counters.increment(SKIP_READ_CTR, entry.records.size() - 1);
			entry.markUsed();
		}
		else {
			entry.add(record);
		}
	}

	public void close() {
		counters.put(SKIP_FRAG_READ_CTR,
				total
				- counters.getCount(WRITE_CTR)
				- counters.getCount(SKIP_READ_CTR));
		writer.close();
	}

	public SAMFileHeader getFileHeader() {
		return writer.getFileHeader();
	}

	public void setProgressLogger(ProgressLoggerInterface progress) {
		this.progressLogger = progress;
	}

}
