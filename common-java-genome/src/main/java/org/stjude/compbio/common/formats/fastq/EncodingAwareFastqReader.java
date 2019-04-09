package org.stjude.compbio.common.formats.fastq;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.LinkedList;

import org.stjude.compbio.util.io.FormatException;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.FastqQualityFormat;

/**
 * Subclass of FastqReader that can auto-detect encoding
 */
public class EncodingAwareFastqReader extends FastqReader {
	/**
	 * Buffer of records that had been read for determining encoding
	 */
	private LinkedList<FastqRecord> buffer;
	/**
	 * Determined format (lazily initialized)
	 */
	private FastqQualityFormat encoding = null;
	
	/**
	 * Opens an EncodingAwareFastqReader on a file
	 *  
	 * @param file the FASTQ file
	 * @throws IOException on read error
	 * @throws FormatException if the format of the first record is bad
	 */
	public EncodingAwareFastqReader(File file)
	throws IOException, FormatException {
		super(file);
		this.buffer = new LinkedList<FastqRecord>();
	}

	/**
	 * Constructs an EncodingAwareFastqReader on the given reader, with the
	 * given number of lines already read from the reader (for reporting)
	 * 
	 * @param reader underlying buffered reader
	 * @param linesRead number of lines already read from it
	 * @throws IOException on read error
	 * @throws FormatException if the format of the first record is bad
	 */
	public EncodingAwareFastqReader(BufferedReader reader, int linesRead)
	throws IOException, FormatException {
		super(reader, linesRead);
		this.buffer = new LinkedList<FastqRecord>();
	}

	@Override
	protected boolean hasNextRecord() {
		return !buffer.isEmpty() || super.hasNextRecord();
	}

	@Override
	public FastqRecord readRecord() throws IOException, FormatException {
		FastqRecord next =
				buffer.isEmpty() ? super.readRecord() : buffer.poll();
		if(encoding == null) determineEncoding(next);
		return next;
	}

	/**
	 * Sets the quality encoding format (there will be no auto-detect in the
	 * future).
	 * 
	 * @param encoding the encoding to set
	 */
	public void setEncoding(FastqQualityFormat encoding) {
		this.encoding = encoding;
	}
	
	/**
	 * Returns quality encoding format, detecting if necessary.
	 * 
	 * If detecting is necessary, it will read ahead until it can determine
	 * the encoding, and if there are problems later in the file, then it will
	 * throw an exception
	 * 
	 * @return FastqQualityFormat the determined encoding, null if unknown
	 * @throws FormatException if a formatting problem is encountered
	 * @throws IOException on read error
	 */
	public FastqQualityFormat getEncoding()
	throws IOException, FormatException {
		while(encoding == null && super.hasNextRecord()) {
			FastqRecord next = super.readRecord();
			determineEncoding(next);
			buffer.add(next);
		}
		return encoding;
	}
	
	/**
	 * Sets encoding according to the given record
	 * 
	 * @param record fastq record to use
	 */
	private void determineEncoding(FastqRecord record) {
		byte[] qualBytes = record.getBaseQualityString().getBytes();
		for(byte qualByte: qualBytes) {
			if(qualByte < 59) {
				encoding = FastqQualityFormat.Standard;
				return;
			}
			else if(qualByte > 75) {
				encoding = FastqQualityFormat.Illumina;
				return;
			}
		}
	}
}
