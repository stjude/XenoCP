package org.stjude.compbio.common.formats.fastq;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import org.stjude.compbio.util.io.FormatException;

import htsjdk.samtools.fastq.FastqConstants;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.IOUtil;

/**
 * FastqReader implementation that supports FASTQs with 0-length sequences and
 * throws FormatExceptions.
 */
public class FastqReader implements Iterable<FastqRecord> {
	/**
	 * Underlying reader
	 */
	private final BufferedReader reader;
	/**
	 * Lines read (also line number of last line read)
	 */
	private int linesRead;
	/**
	 * Buffered record
	 */
	private FastqRecord next;
	
	/**
	 * Opens a FastqReader on a file
	 *  
	 * @param file the FASTQ file
	 * @throws IOException on read error
	 * @throws FormatException if the format of the first record is bad
	 */
	public FastqReader(File file) throws IOException, FormatException {
		this(IOUtil.openFileForBufferedReading(file), 0);
	}

	/**
	 * Constructs a FastqReader on the given reader, with the given number of
	 * lines already read from the reader (for reporting)
	 * 
	 * @param reader underlying buffered reader
	 * @param linesRead number of lines already read from it
	 * @throws IOException on read error
	 * @throws FormatException if the format of the first record is bad
	 */
	public FastqReader(BufferedReader reader, int linesRead)
			throws IOException, FormatException {
		this.reader = reader;
		this.linesRead = linesRead;
		loadBuffer();
	}
	
	/**
	 * Returns true if there is a next record
	 */
	protected boolean hasNextRecord() {
		return next != null;
	}
	
	/**
	 * Reads the next record.
	 * 
	 * This method differs slightly from using the iterable/iterator in that it
	 * throws semantic exceptions, instead of wrapping them in
	 * RuntimeExceptions.
	 * 
	 * @throws FormatException if there is a formatting problem
	 * @throws IOException on read error
	 */
	public FastqRecord readRecord() throws IOException, FormatException {
		FastqRecord out = next;
		loadBuffer();
		return out;
	}

	/**
	 * Reads the next record into the buffer
	 * @throws IOException on read error
	 * @throws FormatException if the format of the next record is bad
	 */
	private void loadBuffer() throws IOException, FormatException {
		// Read the defline
		String seqHeaderLine = reader.readLine();
		if(seqHeaderLine == null) {
			next = null;
			return;
		}
		++linesRead;
		assertPrefix(seqHeaderLine, FastqConstants.SEQUENCE_HEADER);
		
		// Sequence line
		String seq = reader.readLine();
		assertLineRead(seq);
		
		// Qual header line
		String qualHeaderLine = reader.readLine();
		assertLineRead(qualHeaderLine);
		assertPrefix(qualHeaderLine, FastqConstants.QUALITY_HEADER);
		
		// Qual line
		String qual = reader.readLine();
		assertLineRead(qual);
		
		// Add record
		next = new FastqRecord(
				seqHeaderLine.substring(1), seq,
				qualHeaderLine.substring(1), qual);
	}

	private void assertLineRead(String seq) throws FormatException {
		if(seq == null) {
			throw new FormatException("Unexpected EOF after line " + linesRead);
		}
		++linesRead;
	}

	private void assertPrefix(String line, String prefix)
			throws FormatException {
		if(!line.startsWith(prefix)) {
			throw new FormatException(
					"Line must start with " + prefix + " at " + linesRead);
		}
	}

	/**
	 * Gets the number of lines read, which is also the line number of the last
	 * line read.
	 */
	public int getLinesRead() {
		return linesRead;
	}

	public Iterator<FastqRecord> iterator() {
		return new Iterator<FastqRecord>() {
			public boolean hasNext() {
				return hasNextRecord();
			}

			public FastqRecord next() {
				try {
					return readRecord();
				} catch (Exception e) {
					throw new RuntimeException(e);
				}
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}
	
	/**
	 * Closes the underlying reader
	 */
	public void close() throws IOException {
		reader.close();
	}
}
