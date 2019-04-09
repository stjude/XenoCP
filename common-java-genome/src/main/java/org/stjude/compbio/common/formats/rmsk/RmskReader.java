package org.stjude.compbio.common.formats.rmsk;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.StdRefInterval;
import org.stjude.compbio.common.interval.index.RefIntervalIndex;
import org.stjude.compbio.util.io.FormatException;

/**
 * Reads UCSC's rmsk file with RepeatMasker regions
 */
public class RmskReader {
	/**
	 * Reader used to read the file
	 */
	private final BufferedReader reader;
	
	/**
	 * Last line read
	 */
	private String line;

	/**
	 * Opens a RefFlatReader using the underlying BufferedReader
	 */
	public RmskReader(BufferedReader reader) {
		this.reader = reader;
	}
	
	/**
	 * Opens a RefFlatReader on a refFlat-formatted file
	 * @throws FileNotFoundException if the file is not found
	 */
	public RmskReader(File file) throws FileNotFoundException {
		this(new BufferedReader(new FileReader(file)));
	}
	
	/**
	 * Opens a RefFlatReader on a refFlat-formatted input stream
	 */
	public RmskReader(InputStream in) {
		this(new BufferedReader(new InputStreamReader(in)));
	}
	
	/**
	 * Reads all into a RefIntervalIndex.
	 * 
	 * Values in the index are the input lines
	 * 
	 * @param stranded whether or not to build a stranded index
	 * @throws IOException on read error
	 * @throws FormatException on input formatting problem
	 */
	public RefIntervalIndex<String> readAllAsRefIntervalIndex(
			boolean stranded)
	throws IOException, FormatException {
		RefIntervalIndex<String> out = new RefIntervalIndex<String>(stranded);
		RefInterval interval;
		while((interval = readAsRefInterval()) != null) {
			out.add(interval, line);
		}
		return out;
	}
	
	/**
	 * Reads the next record as a RefInterval
	 * @return next RefInterval, or null if at end of file
	 * @throws IOException on read error
	 * @throws FormatException if the input is not correctly formatted
	 */
	public RefInterval readAsRefInterval() throws IOException, FormatException {
		String line = readAsString();
		if(line == null) return null;
		String[] parts = line.split("\t");
		try {
			String refName = parts[5];
			long left = Long.parseLong(parts[6]);
			long right = Long.parseLong(parts[7]);
			Strand strand = Strand.valueOfStr(parts[9]);
			return new StdRefInterval(
					refName, strand, left, right);
		} catch(NumberFormatException e) {
			throw new FormatException(e);
		} catch(IllegalArgumentException e) {
			throw new FormatException(e);
		} catch(IndexOutOfBoundsException e) {
			throw new FormatException(e);
		}
	}

	/**
	 * Reads the next data line of the file
	 * @return next data line, or null if at the end of input
	 * @throws IOException on read error
	 */
	public String readAsString() throws IOException {
		line = reader.readLine();
		while(line != null &&
				(line.equals("") || line.startsWith("#"))) {
			line = reader.readLine();
		}
		return line;
	}
	
	/**
	 * Closes the underlying reader
	 * @throws IOException from underlying reader
	 */
	public void close() throws IOException {
		reader.close();
	}
}
