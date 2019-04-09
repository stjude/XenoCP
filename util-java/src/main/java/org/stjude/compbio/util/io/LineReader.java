package org.stjude.compbio.util.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;

/**
 * Simple class for reading from an underlying reader one line at a time, with
 * support for peeking ahead one line.
 * 
 * Also provides an iterable interface.
 * 
 * When extending this class, if you do not want the Iterable<String>
 * functionality, then extends the base class, AbstractLineReader, which
 * contains all of the functionality except for iterable.  This will allow you
 * to make your class iterable over something else.
 */
public class LineReader extends AbstractLineReader implements Iterable<String> {
	/**
	 * Construct using a BufferedReader
	 */
	public LineReader(BufferedReader reader) {
		super(reader);
	}

	/**
	 * Construct using a File
	 * @throws FileNotFoundException if the file is not found
	 */
	public LineReader(File file) throws FileNotFoundException {
		super(file);
	}

	/**
	 * Construct using an InputStream
	 */
	public LineReader(InputStream stream) {
		super(stream);
	}

	public Iterator<String> iterator() {
		return new AbstractIterator<String>() {
			public String next() {
				try {
					return readLine();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		};
	}
}
