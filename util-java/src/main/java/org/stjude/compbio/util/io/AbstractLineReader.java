package org.stjude.compbio.util.io;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;

/**
 * Abstract base class for LineReader which implements everything except for the
 * iterable functionality.
 */
public abstract class AbstractLineReader implements Closeable {
	/**
	 * Partial iterator implementation
	 * 
	 * @param <T> type of the element being iterated
	 */
	protected abstract class AbstractIterator<T> implements Iterator<T> {
		public boolean hasNext() {
			try {
				return peekLine() != null;
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

		public abstract T next();

		public void remove() {
			throw new UnsupportedOperationException();
		}
	}
	
	/**
	 * The underlying reader
	 */
	protected BufferedReader reader;
	/**
	 * Next-line buffer
	 */
	protected String buffer;

	/**
	 * Construct using a BufferedReader
	 */
	public AbstractLineReader(BufferedReader reader) {
		this.reader = reader;
		this.buffer = null;
	}
	
	/**
	 * Construct using a File
	 * @throws FileNotFoundException if the file is not found
	 */
	public AbstractLineReader(File file) throws FileNotFoundException {
		this(new BufferedReader(new FileReader(file)));
	}
	
	/**
	 * Construct using an InputStream
	 */
	public AbstractLineReader(InputStream stream) {
		this(new BufferedReader(new InputStreamReader(stream)));
	}
	
	/**
	 * Peek ahead at the next line
	 * @throws IOException
	 */
	public String peekLine() throws IOException {
		if(buffer == null) buffer = reader.readLine();
		return buffer;
	}

	/**
	 * Reads the next line
	 * @throws IOException
	 */
	public String readLine() throws IOException {
		if(buffer != null) {
			String out = buffer;
			buffer = null;
			return out;
		}
		else {
			return reader.readLine();
		}
	}
	
	public void close() throws IOException {
		reader.close();
	}

}