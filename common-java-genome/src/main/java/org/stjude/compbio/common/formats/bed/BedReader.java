package org.stjude.compbio.common.formats.bed;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.stjude.compbio.util.io.IoFormat;
import org.stjude.compbio.util.io.IoRecord;
import org.stjude.compbio.util.io.TabularTextReader;
import org.stjude.compbio.util.io.TabularTextReader.Header;

public class BedReader implements Iterable<BedRecord>, Closeable {
	/**
	 * Delimiter for this format
	 */
	private static final String DELIM = "\t";
	/**
	 * The TabularTextReader.Header object for this format.
	 */
	private static final Header HEADER = Header.forbidden();
	
	/**
	 * Underlying reader
	 */
	private TabularTextReader<BedField> reader;
	
	
	public BedReader(BufferedReader reader) throws IOException {
		this(new TabularTextReader<BedField>(reader, DELIM, HEADER));
	}

	public BedReader(File file) throws IOException {
		this(new TabularTextReader<BedField>(file, DELIM, HEADER));
	}

	public BedReader(InputStream stream) throws IOException {
		this(new TabularTextReader<BedField>(stream, DELIM, HEADER));
	}

	public BedReader(TabularTextReader<BedField> reader) {
		this.reader = reader;
		init();
	}

	/**
	 * Performs general initialization that is applicable regardless of the
	 * constructor used.
	 */
	private void init() {
		reader.addKeys(Arrays.asList(BedField.values()), true);
	}

	public Iterator<BedRecord> iterator() {
		final Iterator<IoRecord<BedField>> iterator = reader.iterator();
		final IoFormat<BedField> ioFormat = reader.getIoFormat();
		
		return new Iterator<BedRecord>() {
			public boolean hasNext() {
				return iterator.hasNext();
			}

			public BedRecord next() {
				IoRecord<BedField> record = iterator.next();
				return new BedRecord(ioFormat, record.getData());
			}
			
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}

	public void close() throws IOException {
		reader.close();
	}
	
	public List<BedRecord> readAll() {
		List<BedRecord> out = new LinkedList<BedRecord>();
		for(BedRecord record: this) out.add(record);
		return out;
	}
}
