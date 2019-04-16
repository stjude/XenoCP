package org.stjude.compbio.util.io;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import org.stjude.compbio.util.DelimitedListBuilder;

/**
 * A writer for tabular text that writes IoRecords.
 * 
 * Note that this exposes a writeLine() method, so you can also print raw lines.
 * 
 * The header may optionally be written immediately at construction time using
 * the ioFormat.
 */
public class TabularTextWriter<K> implements Closeable {
	/**
	 * Underlying writer
	 */
	private BufferedWriter writer;
	/**
	 * Delimited list builder.
	 * 
	 * Every operation is done a whole line at a time, so we can store a DLB as
	 * long as it is cleared before every use.
	 */
	private DelimitedListBuilder dlb;
	/**
	 * IoFormat
	 */
	private IoFormat<K> ioFormat;
	
	public TabularTextWriter(BufferedWriter writer,
			String delim, IoFormat<K> ioFormat, boolean writeHeader)
	throws IOException {
		this.writer = writer;
		this.dlb = new DelimitedListBuilder(delim);
		this.ioFormat = ioFormat;
		if(writeHeader) writeHeader();
	}
	public TabularTextWriter(File file,
			String delim, IoFormat<K> ioFormat, boolean writeHeader)
	throws IOException {
		this(new BufferedWriter(new FileWriter(file)),
				delim, ioFormat, writeHeader);
	}
	public TabularTextWriter(OutputStream stream,
			String delim, IoFormat<K> ioFormat, boolean writeHeader)
	throws IOException {
		this(new BufferedWriter(new OutputStreamWriter(stream)),
				delim, ioFormat, writeHeader);
	}

	public IoFormat<K> getIoFormat() {
		return ioFormat;
	}
	
	/**
	 * Writes the header as determined by the ioFormat
	 * @throws IOException 
	 */
	public void writeHeader() throws IOException {
		dlb.clear();
		for(IoField<K> field: ioFormat.getFields()) {
			dlb.append(field.getName());
		}
		writer.write(dlb.toString() + "\n");
	}
	
	/**
	 * Writes a record.
	 * 
	 * @throws IOException
	 */
	public void writeIoRecord(IoRecord<K> ioRecord) throws IOException {
		dlb.clear();
		// Need to loop rather than just appending ioRecord.values() because
		// that will make sure that unset values get an empty column
		for(IoField<K> field: ioFormat.getFields()) {
			dlb.append(ioRecord.valueForField(field));
		}
		writer.write(dlb.toString() + "\n");
	}
	
	/**
	 * Writes a raw line
	 * @throws IOException 
	 */
	public void writeLine(String line) throws IOException {
		writer.write(line + "\n");
	}
	
	/**
	 * Closes the writer
	 * @throws IOException 
	 */
	public void close() throws IOException {
		writer.close();
	}
}
