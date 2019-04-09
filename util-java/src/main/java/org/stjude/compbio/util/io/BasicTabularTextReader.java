package org.stjude.compbio.util.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;

/**
 * Tabular text reader with Object keys, usually used when you do not want to
 * use the key functionality but just want to refer to fields by position and/or
 * name from header.
 */
public class BasicTabularTextReader extends TabularTextReader<Object> {

	public BasicTabularTextReader(BufferedReader reader, String delim,
			Header headerOptions) throws IOException {
		super(reader, delim, headerOptions);
	}

	public BasicTabularTextReader(File file, String delim, Header headerOptions)
			throws IOException {
		super(file, delim, headerOptions);
	}

	public BasicTabularTextReader(InputStream stream, String delim,
			Header headerOptions) throws IOException {
		super(stream, delim, headerOptions);
	}

}
