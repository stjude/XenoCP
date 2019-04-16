package org.stjude.compbio.common.formats.psl;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import org.stjude.compbio.util.io.IoRecord;
import org.stjude.compbio.util.io.TabularTextReader;

public class PslReader extends TabularTextReader<PslField> {
	/**
	 * Delimiter for this format
	 */
	private static final String DELIM = "\t";
	/**
	 * The TabularTextReader.Header object for this format.
	 */
	private static final Header HEADER = Header.skipped(5);
	
	public PslReader(BufferedReader reader) throws IOException {
		super(reader, DELIM, HEADER);
		init();
	}

	public PslReader(File file) throws IOException {
		super(file, DELIM, HEADER);
		init();
	}

	public PslReader(InputStream stream) throws IOException {
		super(stream, DELIM, HEADER);
		init();
	}

	/**
	 * Performs general initialization that is applicable regardless of the
	 * constructor used.
	 */
	private void init() {
		addKeys(Arrays.asList(PslField.values()), true);
	}
	
	/**
	 * Simple scoring function
	 */
	public static int simpleScore(IoRecord<PslField> record) {
		int matches = Integer.parseInt(record.valueForKey(PslField.matches));
		int repMatches = Integer.parseInt(record.valueForKey(PslField.repMatches));
		int misMatches = Integer.parseInt(record.valueForKey(PslField.misMatches));
		int qNumInsert = Integer.parseInt(record.valueForKey(PslField.qNumInsert));
		int tNumInsert = Integer.parseInt(record.valueForKey(PslField.tNumInsert));
		return matches + (repMatches / 2) - misMatches - qNumInsert - tNumInsert;
	}
}
