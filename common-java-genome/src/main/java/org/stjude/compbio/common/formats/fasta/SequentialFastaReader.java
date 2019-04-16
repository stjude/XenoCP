package org.stjude.compbio.common.formats.fasta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.regex.Pattern;

import org.stjude.compbio.util.io.FormatException;

/**
 * Reads a FASTA file sequentially
 */
public class SequentialFastaReader implements Iterable<FastaRecord> {
	/**
	 * Defline pattern
	 */
	private static final Pattern DEFLINE_PATTERN = Pattern.compile("^>");
	
	/**
	 * Underlying reader
	 */
	private BufferedReader reader;
	
	/**
	 * Buffered defline
	 */
	private String nextDefline;
	
	public SequentialFastaReader(BufferedReader reader) throws IOException {
		this.reader = reader;
		init();
	}

	public SequentialFastaReader(File file) throws IOException {
		this(new BufferedReader(new FileReader(file)));
	}

	public SequentialFastaReader(InputStream stream) throws IOException {
		this(new BufferedReader(new InputStreamReader(stream)));
	}

	/**
	 * Initializes by filling the buffer
	 * @throws IOException 
	 */
	private void init() throws IOException {
		// Read defline
		this.nextDefline = reader.readLine();
		// If null, we are done
		if(nextDefline == null) return;
		
		// Otherwise, validate
		if(!DEFLINE_PATTERN.matcher(nextDefline).find()) {
			throw new FormatException(
					"Expected defline as first line in FASTA, got: " +
							nextDefline);
		}
	}
	
	/**
	 * Reads the next record as a FastaRecord.
	 * 
	 * Returns null when exhausted.
	 * @throws IOException 
	 */
	public FastaRecord readNext() throws IOException {
		if(nextDefline == null) return null;
		String defline = nextDefline;
		String seq = readToNextDefline();
		return new FastaRecord(defline, seq);
	}

	/**
	 * Reads until the next defline is found, returning all sequence read, and
	 * storing the defline in nextDefline.
	 * @throws IOException 
	 */
	private String readToNextDefline() throws IOException {
		StringBuffer seq = new StringBuffer();
		String line = null;
		while((line = reader.readLine()) != null) {
			// Break on defline
			if(DEFLINE_PATTERN.matcher(line).find()) break;
			// Append to seq
			seq.append(line);
		}
		// Set nextDefline and return
		this.nextDefline = line;
		return seq.toString();
	}

	public Iterator<FastaRecord> iterator() {
		return new Iterator<FastaRecord>() {
			public boolean hasNext() {
				return nextDefline != null;
			}

			public FastaRecord next() {
				try {
					return readNext();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
	}
}
