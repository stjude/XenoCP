package org.stjude.compbio.common.formats.fasta;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * Class to write fasta data, formatted to a specific width
 */
public class FastaWriter {
	public static int NO_WRAP = 0;
	
	/**
	 * The wrapped writer
	 */
	private PrintWriter writer;
	/**
	 * The desired line length (or NO_WRAP)
	 */
	private int width;
	
	/**
	 * Number of bases written to current line
	 */
	private int curLineSize;

	public FastaWriter(PrintWriter writer, int width) {
		this.writer = writer;
		this.width = width;
		this.curLineSize = 0;
	}
	
	public FastaWriter(File file, int width) throws IOException {
		this(new PrintWriter(new FileWriter(file)), width);
	}
	
	/**
	 * Writes a new defline, ending any previous sequence
	 * 
	 * @param defline the defline to write INCLUDING LEADING '>'
	 */
	public void writeDefline(String defline) {
		if(curLineSize > 0) writer.println();
		writer.println(defline);
		curLineSize = 0;
	}
	
	/**
	 * Writes some sequence
	 * 
	 * @param seq sequence to write
	 */
	public void writeSeq(String seq) {
		// Seq len
		int len = seq.length();
		
		// Easy case: no wrapping
		if(width == NO_WRAP) {
			writer.print(seq);
			curLineSize += len;
			return;
		}
		
		// Read start index
		int start = 0;
		
		// Write full lines
		while(len - start >= width - curLineSize) {
			int toWrite = width - curLineSize;
			writer.println(seq.substring(start, start + toWrite));
			start += toWrite;
			curLineSize = 0;
		}
		
		// Write any partial line
		int remaining = len - start;
		if(remaining > 0) {
			writer.print(seq.substring(start));
			curLineSize += remaining;
		}
	}

	/**
	 * Closes the underlying PrintWriter
	 */
	public void close() {
		if(curLineSize > 0) writer.println();
		writer.close();
	}
}
