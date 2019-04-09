package org.stjude.compbio.common.formats.refFlat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.gene.StdGeneHitLibrary;
import org.stjude.compbio.common.gene.StdGeneModelLibrary;
import org.stjude.compbio.common.gene.StdLibraryBuilder;
import org.stjude.compbio.common.gene.TxModel;
import org.stjude.compbio.common.interval.LongInterval;
import org.stjude.compbio.common.interval.StdLongInterval;
import org.stjude.compbio.util.io.FormatException;

/**
 * Class that can read refFlat-formatted files
 */
public class RefFlatReader {
	/**
	 * Number of parts to the parsable line
	 */
	private static final int NUM_PARTS = 11;

	/**
	 * Reader used to read the file
	 */
	private final BufferedReader reader;

	/**
	 * Opens a RefFlatReader using the underlying BufferedReader
	 */
	public RefFlatReader(BufferedReader reader) {
		this.reader = reader;
	}
	
	/**
	 * Opens a RefFlatReader on a refFlat-formatted file
	 * @throws FileNotFoundException if the file is not found
	 */
	public RefFlatReader(File file) throws FileNotFoundException {
		this(new BufferedReader(new FileReader(file)));
	}
	
	/**
	 * Opens a RefFlatReader on a refFlat-formatted input stream
	 */
	public RefFlatReader(InputStream in) {
		this(new BufferedReader(new InputStreamReader(in)));
	}
	
	/**
	 * Reads everything into a new GeneModelLibrary
	 * @throws IOException on read error
	 * @throws FormatException if there is a data formatting problem
	 */
	public StdGeneModelLibrary readAllGeneModels()
	throws IOException, FormatException {
		StdLibraryBuilder builder = new StdLibraryBuilder(true);
		builder.startGeneModelLibrary();
		readIntoBuilder(builder);
		return builder.getGeneModelLibrary();
	}
	
	/**
	 * Reads everything into a new GeneHitLibrary
	 * @throws IOException on read error
	 * @throws FormatException if there is a data formatting problem
	 */
	public StdGeneHitLibrary readAllGeneHits()
	throws IOException, FormatException {
		StdLibraryBuilder builder = new StdLibraryBuilder(true);
		builder.startGeneHitLibrary();
		readIntoBuilder(builder);
		return builder.getGeneHitLibrary();
	}
	
	/**
	 * Reads everything into the library builder.
	 * 
	 * If at any point it is detected that the builder has become degenerate,
	 * then reading will stop.
	 * 
	 * @param builder pre-configured StdLibraryBuilder
	 * @return true if the entire file was read, false otherwise
	 * @throws IOException on read error
	 * @throws FormatException if there is a data formatting problem
	 */
	public boolean readIntoBuilder(StdLibraryBuilder builder)
	throws IOException, FormatException {
		String line;
		while(!builder.isDegenerate() && (line = reader.readLine()) != null) {
			// Skip comments/headers
			if(line.startsWith("#")) continue;
			
			// Split
			String[] parts = line.split("\t");
			if(parts.length < NUM_PARTS) {
				throw new FormatException("Expected " + NUM_PARTS +
						" fields, but only found " + parts.length + " in " +
						line);
			}
			
			// Parse
			int field = 0;
			
			// Easy ones
			String geneName = parts[field++];
			String name = parts[field++];
			String chrom = parts[field++];
			
			// Strand
			String strandStr = parts[field++];
			Strand strand;
			try {
				strand = Strand.valueOfStr(strandStr);
			} catch(Exception e) {
				throw new FormatException("Illegal strand: " + strandStr +
						" in line " + line, e);
			}
			
			// tx and cds intervals
			LongInterval tx = newStdLongInterval(
					parts[field++], parts[field++]);
			LongInterval cds = newStdLongInterval(
					parts[field++], parts[field++]);
			if(cds.length() == 0) cds = null;
			
			// Exons
			int exonCount;
			try {
				exonCount = Integer.parseInt(parts[field++]);
			} catch(NumberFormatException e) {
				throw new FormatException(e);
			}
			List<LongInterval> exons = new ArrayList<LongInterval>(exonCount);
			try {
				String[] lefts = parts[field++].split(",");
				String[] rights = parts[field++].split(",");
				for(int i = 0; i < exonCount; i++) {
					exons.add(newStdLongInterval(lefts[i], rights[i]));
				}
			} catch(ArrayIndexOutOfBoundsException e) {
				throw new FormatException(e);
			}
			
			// Build the TxModel
			TxModel txModel = new TxModel(chrom, strand, exons, cds);
			// Double-check tx
			if(!tx.equals(txModel.getTxInterval())) {
				throw new FormatException("Tx interval " + tx +
						" is inconsistent with exons");
			}
			
			// Add to builder
			builder.addTxModel(geneName, name, txModel);
		}
		return !builder.isDegenerate();
	}

	private LongInterval newStdLongInterval(String leftStr, String rightStr) 
	throws FormatException {
		try {
			return new StdLongInterval(
					Long.parseLong(leftStr), Long.parseLong(rightStr));
		} catch(NumberFormatException e) {
			throw new FormatException(e);
		}
	}
}
