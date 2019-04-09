package org.stjude.compbio.common;

import java.io.File;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.stjude.compbio.util.MainHelper;

public class GenomeMainHelper extends MainHelper {

	/**
	 * Gene annotation option name
	 */
	public static String GENE_ANNO_OPT = "g";
	/**
	 * Reference sequence file option name
	 */
	public static String REF_FASTA_OPT = "f";
	/**
	 * Output FASTA width option name
	 */
	public static String FASTA_WIDTH_OPT = "w";
	/**
	 * Default FASTA width
	 */
	public static int DEFAULT_FASTA_WIDTH = 72;

	public GenomeMainHelper(Options options, String[] args, boolean strict) {
		super(options, args, strict);
	}

	public GenomeMainHelper(Options options) {
		super(options);
	}

	/**
	 * Adds the gene annotation option
	 */
	public void addGeneAnnoOption() {
		getOptions().addOption(new Option(GENE_ANNO_OPT, true,
				"gene annotation file in refFlat format"));
	}

	/**
	 * Gets the user-selected gene annotation file.
	 */
	public File getGeneAnnoFile() {
		return getArgFile(GENE_ANNO_OPT);
	}

	/**
	 * Adds the reference genome FASTA option
	 */
	public void addRefFastaOption() {
		getOptions().addOption(new Option(REF_FASTA_OPT, true,
				"reference genome sequence in FASTA format"));
	}

	/**
	 * Gets the user-selected reference genome FASTA file.
	 */
	public File getRefFastaFile() {
		return getArgFile(REF_FASTA_OPT);
	}

	/**
	 * Adds the FASTA width
	 */
	public void addFastaWidthOption() {
		getOptions().addOption(new Option(FASTA_WIDTH_OPT, true,
				"output FASTA line length (default " + 
						DEFAULT_FASTA_WIDTH + ")"));
	}

	/**
	 * Gets the user-selected output FASTA width
	 */
	public int getFastaWidth() {
		Integer fastaWidth = getArgInt(FASTA_WIDTH_OPT);
		return fastaWidth == null ? DEFAULT_FASTA_WIDTH : fastaWidth;
	}

}