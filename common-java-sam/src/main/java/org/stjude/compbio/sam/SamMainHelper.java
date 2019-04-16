package org.stjude.compbio.sam;


import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.FileNotFoundException;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.stjude.compbio.common.GenomeMainHelper;

/**
 * Wrapper class for MainHelper that supports common SAM-specific command-line
 * options.
 */
public class SamMainHelper extends GenomeMainHelper {
	/**
	 * Validation stringency option name
	 */
	public static String VALIDATION_STRINGENCY_OPT = "V";
	/**
	 * Sort order option name
	 */
	public static String SORT_ORDER_OPT = "O";
	/**
	 * Default validation stringency
	 */
	private ValidationStringency defaultVs = ValidationStringency.SILENT;

	public SamMainHelper(Options options) {
		super(options);
	}

	/**
	 * Adds the validation stringency option
	 */
	public void addValidationStringencyOption() {
		getOptions().addOption(new Option(VALIDATION_STRINGENCY_OPT, true,
				"validation stringency: STRICT, LENIENT, or SILENT (default: " +
				defaultVs + ")"));
	}

	/**
	 * Gets the user-selected ValidationStringency (or default).
	 * 
	 * This will not throw a CommandLineException for missing option.
	 */
	public ValidationStringency getValidationStringency() {
		String argStr = null;
		try {
			argStr = getArgString(VALIDATION_STRINGENCY_OPT);
		} catch(CommandLineException e) {
			// Drop through to use default
		}
		if(argStr == null) return defaultVs;
		try {
			return ValidationStringency.valueOf(argStr);
		} catch(IllegalArgumentException e) {
			throw new IllegalArgumentException(
					"Invalid validation stringency option: " + argStr, e); 
		}
	}
	
	/**
	 * Sets the user-selected validation stringency as the SAMReader default
	 */
	public void setDefaultValidationStringency() {
		SamReaderFactory.setDefaultValidationStringency(
				getValidationStringency());
	}
	
	/**
	 * Opens a SamReader using the file from the command-line.
	 * @param argName name of argument holding input file
	 * @throws FileNotFoundException if the input file is not found
	 */
	public SamReader openSamReader(String argName) throws CommandLineException, IllegalStateException, FileNotFoundException {
		return openSamReader(getInputFile(argName));
	}
	
	/**
	 * Opens a SamReader using the options from the command-line.
	 * @param file file to open, or null for stdin
	 */
	public SamReader openSamReader(File file) {
		SamReaderFactory factory = SamReaderFactory.makeDefault()
			.validationStringency(getValidationStringency());
		return file == null 
				? factory.open(SamInputResource.of(System.in))
				: factory.open(file);
	}
	
	/**
	 * Adds the sort order option -O
	 */
	public void addSortOrderOption() {
		getOptions().addOption(new Option(SORT_ORDER_OPT, true,
				"sort order: coordinate, queryname, or unsorted"));
	}

	/**
	 * Gets the user-selected SortOrder.
	 */
	public SortOrder getSortOrder() {
		String argStr = getArgString(SORT_ORDER_OPT);
		if(argStr == null) return null;
		try {
			return SortOrder.valueOf(argStr);
		} catch(IllegalArgumentException e) {
			throw new IllegalArgumentException(
					"Invalid sort order option: " + argStr, e); 
		}
	}

}
