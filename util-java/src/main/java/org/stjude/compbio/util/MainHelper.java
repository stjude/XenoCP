package org.stjude.compbio.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Helper class for main methods
 */
public class MainHelper  {
	/**
	 * Command-line options, to be set by the implementation in the constructor
	 */
	private Options options;
	/**
	 * CommandLine parsed from args using options
	 */
	private CommandLine commandLine;
	/**
	 * Whether to throw an exception when an argument that does not exist is
	 * requested.
	 */
	private boolean strict;
	/**
	 * Whether to automatically make parent directories of output files when
	 * they don't exist.
	 */
	private boolean autoMakeOutputFileParents;
	
	/**
	 * Full constructor, using options, command line, and strictness
	 * 
	 * @param options command line options
	 * @param args actual arguments
	 * @param strict
	 */
	public MainHelper(Options options, String[] args, boolean strict) {
		this(options);
		setArgs(args);
		this.strict = strict;
		this.autoMakeOutputFileParents = false;
	}

	/**
	 * The normal constructor, using Options, defaults to strict
	 * 
	 * @param options command-line options this command accepts
	 */
	public MainHelper(Options options) {
		this.options = options;
		this.commandLine = null;
		this.strict = true;
	}

	/**
	 * Gets acceptable command-line options (not actual values)
	 * @return Options object
	 */
	public Options getOptions() {
		return options;
	}

	/**
	 * Sets acceptable command-line options
	 * @param options Options object
	 */
	public void setOptions(Options options) {
		this.options = options;
	}
	
	/**
	 * Whether arg getters/setters should throw an exception when the argument
	 * was not passed.
	 * @return value of current setting.
	 */
	public boolean isStrict() {
		return strict;
	}

	/**
	 * Sets whether arg getters/setters should throw an exception when the argument
	 * was not passed.
	 * @param strict new value (true for exception, false for return null)
	 */
	public void setStrict(boolean strict) {
		this.strict = strict;
	}

	/**
	 * Whether parent directories of output files should be created
	 * automatically
	 * @return value of current setting (default false) 
	 */
	public boolean isAutoMakeOutputFileParents() {
		return autoMakeOutputFileParents;
	}

	/**
	 * Sets whether parent directories of output files should be created
	 * automatically.
	 * @param autoMakeOutputFileParents new value (true for mkdir, flase for
	 * 			FileNotFoundException)
	 */
	public void setAutoMakeOutputFileParents(
			boolean autoMakeOutputFileParents) {
		this.autoMakeOutputFileParents = autoMakeOutputFileParents;
	}

	/**
	 * Gets CommandLine parsed from args using options
	 * @return CommandLine, or null if args not set
	 */
	public CommandLine getCommandLine() {
		return commandLine;
	}
	
	/**
	 * Gets the non-option argument list
	 * @return list of string non-option arguments
	 */
	@SuppressWarnings("unchecked")
	public List<String> getArgs() {
		return commandLine.getArgList();
	}

	/**
	 * Sets args, which makes CommandLine and other related methods available.
	 * 
	 * @param args specified args (from main)
	 * @throws CommandLineException if there is a parsing problem
	 */
	public void setArgs(String[] args) throws CommandLineException {
		CommandLineParser parser = new GnuParser();
		try {
			this.commandLine = parser.parse(options, args);
		} catch (ParseException e) {
			throw new CommandLineException(e);
		}
	}
	
	/**
	 * Returns true if an option argument is set
	 * @param optName option name
	 */
	public boolean hasArg(String optName) {
		return commandLine.hasOption(optName);
	}
	
	/**
	 * Gets a numbered non-option argument as a string
	 * 
	 * @param int pos 0-indexed position
	 * @return argument as a string, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 */
	public String getArgString(int pos) 
	throws CommandLineException, IllegalStateException {
		if(commandLine == null) throw new IllegalStateException("No args set");
		List<?> nonOptArgs = commandLine.getArgList();
		if(nonOptArgs.size() <= pos) {
			if(strict) throw new CommandLineException(
					"Not enough arguments to read argument # " + (pos + 1));
			return null;
		}
		return (String)nonOptArgs.get(pos);
	}
	
	/**
	 * Gets a named option argument as a string
	 * 
	 * @param optName option name
	 * @return argument as a string, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 */
	public String getArgString(String optName) 
	throws CommandLineException, IllegalStateException {
		if(commandLine == null) throw new IllegalStateException("No args set");
		if(strict && !commandLine.hasOption(optName)) {
			throw new CommandLineException("Option " + optName +
					" was not specified.");
		}
		return commandLine.getOptionValue(optName);
	}
	
	/**
	 * Gets a named option argument as a string, with a default if not strict
	 * 
	 * @param optName option name
	 * @param def default return value
	 * @return argument as a string, or def if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 */
	public String getArgString(String optName, String def) 
	throws CommandLineException, IllegalStateException {
		String arg = getArgString(optName);
		return (arg == null) ? def : arg;
	}
	
	/**
	 * Gets a numbered non-option argument as an integer
	 * 
	 * @param int pos 0-indexed position
	 * @return argument as an Integer, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 * @throws NumberFormatException if arg is not an integer
	 */
	public Integer getArgInt(int pos) 
	throws CommandLineException, IllegalStateException {
		String arg = getArgString(pos);
		if(arg == null) return null;
		return Integer.valueOf(arg);
	}
		
	/**
	 * Gets a named option argument as an integer
	 * @param optName option name
	 * @return argument as an Integer, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 * @throws NumberFormatException if arg is not an integer
	 */
	public Integer getArgInt(String optName) {
		String arg = getArgString(optName);
		if(arg == null) return null;
		return Integer.valueOf(arg);
	}
	
	/**
	 * Gets a named option argument as an int, with a default if not strict
	 * 
	 * @param optName option name
	 * @param def default return value
	 * @return argument as an int, or def if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 */
	public int getArgInt(String optName, int def) 
	throws CommandLineException, IllegalStateException {
		Integer arg = getArgInt(optName);
		return (arg == null) ? def : arg;
	}

	/**
	 * Gets a numbered non-option argument as a double
	 * 
	 * @param int pos 0-indexed position
	 * @return argument as a Double, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 * @throws NumberFormatException if arg is not a double
	 */
	public Double getArgDouble(int pos) 
	throws CommandLineException, IllegalStateException {
		String arg = getArgString(pos);
		if(arg == null) return null;
		return Double.valueOf(arg);
	}
		
	/**
	 * Gets a named option argument as an double
	 * @param optName option name
	 * @return argument as a Double, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 * @throws NumberFormatException if arg is not an double
	 */
	public Double getArgDouble(String optName) {
		String arg = getArgString(optName);
		if(arg == null) return null;
		return Double.valueOf(arg);
	}
	
	/**
	 * Gets a named option argument as an int, with a default if not strict
	 * 
	 * @param optName option name
	 * @param def default return value
	 * @return argument as a double, or def if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 */
	public double getArgDouble(String optName, double def) 
	throws CommandLineException, IllegalStateException {
		Double arg = getArgDouble(optName);
		return (arg == null) ? def : arg;
	}
	
	/**
	 * Gets a named option argument as an enum
	 * @param optName option name
	 * @param type enum type
	 * @return argument as a string, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 * @throws IllegalArgumentException if arg is not a valid enum value
	 */
	public <T extends Enum<T>> T getArgEnum(String optName, Class<T> enumType) {
		String arg = getArgString(optName);
		if(arg == null) return null;
		return Enum.valueOf(enumType, arg);
	}

	/**
	 * Gets a numbered non-option argument as a File
	 * 
	 * @param int pos 0-indexed position
	 * @return argument as a File, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 */
	public File getArgFile(int pos) 
	throws CommandLineException, IllegalStateException {
		String filename = getArgString(pos);
		if(filename == null) return null;
		return new File(filename);
	}
	
	/**
	 * Gets a named option argument as a File
	 * 
	 * @param optName option name
	 * @return argument as a File, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 */
	public File getArgFile(String optName) 
	throws CommandLineException, IllegalStateException {
		String filename = getArgString(optName);
		if(filename == null) return null;
		return new File(filename);
	}
	
	/**
	 * Validates an input file argument
	 * @param file file to validate
	 * @return processed argument
	 * @throws FileNotFoundException if file does not exist
	 */
	private File validateInputFile(File file)
			throws FileNotFoundException {
		if(file != null && !file.exists()) {
			throw new FileNotFoundException( 
					"Input file does not exist: " + file);
		}
		return file;
	}
	
	/**
	 * Gets and validates a non-option input file from the command line.
	 * 
	 * @param int pos 0-indexed position
	 * @return argument as a File, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 * @throws FileNotFoundException if the file does not exist
	 */
	public File getInputFile(int pos)
	throws CommandLineException, IllegalStateException, FileNotFoundException {
		return validateInputFile(getArgFile(pos));
	}
	
	/**
	 * Gets and validates an option input file from the command line.
	 * 
	 * @param optName option name
	 * @return argument as a File, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 * @throws FileNotFoundException if the file does not exist
	 */
	public File getInputFile(String optName)
	throws CommandLineException, IllegalStateException, FileNotFoundException {
		return validateInputFile(getArgFile(optName));
	}
	
	/**
	 * Performs output file arg validation and processing
	 * 
	 * @param file raw file arg (may be null)
	 * @return processed argument
	 * @throws FileNotFoundException if the file does not exist
	 */
	private File validateOutputFile(File file) throws FileNotFoundException {
		if(file != null) {
			File parent = file.getParentFile();
			if(parent != null && !parent.exists()) {
				if(autoMakeOutputFileParents) {
					if(parent.mkdirs()) return file;
				}
				throw new FileNotFoundException( 
						"Parent directory does not exist for output file: " +
								file);
			}
		}
		return file;
	}
	
	/**
	 * Gets and validates a non-option output file from the command line.
	 * 
	 * @param int pos 0-indexed position
	 * @return argument as a File, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 * @throws FileNotFoundException if the file does not exist
	 */
	public File getOutputFile(int pos)
	throws CommandLineException, IllegalStateException, FileNotFoundException {
		return validateOutputFile(getArgFile(pos));
	}
	
	/**
	 * Gets and validates an option output file from the command line.
	 * 
	 * @param optName option name
	 * @return argument as a File, or null if not strict and not set
	 * @throws CommandLineException if strict and that arg does not exist
	 * @throws IllegalStateException if args not set
	 * @throws FileNotFoundException if the file does not exist
	 */
	public File getOutputFile(String optName)
	throws CommandLineException, IllegalStateException, FileNotFoundException {
		return validateOutputFile(getArgFile(optName));
	}

	/**
	 * Print usage instructions to stderr, automatically showing all parameters
	 * using options.
	 * 
	 * @param pre text that should appear before the automatically-generated
	 *            options part
	 */
	public void usage(String pre) {
		usage(System.err, pre);
	}

	/**
	 * Print usage instructions, automatically showing all parameters
	 * using options.
	 * 
	 * @param pre text that should appear before the automatically-generated
	 *            options part
	 */
	public void usage(OutputStream stream, String pre) {
		usage(new PrintStream(stream), pre);
	}

	/**
	 * Helper method that prints usage instructions, automatically showing all
	 * parameters using options.
	 * 
	 * @param out the PrintStream to which to write help
	 * @param pre text that should appear before the automatically-generated
	 *            options part
	 */
	public void usage(PrintStream out, String pre) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.printHelp(pre, options);
	}

	/**
	 * Exception that wraps common command-line problems, allowing calling
	 * classes to handle only one exception, or to not handle it at all (it is
	 * a RuntimeException).
	 */
	public static final class CommandLineException extends RuntimeException {
		private static final long serialVersionUID = 1L;

		public CommandLineException() {
			super();
		}

		public CommandLineException(String message, Throwable cause) {
			super(message, cause);
		}

		public CommandLineException(String message) {
			super(message);
		}

		public CommandLineException(Throwable cause) {
			super(cause);
		}
		
	}

	/**
	 * Adds a long option to an Options object.
	 * 
	 * @param options Options object to add to
	 * @param longOpt long option name (e.g. my-option for --my-option)
	 * @param argLabel argument label, or null if it takes no arg
	 * @param description description
	 */
	@SuppressWarnings("static-access")
	public static void addLongOpt(Options options, String longOpt,
			String argLabel, String description) {
		OptionBuilder ob =
				OptionBuilder.withLongOpt(longOpt).withDescription(description);
		if(argLabel != null) {
			ob = ob.hasArg().withArgName(argLabel);
		}
		options.addOption(ob.create());
	}
}