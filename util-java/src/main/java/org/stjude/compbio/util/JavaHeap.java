package org.stjude.compbio.util;

import java.io.PrintStream;

/**
 * Class with convenience methods for interacting with the Java heap
 */
public class JavaHeap {
	/**
	 * The Runtime object
	 */
	private static Runtime runtime = Runtime.getRuntime();
	
	/**
	 * Prints a standardized message to stderr
	 * @param label line label (may be null to omit)
	 */
	public static void report(String label) {
		report(System.err, label);
	}

	/**
	 * Prints a standardized message to a PrintStream
	 * 
	 * @param ps print stream to which to write
	 * @param label line label (may be null to omit)
	 */
	private static void report(PrintStream ps, String label) {
		long max = runtime.maxMemory() / 1024;
		long cur = runtime.totalMemory() / 1024;
		long curFree = runtime.freeMemory() / 1024;
		ps.printf("Used:%8d Free:%8d Maxi:%8d CrSz:%8d CrFr:%8d (K) %s\n",
				cur - curFree,
				curFree + max - cur,
				max,
				cur,
				curFree, 
				label == null ? "" : label);
	}
}
