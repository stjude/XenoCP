package org.stjude.compbio.util;

import java.util.Arrays;

/**
 * Array helper functions not covered in java.util.Arrays
 * 
 * Concat function(s) courtesy of Joachim Sauer:
 * http://stackoverflow.com/questions/80476/how-to-concatenate-two-arrays-in-java
 */
public class ArrayHelper {
	/**
	 * Concatenates two arrays
	 * 
	 * @param first first array
	 * @param second second array
	 * @return concatenated array, or null if either array is null
	 */
	public static byte[] concat(byte[] first, byte[] second) {
		if(first == null || second == null) return null;
		byte[] result = Arrays.copyOf(first, first.length + second.length);
		System.arraycopy(second, 0, result, first.length, second.length);
		return result;
	}
}
