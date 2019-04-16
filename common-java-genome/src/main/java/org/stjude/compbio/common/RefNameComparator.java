package org.stjude.compbio.common;

import java.util.Comparator;

/**
 * Compare implementation for reference names, which sorts numbers as
 * numbers, before X, before Y, before all other letters, alphabetically.
 */
public class RefNameComparator implements Comparator<String> {
	public int compare(String o1, String o2) {
		// Put in array to allow easily repeating manipulations of both
		// values
		String[] chrs = { o1, o2 };
		int[] chrNums = new int[2];
		
		// Pre-process
		for(int i = 0; i < chrs.length; i++) {
			// Strip off leading "chr"
			if(chrs[i].startsWith("chr")) chrs[i] = chrs[i].substring(3);
		
			// Compare as numbers if possible
			try {
				chrNums[i] = Integer.parseInt(chrs[i]);
			} catch(NumberFormatException e) {
				chrNums[i] = Integer.MAX_VALUE;
			}
		}
		
		// Compare as numbers if possible
		int cmp = chrNums[0] - chrNums[1];
		if(cmp != 0) return cmp;
		
		// Otherwise, sort X before Y before other
		for(int i = 0; i < chrs.length; i++) {
			if(chrs[i] == "X" || chrs[i] == "x") {
				chrNums[i] = 'X';
			}
			else if(chrs[i] == "Y" || chrs[i] == "y") {
				chrNums[i] = 'Y';
			}
			else {
				chrNums[i] = Integer.MAX_VALUE;
			}
		}
		cmp = chrNums[0] - chrNums[1];
		if(cmp != 0) return cmp;
		
		// Finally, just do alphabetical
		return chrs[0].compareTo(chrs[1]);
	}
	
}