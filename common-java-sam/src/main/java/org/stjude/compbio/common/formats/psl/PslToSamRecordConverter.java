package org.stjude.compbio.common.formats.psl;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;

import org.stjude.compbio.sam.block.Blocklist;
import org.stjude.compbio.sam.block.BlocklistFactory;
import org.stjude.compbio.util.io.IoRecord;

/**
 * Class with helper functions for converting PSL to SAM records.
 *
 */
public class PslToSamRecordConverter {
	/**
	 * BlocklistFactory to use
	 */
	private static BlocklistFactory factory = new BlocklistFactory(false, null);

	public static Cigar createCigar(IoRecord<PslField> pslRecord) {
		// Parse relevant fields
		Iterator<Integer> blockSizes = parseIntList(
				pslRecord.valueForKey(PslField.blockSizes)).iterator();
		Iterator<Integer> readLefts = parseIntList(
				pslRecord.valueForKey(PslField.qStarts)).iterator();
		Iterator<Integer> refLefts = parseIntList(
				pslRecord.valueForKey(PslField.tStarts)).iterator();
		
		// Create the empty blocklist
		int refLeft = refLefts.next();
		Blocklist blocklist = factory.newEmptyBlocklist(refLeft, 0);
		
		// Add left soft-clip if necessary
		int readLeft = readLefts.next();
		if(readLeft > 0) {
			blocklist.add(readLeft, CigarOperator.S);
		}
		
		// Add first match
		blocklist.add(blockSizes.next(), CigarOperator.M);
		
		// Add (I/D, M) pairs until the end
		while(refLefts.hasNext()) {
			// Add I or D
			int nextRef = refLefts.next();
			int nextRead = readLefts.next();
			if(nextRef > refLeft) {
				blocklist.add(nextRef - refLeft, CigarOperator.D);
				refLeft = nextRef;
			}
			else {
				blocklist.add(nextRead - readLeft, CigarOperator.I);
				readLeft = nextRead;
			}
			// Add M
			blocklist.add(blockSizes.next(), CigarOperator.M);
		}
		
		// Convert to Cigar and return
		return blocklist.toCigar();
	}

	/**
	 * Converts a string containing a comma-separated list to a list of Integers
	 * 
	 * @param listStr list in string form
	 * @return list of Integers
	 */
	public static List<Integer> parseIntList(String listStr) {
		String[] strs = listStr.split(",");
		List<Integer> out = new ArrayList<Integer>(strs.length);
		for(String str: strs) {
			out.add(Integer.parseInt(str));
		}
		return out;
	}
}
