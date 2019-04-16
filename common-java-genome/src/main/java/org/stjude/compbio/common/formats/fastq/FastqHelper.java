package org.stjude.compbio.common.formats.fastq;

import htsjdk.samtools.fastq.FastqRecord;

/**
 * Helper class for FASTQ files
 */
public class FastqHelper {
	/**
	 * Gets the name of a FastqRecord
	 */
	public static String getName(FastqRecord record) {
		return record.getReadHeader().split("[ \t]")[0];
	}
}
