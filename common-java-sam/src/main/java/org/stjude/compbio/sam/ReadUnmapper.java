package org.stjude.compbio.sam;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import org.stjude.compbio.common.seq.SequenceHelper;

import htsjdk.samtools.SAMRecord;

/**
 * Provides functionality for setting reads to unmapped.
 * 
 * In the future, we may want to use a better Set implementation, so that we
 * can scale up to larger name sets.
 */
public class ReadUnmapper {
	/**
	 * Read names to set unmapped
	 */
	private Set<String> names;
	
	public ReadUnmapper() {
		this.names = new HashSet<String>();
	}
	
	/**
	 * Gets the names of reads to be set to unmapped.
	 * @return live set
	 */
	public Set<String> getNames() {
		return names;
	}

	/**
	 * Sets the names of reads to be set to unmapped.
	 * @param names set of read names
	 */
	public void setNames(Set<String> names) {
		this.names = names;
	}
	
	/**
	 * Reads read names from a file into the set of names to set to unmapped.
	 * 
	 * If the set is currently non-empty then it is not replaced.
	 * 
	 * @param namesFile file containing the read names, one per line
	 * @throws IOException on any read problem
	 */
	public void readNamesFile(File namesFile) throws IOException {
		BufferedReader reader = new BufferedReader(
				new FileReader(namesFile));
		String name;
		while((name = reader.readLine()) != null) {
			names.add(name);
		}
		reader.close();
	}
	
	/**
	 * Returns true if the SAM record's read name is in the set
	 * @param record record to check
	 */
	public boolean hits(SAMRecord record) {
		return names.contains(record.getReadName());
	}
	
	/**
	 * Sets a read to unmapped (including mate information) and adds the XU tag.
	 * @param record record to update
	 */
	public static void setUnmapped(SAMRecord record) {
		// Set flags
		record.setReadUnmappedFlag(true);
		record.setMateUnmappedFlag(true);
		
		// Move mapping information to xu tag
		StringBuilder xu = new StringBuilder();
		// MAPQ
		xu.append("MAPQ=");
		xu.append(record.getMappingQuality());
		record.setMappingQuality(0);
		// CIGAR
		xu.append(";CIGAR=");
		xu.append(record.getCigarString());
		record.setCigarString(SAMRecord.NO_ALIGNMENT_CIGAR);
		// ISIZE
		xu.append(";ISIZE=");
		xu.append(record.getInferredInsertSize());
		record.setInferredInsertSize(0);
		// Add tag
		record.setAttribute("XU", xu.toString());
	}

	/**
	 * Creates a copy of the record that only has seq, qual, RG, and read
	 * pairing flags
	 * 
	 * @param record the record to strip
	 * @param applyReverseComplement if true and the record's reverse complement
	 * 		  flag is set, then the read bases will be reverse complemented and
	 *        qualities will be reversed.
	 */
	public static SAMRecord strip(SAMRecord record, boolean applyReverseComplement) {
		SAMRecord uRecord = new SAMRecord(record.getHeader());
		uRecord.setReadName(record.getReadName());
		uRecord.setFlags(record.getFlags() & 0x00c1);
		uRecord.setReadUnmappedFlag(true);
		uRecord.setMateUnmappedFlag(true);

		if(applyReverseComplement && record.getReadNegativeStrandFlag()) {
			uRecord.setReadBases(
					SequenceHelper.reverseComplement(record.getReadBases()));
			uRecord.setBaseQualityString(
					new StringBuilder(record.getBaseQualityString()).reverse().toString());
		}
		else {
			uRecord.setReadBases(record.getReadBases());
			uRecord.setBaseQualities(record.getBaseQualities());
		}

		String[] tagsToCopy = { "RG", "PG" };
		for(String tag: tagsToCopy) {
			Object attr = record.getAttribute(tag); 
			if(attr != null) uRecord.setAttribute(tag, attr);
		}
		
		return uRecord;
	}

}
