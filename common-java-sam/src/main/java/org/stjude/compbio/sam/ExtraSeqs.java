package org.stjude.compbio.sam;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.stjude.compbio.util.DelimitedListBuilder;

import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.SAMRecord;

/**
 * Support for the extra seqs SAMRecord attribute
 */
public class ExtraSeqs {
	public static final String ATTR_NAME = MoreSAMTag.EXTRA_SEQS.name();
	public static final String ENTRY_DELIM = "~";
	public static final String PART_DELIM = "|";
	
	private Map<String, FastqRecord> seqs;
	
	public ExtraSeqs() {
		seqs = new HashMap<String, FastqRecord>();
	}
	
	public Map<String, FastqRecord> getSeqs() {
		return Collections.unmodifiableMap(seqs);
	}

	public void add(String label, String seq, String qual) {
		add(label, new FastqRecord("", seq, "", qual));
	}

	public void add(String label, FastqRecord record) {
		seqs.put(label, record);
	}
	
	public boolean addIfNotEmpty(String label, FastqRecord record) {
		if(record.getReadString().length() == 0) return false;
		add(label, record);
		return true;
	}
	
	/**
	 * Returns the extra seqs as a single string value
	 */
	@Override
	public String toString() {
		DelimitedListBuilder dlb = new DelimitedListBuilder(ENTRY_DELIM);
		for(String label: seqs.keySet()) {
			String seq = seqs.get(label).getReadString();
			String qual = seqs.get(label).getBaseQualityString();
			dlb.append(label + PART_DELIM + seq + PART_DELIM + qual);
		}
		return dlb.toString();
	}
	
	/**
	 * Parses a string representation
	 * @param string
	 * @return
	 */
	public static ExtraSeqs parse(String string) {
		ExtraSeqs out = new ExtraSeqs();
		String[] entries = string.split("[" + ENTRY_DELIM + "]");
		for(String entry: entries) {
			if(entry.equals("")) continue;
			String[] parts = entry.split("[" + PART_DELIM + "]");
			out.add(parts[0],
					(parts.length > 1 ? parts[1] : ""),
					(parts.length > 2 ? parts[2] : ""));
		}
		return out;
	}
	
	/**
	 * Sets attribute on a SAMRecord
	 */
	public void setSAMRecordAttribute(SAMRecord record) {
	    record.setAttribute(ATTR_NAME, toString());
	}
	
	/**
	 * Parses from a SAMRecord
	 */
	public static ExtraSeqs parse(SAMRecord record) {
		Object attrib = record.getAttribute(ATTR_NAME);
		return parse(attrib == null ? "" : attrib.toString());
	}
	
	public static Map<String, FastqRecord> decode(SAMRecord record) {
		return parse(record).getSeqs();
	}

	/**
	 * Decodes an XX attribute
	 * 
	 * @param encoded XX attribute
	 * @return Map from label to FastqRecord holding seq and qual
	 */
	public static Map<String, FastqRecord> decode(String encoded) {
		return parse(encoded).getSeqs();
	}
}
