package org.stjude.compbio.common.formats.seqdict;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import org.stjude.compbio.common.RefNameMap;


public class SeqDictRefNameMapFactory {
	/**
	 * Constructs a RefNameMap to map alternate forms to the form used in a
	 * SAMSequenceDictionary
	 */
	public static RefNameMap<String> newInstance(
			SAMSequenceDictionary seqDict) {
		RefNameMap<String> map = new RefNameMap<String>();
		for(SAMSequenceRecord seqRec: seqDict.getSequences()) {
			map.add(seqRec.getSequenceName());
		}
		return map;
	}
}
