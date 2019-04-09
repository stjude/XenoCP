package org.stjude.compbio.sam.block;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.stjude.compbio.common.RefNameMap;
import org.stjude.compbio.common.formats.fasta.IndexedFastaFetcher;
import org.stjude.compbio.common.interval.RefInterval;

/**
 * RefBasesFetcher that uses an IndexedFastaFetcher
 */
public class IndexedFastaRefBasesFetcher implements RefBasesFetcher {
	private IndexedFastaFetcher fetcher;

	public IndexedFastaRefBasesFetcher(IndexedFastaFetcher fetcher) {
		this.fetcher = fetcher;
	}
	
	/**
	 * Opens a new IndexedFastaRefBasesFetcher on an indexed FASTA file.
	 * 
	 * @param fasta IndexedFastaSequenceFile object
	 * @param refNameMap RefNameMap used for canonicalization of ref names
	 */
	public IndexedFastaRefBasesFetcher(IndexedFastaSequenceFile fasta,
			RefNameMap<String> refNameMap) {
		this(new IndexedFastaFetcher(fasta, refNameMap));
	}
	
	/**
	 * Opens a new IndexedFastaRefBasesFetcher on an indexed FASTA file.
	 * 
	 * The RefNameMap is built from the corresponding sequence dictionary file,
	 * if found.
	 * 
	 * @param fasta IndexedFastaSequenceFile object
	 */
	public IndexedFastaRefBasesFetcher(IndexedFastaSequenceFile fasta) {
		this(new IndexedFastaFetcher(fasta));
	}

	public byte[] getRefBases(RefInterval interval) {
		return fetcher.fetchBytes(interval);
	}

}
