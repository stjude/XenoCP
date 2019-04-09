package org.stjude.compbio.common.formats.fasta;

import java.io.File;
import java.io.FileNotFoundException;

import org.stjude.compbio.common.RefNameMap;
import org.stjude.compbio.common.formats.seqdict.SeqDictRefNameMapFactory;
import org.stjude.compbio.common.interval.RefInterval;

import picard.PicardException;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.SequenceUtil;

/**
 * Wrapper around an IndexedFastaSequenceFile that uses RefIntervals and handles
 * chromosome naming convention inconsistencies.
 * 
 */
public class IndexedFastaFetcher {
	/**
	 * The underlying IndexedFastaSequenceFile object
	 */
	private final IndexedFastaSequenceFile fasta;
	/**
	 * The RefNameMap used for canonicalization of ref names
	 */
	private final RefNameMap<String> refNameMap;
	
	/**
	 * Gets a RefNameMap from the sequence dictionary, if possible
	 * 
	 * @param fasta IndexedFastaSequenceFile object
	 * @return RefNameMap, which will be filled if there is a seq dict
	 */
	private static RefNameMap<String> getRefNameMapFromSeqDict(
			IndexedFastaSequenceFile fasta) {
		SAMSequenceDictionary seqDict = fasta.getSequenceDictionary();
		if(seqDict == null) return new RefNameMap<String>();
		return SeqDictRefNameMapFactory.newInstance(seqDict);
	}

	/**
	 * @param fasta IndexedFastaSequenceFile object
	 * @param refNameMap RefNameMap used for canonicalization of ref names
	 */
	public IndexedFastaFetcher(IndexedFastaSequenceFile fasta,
			RefNameMap<String> refNameMap) {
		this.fasta = fasta;
		this.refNameMap = refNameMap;
	}
	
	/**
	 * Opens a new IndexedFastaFetcher on an indexed FASTA file.
	 * 
	 * The RefNameMap is built from the corresponding sequence dictionary file,
	 * if found.
	 * 
	 * @param fastaFile FASTA file, must be indexed
	 */
	public IndexedFastaFetcher(IndexedFastaSequenceFile fasta) {
		this(fasta, getRefNameMapFromSeqDict(fasta));
		
	}

	/**
	 * Opens a new IndexedFastaFetcher on an indexed FASTA file.
	 * 
	 * The RefNameMap is built from the corresponding sequence dictionary file,
	 * if found.
	 * 
	 * @param fastaFile FASTA file, must be indexed
	 * @throws FileNotFoundException if the FASTA or index file is not found
	 */
	public IndexedFastaFetcher(File fastaFile) throws FileNotFoundException {
		this(new IndexedFastaSequenceFile(fastaFile));
	}
	
	/**
	 * Constructs an IndexedFastaFetcher from the given seqDb, which is assumed
	 * to be an indexed FASTA file.
	 * 
	 * @param seqDb indexed FASTA file
	 * @return new IndexedFastaFetcher
	 * @throws IllegalArgumentException if the seqDb is not an indexed FASTA
	 */
	public static IndexedFastaFetcher newInstance(ReferenceSequenceFile seqDb)
	throws IllegalArgumentException {
		try {
			return new IndexedFastaFetcher((IndexedFastaSequenceFile)seqDb);
		} catch(ClassCastException e) {
			throw new IllegalArgumentException(
					"Reference sequence file is not an indexed FASTA file", e);
		}
	}

	/**
	 * Fetches sequence using a RefInterval
	 * 
	 * If the interval is on the reverse strand, then the sequence if reverse
	 * complemented.
	 * 
	 * @param refInterval the reference interval
	 * @return sequence as a string
	 */
	public String fetchString(RefInterval refInterval) {
		return fetchString(
				refInterval.getRefName(),
				refInterval.getLeftInBase(),
				refInterval.getRight(),
				refInterval.getStrand().isReverse());
	}

	/**
	 * Fetches sequence using ref name, coordinates, and reverse complement
	 * option.
	 * 
	 * The reference name is canonicalized, if possible, to make it match the
	 * FASTA names.
	 * 
	 * @param refName reference name
	 * @param left in-base left end
	 * @param right in-base right end
	 * @param rc whether to reverse complement
	 */
	public String fetchString(String refName, long left, long right,
			boolean rc) {
		// Canonicalize ref name using map if possible
		if(refNameMap.containsKey(refName)) {
			refName = refNameMap.get(refName);
		}

		// Pull out the subsequence, in reference orientation
		ReferenceSequence subseqObj;
		try {
			// Try to get the sequence using the indicated refName
			subseqObj = fasta.getSubsequenceAt(refName, left, right);
		} catch(PicardException e) {
			// If that fails, then try the alternate version (w/ or w/out "chr")
			String altName;
			try {
				altName = RefNameMap.altForm(refName);
			} catch(PicardException f) {
				throw new PicardException(
						"Second sequence fetch attempt failed for " + refName +
						":" + left + "-" + right + "; first attempt message: " +
						e.getMessage(), f);
			}
			subseqObj = fasta.getSubsequenceAt(altName, left, right);
			// If that fails, then the exception will be thrown, and we won't
			// get this far...
			
			// Since it worked with the alternate version, remember in the
			// RefNameMap
			refNameMap.add(altName);
		}
		
		// Convert sequence to string
		String subseq = new String(subseqObj.getBases());
		
		// Reverse comp if specified
		if(rc) subseq = SequenceUtil.reverseComplement(subseq);
		
		return subseq;
	}

	/**
	 * Fetches sequence using a RefInterval
	 * 
	 * If the interval is on the reverse strand, then the sequence if reverse
	 * complemented.
	 * 
	 * @param refInterval the reference interval
	 * @return sequence as a byte array
	 */
	public byte[] fetchBytes(RefInterval refInterval) {
		return fetchString(refInterval).getBytes();
	}

	/**
	 * Fetches sequence using ref name, coordinates, and reverse complement
	 * option.
	 * 
	 * The reference name is canonicalized, if possible, to make it match the
	 * FASTA names.
	 * 
	 * @param refName reference name
	 * @param left in-base left end
	 * @param right in-base right end
	 * @param rc whether to reverse complement
	 */
	public byte[] fetchBytes(String refName, long left, long right,
			boolean rc) {
		return fetchString(refName, left, right, rc).getBytes();
	}
}
