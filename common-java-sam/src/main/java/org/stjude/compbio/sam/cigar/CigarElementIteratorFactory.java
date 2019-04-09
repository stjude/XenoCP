package org.stjude.compbio.sam.cigar;

import java.util.Iterator;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;

/**
 * Class for constructing iterators over Cigar with elements returned at varying
 * levels of granularity.
 * 
 * Only COARSE and SOURCE granularity are supported.
 * 
 * @see CigarElementGranularity
 */
public class CigarElementIteratorFactory {
	/**
	 * Level of granularity
	 */
	private CigarElementGranularity granularity;

	private CigarElementIteratorFactory(CigarElementGranularity granularity) {
		this.granularity = granularity;
	}

	/**
	 * Constructs a factory with SOURCE granularity
	 */
	public CigarElementIteratorFactory() {
		this(CigarElementGranularity.SOURCE);
	}

	public CigarElementGranularity getGranularity() {
		return granularity;
	}

	public void setGranularity(CigarElementGranularity granularity) {
		this.granularity = granularity;
	}
	
	/**
	 * Constructs a new iterator over a Cigar with the proper level of 
	 * granularity
	 * 
	 * @param cigar input Cigar
	 * @return iterator over its elements, possibly split or merged
	 */
	public Iterator<CigarElement> newIterator(Cigar cigar) {
		// Get source iterator
		Iterator<CigarElement> iterator = cigar.getCigarElements().iterator();
		
		switch(granularity) {
		case COARSE:
			return new CoarseCigarElementIterator(iterator);
			
		case FINE:
			throw new IllegalArgumentException(
					"Cannot create fine granularity iterator; need more " +
					"information (consider converting to Blocklist first).");
			
		case SOURCE:
		default:
			return iterator;
		}
	}
}
