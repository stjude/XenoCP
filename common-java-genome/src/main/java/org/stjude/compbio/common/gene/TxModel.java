package org.stjude.compbio.common.gene;

import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.RefStranded;
import org.stjude.compbio.common.RefStrand;
import org.stjude.compbio.common.interval.Interval;
import org.stjude.compbio.common.interval.LongInterval;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.SubRefInterval;
import org.stjude.compbio.util.DelimitedListBuilder;

/**
 * Represents a single transcript model.
 * 
 * Transcripts are defined by their reference/strand/intervals, which are
 * unmodifiable.  Other transcript properties are modifiable.
 */
public class TxModel implements RefStranded {
	/**
	 * Reference name/strand
	 */
	private final RefStranded refStrand;
	/**
	 * Exon intervals
	 */
	private final LinkedList<SubRefInterval> exons;
	/**
	 * CDS interval
	 */
	private final RefInterval cdsInterval;
	
	// Lazily-initialized fields
	private RefInterval txInterval = null;
	private LinkedList<SubRefInterval> cdsIntervals = null;
	private LinkedList<SubRefInterval> cdsIntervals5to3 = null;

	/**
	 * Constructor using all fields
	 */
	public TxModel(String refName, Strand strand,
			List<LongInterval> exons,
			Interval<? extends Long> cdsInterval) {
		// Assign name/strand
		this.refStrand = new RefStrand(refName, strand);
		
		// Build up SubRefInterval list of exons by copying/converting from
		// input exon interval list
		this.exons = new LinkedList<SubRefInterval>();
		int ordinal = 0;
		for(Interval<? extends Long> exon: exons) {
			this.exons.add(new SubRefInterval(this, ordinal++, exon));
		}
		
		// Normalize 0-length CDS interval to null
		if(cdsInterval == null) {
			this.cdsInterval = null;
		}
		else if(cdsInterval.getLeft().compareTo(cdsInterval.getRight()) >= 0) {
			this.cdsInterval = null;
		}
		else {
			this.cdsInterval = new SubRefInterval(this, null, cdsInterval);
		}
	}

	public String getRefName() {
		return refStrand.getRefName();
	}

	public Strand getStrand() {
		return refStrand.getStrand();
	}

	/**
	 * Returns the exons as an unmodifiable list
	 */
	public List<SubRefInterval> getExons() {
		return Collections.unmodifiableList(exons);
	}

	/**
	 * Returns a copy of the exon intervals in an array
	 */
	public RefInterval[] getExonsArray() {
		return exons.toArray(new RefInterval[0]);
	}
	
	/**
	 * Returns a numbered exon
	 * 
	 * @param index exon index, where 0 is the leftmost
	 */
	public RefInterval getExon(int index) {
		return exons.get(index);
	}

	/**
	 * Returns an interval covering the total span of the CDS.
	 */
	public RefInterval getCdsInterval() {
		return cdsInterval;
	}

	/**
	 * Returns an interval covering the total span of the transcript
	 */
	public RefInterval getTxInterval() {
		// Lazy initialization
		if(txInterval == null) {
			txInterval = new SubRefInterval(this,
					exons.getFirst().getLeft(), exons.getLast().getRight());
		}
		return txInterval;
	}
	
	/**
	 * Determines and returns the index of the exon that contains the CDS start
	 * 
	 * @return index of exon containing cds start, or null if no cds
	 */
	public Integer getCdsStartExon() {
		getCdsIntervals5to3();
		if(cdsIntervals5to3 == null) return null;
		return cdsIntervals5to3.getFirst().getOrdinal();
	}
	
	/**
	 * Determines and returns the index of the exon that contains the CDS end
	 * 
	 * @return index of exon containing cds end, or null if no cds
	 */
	public Integer getCdsEndExon() {
		getCdsIntervals5to3();
		if(cdsIntervals5to3 == null) return null;
		return cdsIntervals5to3.getLast().getOrdinal();
	}
	
	/**
	 * Returns a list of CDS intervals from left to right
	 */
	public List<SubRefInterval> getCdsIntervals() {
		// If no cds, leave at null
		if(cdsInterval == null) return null;
		
		// If not yet determined, lazily initialize
		if(cdsIntervals == null) {
			cdsIntervals = new LinkedList<SubRefInterval>();
			for(SubRefInterval exon: exons) {
				if(exon.overlaps(cdsInterval)) {
					cdsIntervals.add(exon.intersect(cdsInterval));
				}
			}
		}
		
		return Collections.unmodifiableList(cdsIntervals);
	}
	
	/**
	 * Returns a list of CDS intervals from left to right
	 */
	public List<SubRefInterval> getCdsIntervals5to3() {
		// Lazy initialization
		if(cdsIntervals5to3 == null) {
			// Return null if no CDS
			List<SubRefInterval> lToR = getCdsIntervals();
			if(lToR == null) return null;

			// Reverse if necessary
			if(refStrand.getStrand().isReverse()) {
				cdsIntervals5to3 = new LinkedList<SubRefInterval>(lToR);
				Collections.reverse(cdsIntervals5to3);
			}
			else {
				cdsIntervals5to3 = cdsIntervals;
			}
		}
		
		return cdsIntervals5to3;
	}

	/**
	 * For a given 0-based left-to-right ordinal, get the corresponding 1-based
	 * start-to-stop ordinal
	 * 
	 * @param ordinal 0-based ordinal as used internally
	 * @return 1-based 5'-to-3' ordinal for human use
	 */
	public int to5to3Ordinal(int ordinal) {
		if(refStrand.getStrand().isReverse()) {
			return exons.size() - ordinal;
		}
		else {
			return ordinal + 1;
		}
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((refStrand == null) ? 0 : refStrand.hashCode());
		result = prime * result + ((exons == null) ? 0 : exons.hashCode());
		result = prime * result
				+ ((cdsInterval == null) ? 0 : cdsInterval.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if(this == obj) return true;
		if(obj == null) return false;
		if(getClass() != obj.getClass()) return false;
		
		TxModel other = (TxModel)obj;
		
		if(refStrand == null) {
			if(other.refStrand != null) return false;
		}
		else if(!refStrand.equals(other.refStrand)) {
			return false;
		}
		
		if(exons == null) {
			if(other.exons != null) return false;
		}
		else if(!exons.equals(other.exons)) {
			return false;
		}
		
		if(cdsInterval == null) {
			if(other.cdsInterval != null) return false;
		}
		else if(!cdsInterval.equals(other.cdsInterval)) {
			return false;
		}
		
		return true;
	}

	@Override
	public String toString() {
		DelimitedListBuilder dlb = new DelimitedListBuilder(
				",", refStrand.toString() + ":", false);
		for(RefInterval exon: exons) {
			dlb.append(exon.getLeft() + "-" + exon.getRight());
		}
		return dlb.toString();
	}
}
