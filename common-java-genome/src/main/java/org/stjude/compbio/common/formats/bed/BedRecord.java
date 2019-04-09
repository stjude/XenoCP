package org.stjude.compbio.common.formats.bed;

import java.util.SortedMap;

import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.interval.Interval;
import org.stjude.compbio.common.interval.LongInterval;
import org.stjude.compbio.common.interval.NamedRefInterval;
import org.stjude.compbio.common.interval.StdNamedRefInterval;
import org.stjude.compbio.util.io.IoField;
import org.stjude.compbio.util.io.IoFormat;
import org.stjude.compbio.util.io.IoRecord;

public class BedRecord extends IoRecord<BedField> implements NamedRefInterval {

	private NamedRefInterval interval;
	
	public BedRecord(IoFormat<BedField> format,
			SortedMap<IoField<BedField>, String> data) {
		super(format, data);
		initRefInterval();
	}

	/**
	 * Initializes RefInterval
	 */
	private void initRefInterval() {
		// Extract information
		String refName = valueForKey(BedField.chrom);
		long left = Long.parseLong(valueForKey(BedField.chromStart));
		long right = Long.parseLong(valueForKey(BedField.chromEnd));
		String strandPlusMinus = valueForKey(BedField.strand);
		
		// Parse strand
		Strand strand;
		try {
			strand = Strand.valueOfStr(strandPlusMinus);
		}
		catch(IllegalArgumentException e) {
			strand = Strand.NA;
		}
		
		// Construct the RefInterval
		this.interval = new StdNamedRefInterval(refName, strand, left, right, getName());
	}
	
	/**
	 * Gets the name
	 */
	public String getName() {
		return valueForKey(BedField.name);
	}
	
	

	
	/* Delegate methods to "interval", used to implement RefInterval */
	
	/**
	 * @param o
	 * @return
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	public int compareTo(Interval<Long> o) {
		return interval.compareTo(o);
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.LongInterval#getLeftInt()
	 */
	public int getLeftInt() {
		return interval.getLeftInt();
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.RefStranded#getRefName()
	 */
	public String getRefName() {
		return interval.getRefName();
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.LongInterval#getRightInt()
	 */
	public int getRightInt() {
		return interval.getRightInt();
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.RefStranded#getStrand()
	 */
	public Strand getStrand() {
		return interval.getStrand();
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.LongInterval#length()
	 */
	public long length() {
		return interval.length();
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.LongInterval#getMidpoint()
	 */
	public long getMidpoint() {
		return interval.getMidpoint();
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#getLeftInBase()
	 */
	public long getLeftInBase() {
		return interval.getLeftInBase();
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.Interval#getLeft()
	 */
	public Long getLeft() {
		return interval.getLeft();
	}

	/**
	 * @param other
	 * @return
	 * @see org.stjude.compbio.common.interval.LongInterval#overlapLength(org.stjude.compbio.common.interval.Interval)
	 */
	public long overlapLength(Interval<? extends Long> other) {
		return interval.overlapLength(other);
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#getLeftInBaseInt()
	 */
	public int getLeftInBaseInt() {
		return interval.getLeftInBaseInt();
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.Interval#getRight()
	 */
	public Long getRight() {
		return interval.getRight();
	}

	/**
	 * @param point
	 * @return
	 * @see org.stjude.compbio.common.interval.Interval#isEndpoint(java.lang.Comparable)
	 */
	public boolean isEndpoint(Long point) {
		return interval.isEndpoint(point);
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#getStart()
	 */
	public long getStart() {
		return interval.getStart();
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#getStartInBase()
	 */
	public long getStartInBase() {
		return interval.getStartInBase();
	}

	/**
	 * @param point
	 * @return
	 * @see org.stjude.compbio.common.interval.Interval#containsClosed(java.lang.Comparable)
	 */
	public boolean containsClosed(Long point) {
		return interval.containsClosed(point);
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#getEnd()
	 */
	public long getEnd() {
		return interval.getEnd();
	}

	/**
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#getEndInBase()
	 */
	public long getEndInBase() {
		return interval.getEndInBase();
	}

	/**
	 * @param point
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#containsInBase(long)
	 */
	public boolean containsInBase(long point) {
		return interval.containsInBase(point);
	}

	/**
	 * @param point
	 * @return
	 * @see org.stjude.compbio.common.interval.Interval#containsOpen(java.lang.Comparable)
	 */
	public boolean containsOpen(Long point) {
		return interval.containsOpen(point);
	}

	/**
	 * @param point
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#getStartInterval(java.lang.Long)
	 */
	public NamedRefInterval getStartInterval(Long point) {
		return interval.getStartInterval(point);
	}

	/**
	 * @param other
	 * @return
	 * @see org.stjude.compbio.common.interval.Interval#contains(org.stjude.compbio.common.interval.Interval)
	 */
	public boolean contains(Interval<? extends Long> other) {
		return interval.contains(other);
	}

	/**
	 * @param point
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#getEndInterval(java.lang.Long)
	 */
	public NamedRefInterval getEndInterval(Long point) {
		return interval.getEndInterval(point);
	}

	/**
	 * @param amount
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#shift(long)
	 */
	public NamedRefInterval shift(long amount) {
		return interval.shift(amount);
	}

	/**
	 * @param amount
	 * @return
	 * @see org.stjude.compbio.common.interval.LongInterval#expandLeft(long)
	 */
	public LongInterval expandLeft(long amount) {
		return interval.expandLeft(amount);
	}

	/**
	 * @param amount
	 * @return
	 * @see org.stjude.compbio.common.interval.LongInterval#expandRight(long)
	 */
	public LongInterval expandRight(long amount) {
		return interval.expandRight(amount);
	}

	/**
	 * @param amount
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#expandStart(long)
	 */
	public NamedRefInterval expandStart(long amount) {
		return interval.expandStart(amount);
	}

	/**
	 * @param amount
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#expandEnd(long)
	 */
	public NamedRefInterval expandEnd(long amount) {
		return interval.expandEnd(amount);
	}

	/**
	 * @param other
	 * @return
	 * @see org.stjude.compbio.common.interval.Interval#overlaps(org.stjude.compbio.common.interval.Interval)
	 */
	public boolean overlaps(Interval<? extends Long> other) {
		return interval.overlaps(other);
	}

	/**
	 * @param other
	 * @return
	 * @throws IllegalArgumentException
	 * @see org.stjude.compbio.common.interval.RefInterval#intersect(org.stjude.compbio.common.interval.Interval)
	 */
	public NamedRefInterval intersect(Interval<? extends Long> other) throws IllegalArgumentException {
		return interval.intersect(other);
	}

	/**
	 * @param other
	 * @return
	 * @see org.stjude.compbio.common.interval.Interval#overlapsOrAbuts(org.stjude.compbio.common.interval.Interval)
	 */
	public boolean overlapsOrAbuts(Interval<? extends Long> other) {
		return interval.overlapsOrAbuts(other);
	}

	/**
	 * @param other
	 * @return
	 * @see org.stjude.compbio.common.interval.RefInterval#union(org.stjude.compbio.common.interval.Interval)
	 */
	public NamedRefInterval union(Interval<? extends Long> other) {
		return interval.union(other);
	}	
}
