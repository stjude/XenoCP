package org.stjude.compbio.common;

/**
 * Interface for objects that have a reference name and Strand
 */
public interface RefStranded {
	/**
	 * Gets the reference name
	 */
	String getRefName();
	
	/**
	 * Gets the strand/orientation relative to the reference.
	 * 
	 * By convention, this should not be null, but it may be RefStranded.NA
	 */
	Strand getStrand();
}
