package org.stjude.compbio.common.formats.bed;

/**
 * Represents a field in a BED file
 *
 * Names and descriptions from https://genome.ucsc.edu/FAQ/FAQformat.html#format1
 */
public enum BedField {

	// Required fields
	
	chrom, // The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
	chromStart, // The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
	chromEnd, // The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 

	// Optional BED fields are:

	name, // Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
	score, // A score between 0 and 1000.
	strand, // Defines the strand, - either '+' or '-'.
	thickStart, // The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
	thickEnd, // The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
	itemRgb, // An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
	blockCount, // The number of blocks (exons) in the BED line.
	blockSizes, // A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
	blockStarts, // A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
}
