package org.stjude.compbio.common.formats.psl;

/**
 * Represents a field in a psl file
 *
 * Names and descriptions from http://useast.ensembl.org/info/website/upload/psl.html
 * 
 * TODO rename and/or change descriptions
 */
public enum PslField {
	matches, // Number of matching bases that aren't repeats.
	misMatches, // Number of bases that don't match.
	repMatches, // Number of matching bases that are part of repeats.
	nCount, // Number of 'N' bases.
	qNumInsert, // Number of inserts in query.
	qBaseInsert, // Number of bases inserted into query.
	tNumInsert, // Number of inserts in target.
	tBaseInsert, // Number of bases inserted into target.
	strand, // defined as + (forward) or - (reverse) for query strand. In mouse,
			// a second '+' or '-' indicates genomic strand.
	qName, // Query sequence name.
	qSize, // Query sequence size.
	qStart, // Alignment start position in query.
	qEnd, // Alignment end position in query.
	tName, // Target sequence name.
	tSize, // Target sequence size.
	tStart, // Alignment start position in query.
	tEnd, // Alignment end position in query.
	blockCount, // Number of blocks in the alignment.
	blockSizes, // Comma-separated list of sizes of each block.
	qStarts, // Comma-separated list of start position of each block in query.
	tStarts // Comma-separated list of start position of each block in target.
}
