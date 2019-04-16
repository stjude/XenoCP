package org.stjude.compbio.sam.block;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

import org.stjude.compbio.common.formats.fasta.IndexedFastaFetcher;
import org.stjude.compbio.common.interval.StdLongInterval;
import org.stjude.compbio.sam.MdElt;
import org.stjude.compbio.sam.MdElts;

import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

public class BlockTestHelper {
	private static final File FASTA_FILE = new File(
			BlockTestHelper.class.getResource(
					"hg19_chr1_1_100120.fa").getPath());
	
	/**
	 * The test record
	 */
	public SAMRecord record;
	
	public BlockTestHelper() {
		resetRecord();
	}
	
	/**
	 * Resets the SAMRecord used for testing to a new record.
	 * 
	 * It is set to aligned at 1:100001 in-base, which means interbase left is
	 * 1000000.  You need to set the Cigar, MD, and readBases.
	 * 
	 * Sequence around this position is as follows:
	 * 
	 *  (interbase)
	 *  99960-100000 AGAAAGGATG AATCTTTCTG AAGGTTATGT CATCACACTC
	 * 100000-100040 ACTAAGCACA CAGAGAATAA TGTCTAGAAT CTGAGTGCCA
	 * 100040-100080 TGTTATCAAA TTGTACTGAG ACTCTTGCAG TCACACAGGC
	 * 100080-100120 TGACATGTAA GCATCGCCAT GCCTAGTACA GACTCTCCCT
	 */
	public void resetRecord() {
		record = new SAMRecord(null);
		record.setReferenceName("1");
		record.setAlignmentStart(100001);
	}
	
	/**
	 * Sets up a record that looks like this:
	 * 
	 * CIGAR: 10H    10M    2I   8M   2D   10M    50N   10M      10S
	 * Read:      ACTAAGCACATGCAGCGAAT**TGTCCTGAAT   TGTCAAGTAAACGTACGTAC
	 * Mismatch:                 X          XX         X  X
	 * Ref:       ACTAAGCACA**CAGAGAATAATGTCTAGAAT   TGACATGTAA
	 * MD:               13      A 4 ^AA  4 TA   6     A2 T 4
	 */
	public void setStandardRecordInfo() {
		setRecordInfo(
				"10H10M2I8M2D10M50N10M10S",
				"ACTAAGCACATGCAGCGAATTGTCCTGAATTGTCAAGTAAACGTACGTAC",
				"13A4^AA4TA6A2T4");
	}

	/**
	 * Sets CIGAR, read bases, and MD for the record
	 * 
	 * @param cigarStr CIGAR as a string
	 * @param readBasesStr read bases as a string
	 * @param mdStr MD as a string, or null for no MD
	 */
	public void setRecordInfo(
			String cigarStr, String readBasesStr, String mdStr) {
		record.setCigarString(cigarStr);
		record.setReadBases(readBasesStr.getBytes());
		if(mdStr != null) record.setAttribute(SAMTag.MD.toString(), mdStr);
	}

	/**
	 * Gets the MD string from the record
	 */
	public String getMdString() {
		return (String)record.getAttribute("MD");
	}
	
	/**
	 * Gets MdElts parsed from record
	 */
	public List<MdElt> getMdElts() {
		return MdElts.parse(getMdString());
	}
	
	/**
	 * Gets MdElts parsed from record, and after consuming some number of bases
	 * 
	 * @param from starting position in terms of ref bases
	 */
	public List<MdElt> getMdElts(int from) {
		return MdElts.sublist(getMdElts(), from);
	}
	
	/**
	 * Gets a CIGAR element (convenience method)
	 */
	public CigarElement getCigarElement(int i) {
		return record.getCigar().getCigarElement(i);
	}

	/**
	 * Asserts CIGAR info about a block
	 */
	public void assertBlockCigar(String expCigarEltStr, Block block) {
		String cigarEltStr =
				"" + block.getCigarElement().getLength() + block.getOperator();
		assertEquals(expCigarEltStr, cigarEltStr);
	}
	
	/**
	 * Asserts all info about a block
	 */
	public void assertBlock(Block block, CigarElement ce,
			long readL, long readR, long refL, long refR,
			String readBases, String refBases, BlockForm covers,
			int len, Integer matches) {
		assertEquals(ce, block.getCigarElement());
		assertEquals(new StdLongInterval(readL, readR), block.getReadInterval());
		assertEquals(new StdLongInterval(refL, refR), block.getRefInterval());
		assertNotNull(block.getReadBases());
		assertEquals(readBases, new String(block.getReadBases()));
		if(refBases == null) {
			assertNull(block.getRefBases());
		}
		else {
			assertNotNull(block.getRefBases());
			assertEquals(refBases, new String(block.getRefBases()));
		}
		assertSame(covers, block.getForm());
		assertEquals(len, block.length());
		if(matches != null) {
			assertEquals(matches.intValue(), block.getMatchCount().count);
		}
	}

	/**
	 * @return
	 */
	public static ReferenceSequenceFile buildSeqDb() {
		return ReferenceSequenceFileFactory.getReferenceSequenceFile(
				FASTA_FILE);
	}
	
	public static IndexedFastaFetcher buildFetcher()
	throws FileNotFoundException {
		return new IndexedFastaFetcher(FASTA_FILE);
	}
}
