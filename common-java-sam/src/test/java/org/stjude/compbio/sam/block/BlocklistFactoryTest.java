package org.stjude.compbio.sam.block;

import static org.junit.Assert.*;

import java.util.Iterator;

import htsjdk.samtools.CigarElement;

import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.Before;

import org.stjude.compbio.common.formats.fasta.IndexedFastaFetcher;
import org.stjude.compbio.sam.MdElts;

public class BlocklistFactoryTest {
	/**
	 * Reference sequence (set up one time)
	 */
	private static IndexedFastaFetcher fastaFetcher;
	
	/**
	 * The factory, which will be created but must be configured
	 */
	private BlocklistFactory factory;
	
	/**
	 * Test data object that holds the SAMdata.record
	 */
	private BlockTestHelper helper;
	
	@BeforeClass
	public static void setUpClass() throws Exception {
		fastaFetcher = BlockTestHelper.buildFetcher();
	}
	
	@Before
	public void setUp() throws Exception {
		helper = new BlockTestHelper();
		factory = new BlocklistFactory(false, null);
	}
	
	@Test
	public void testStandardUsingSeqDbOnly() {
		helper.setStandardRecordInfo();
		
		// Seq db only
		factory.setRefFetcher(fastaFetcher);
		factory.setMdParsingEnabled(false);
		
		// Make the blocklist
		Blocklist blocklist = factory.newBlocklist(helper.record);
		
		// Make sure it corresponds to what we'd expect from the standard
		assertStandardBlocklist(blocklist);
	}
	
	@Test
	public void testStandardUsingMDOnly() {
		helper.setStandardRecordInfo();
		
		// MD only
		factory.setMdParsingEnabled(true);
		
		// Make the blocklist
		Blocklist blocklist = factory.newBlocklist(helper.record);
		
		// Make sure it corresponds to what we'd expect from the standard
		assertStandardBlocklist(blocklist);
	}
	
	private void assertStandardBlocklist(Blocklist blocklist) {
		Iterator<Block> iter = blocklist.iterator();
		Iterator<CigarElement> ceIter =
			helper.record.getCigar().getCigarElements().iterator();

		// 0 10H 
		assertTrue(iter.hasNext());
		helper.assertBlock(iter.next(), ceIter.next(), 0, 0, 100000, 100000,
				"", "", BlockForm.NEITHER, 0, null);

		// 1 10M 
		assertTrue(iter.hasNext());
		helper.assertBlock(iter.next(), ceIter.next(), 0, 10, 100000, 100010,
				"ACTAAGCACA", "ACTAAGCACA", BlockForm.BOTH, 10, 10);

		// 2 2I
		assertTrue(iter.hasNext());
		helper.assertBlock(iter.next(), ceIter.next(), 10, 12, 100010, 100010,
				"TG", "", BlockForm.READ_ONLY, 2, null);
		
		// 3 8M
		assertTrue(iter.hasNext());
		helper.assertBlock(iter.next(), ceIter.next(), 12, 20, 100010, 100018,
				"CAGCGAAT", "CAGAGAAT", BlockForm.BOTH, 8, 7);

		// 4 2D
		assertTrue(iter.hasNext());
		helper.assertBlock(iter.next(), ceIter.next(), 20, 20, 100018, 100020,
				"", "AA", BlockForm.REF_ONLY, 2, null);

		// 5 10M
		assertTrue(iter.hasNext());
		helper.assertBlock(iter.next(), ceIter.next(), 20, 30, 100020, 100030,
				"TGTCCTGAAT", "TGTCTAGAAT", BlockForm.BOTH, 10, 8);
		
		// 6 50N
		assertTrue(iter.hasNext());
		helper.assertBlock(iter.next(), ceIter.next(), 30, 30, 100030, 100080,
				"", null, BlockForm.REF_SKIP, 50, null);
		
		// 7 10M
		assertTrue(iter.hasNext());
		helper.assertBlock(iter.next(), ceIter.next(), 30, 40, 100080, 100090,
				"TGTCAAGTAA", "TGACATGTAA", BlockForm.BOTH, 10, 8);

		// 8 10S
		assertTrue(iter.hasNext());
		helper.assertBlock(iter.next(), ceIter.next(), 40, 50, 100090, 100090,
				"ACGTACGTAC", "", BlockForm.READ_ONLY, 10, null);
		
		// End
		assertFalse(iter.hasNext());
		
		// CIGAR
		assertEquals(helper.record.getCigarString(),
				blocklist.toCigar().toString());
		
		// MD
		assertEquals(helper.getMdString(),
				MdElts.listToString(blocklist.toMdEltList()));
	}
	
	/**
	 * Sets the record info to same as standard, but using = and X in CIGAR
	 */
	private void setFineStandardRecordInfo() {
		helper.setRecordInfo(
				"10H10=2I3=1X4=2D4=2X4=50N2=1X2=1X4=10S",
				"ACTAAGCACATGCAGCGAATTGTCCTGAATTGTCAAGTAAACGTACGTAC",
				"13A4^AA4TA6A2T4");
	}

	@Test
	public void testFine() {
		// Same as standard, but fine-granularity
		setFineStandardRecordInfo();
		
		// MD only (and source granularity by default)
		factory.setMdParsingEnabled(true);
		
		// Make the blocklist
		Blocklist blocklist = factory.newBlocklist(helper.record);
		
		// Spot-check last two alignment blocks (1X4= at 15,16)
		helper.assertBlock(
				blocklist.get(14), helper.getCigarElement(14),
				35, 36, 100085, 100086, "A", "T", BlockForm.BOTH, 1, 0);
		helper.assertBlock(
				blocklist.get(15), helper.getCigarElement(15),
				36, 40, 100086, 100090, "GTAA", "GTAA", BlockForm.BOTH, 4, 4);
	}
	
	@Test
	public void testDND() {
		helper.setRecordInfo("8M2D10N2D8M", "ACTAAGCATCTAGAAT", "8^CATG8");
		
		// Seq db only
		factory.setRefFetcher(fastaFetcher);
		factory.setMdParsingEnabled(false);
		
		// Convert back to MD
		Blocklist blocklist = factory.newBlocklist(helper.record);
		
		String newMd = MdElts.listToString(blocklist.toMdEltList());
		assertEquals("8^CA^TG8", newMd);
	}
	
/*
	@Test
	public void testCoarse() {
		// Same as standard, but fine-granularity
		setFineStandardRecordInfo();
		
		// MD only, coarse
		factory.setMdParsingEnabled(true);
		factory.setGranularity(CigarElementGranularity.COARSE);
		
		// Make the blocklist
		Blocklist blocklist = factory.newBlocklist(helper.record);
		
		// Make sure it corresponds to what we'd expect from the standard
		helper.resetRecord();
		helper.setStandardRecordInfo();
		assertStandardBlocklist(blocklist);
	}
	*/
}
