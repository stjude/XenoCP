package org.stjude.compbio.sam.block;

import static org.junit.Assert.*;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.stjude.compbio.sam.MdElts;

public class BlocklistTest {
	/**
	 * Reference sequence (set up one time)
	 */
	//private static IndexedFastaFetcher fastaFetcher;
	
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
		//fastaFetcher = BlockTestHelper.buildFetcher();
	}
	
	@Before
	public void setUp() throws Exception {
		helper = new BlockTestHelper();
		factory = new BlocklistFactory(false, null);
	}
	
	@Test
	public void testGet() {
		helper.setStandardRecordInfo();
		factory.setMdParsingEnabled(true);
		
		// Make the blocklist
		Blocklist blocklist = factory.newBlocklist(helper.record);

		Block block = blocklist.get(5);
		helper.assertBlock(block, new CigarElement(10, CigarOperator.M),
				20, 30, 100020, 100030,
				"TGTCCTGAAT", "TGTCTAGAAT", BlockForm.BOTH, 10, 8);
	}

	@Test
	public void testToCigar() {
		helper.setStandardRecordInfo();
		factory.setMdParsingEnabled(true);
		
		// Make the blocklist
		Blocklist blocklist = factory.newBlocklist(helper.record);
		// Get CIGAR from blocklist, and make sure it matches record
		assertEquals(helper.record.getCigarString(),
				blocklist.toCigar().toString());
	}

	@Test
	public void testToMdEltList() {
		helper.setStandardRecordInfo();
		factory.setMdParsingEnabled(true);
		
		// Make the blocklist
		Blocklist blocklist = factory.newBlocklist(helper.record);
		// Get CIGAR from blocklist, and make sure it matches record
		assertEquals(helper.record.getAttribute("MD"),
				MdElts.listToString(blocklist.toMdEltList()));
	}

	@Test
	public void testReplace() {
		helper.setStandardRecordInfo();
		factory.setMdParsingEnabled(true);
		
		// Make the blocklist
		Blocklist blocklist = factory.newBlocklist(helper.record);
	
		Block block = blocklist.get(5);
		Block next = blocklist.remove(block);
		blocklist.insert(next, 10, CigarOperator.M);
		
		assertEquals(helper.record.getCigarString(),
				blocklist.toCigar().toString());		
		assertEquals(helper.record.getAttribute("MD"),
				MdElts.listToString(blocklist.toMdEltList()));
	}

	@Test
	public void testFreezeAndTranslate() {
		helper.setStandardRecordInfo();
		factory.setMdParsingEnabled(true);
		
		// Make the blocklist
		Blocklist blocklist = factory.newBlocklist(helper.record);
		blocklist.setStickyRefBases(true);
		blocklist.translate("foo", 0);

		Block block = blocklist.get(5);
		assertEquals("foo", block.getRefInterval().getRefName());
		helper.assertBlock(block, new CigarElement(10, CigarOperator.M),
				20, 30, 20, 30,
				"TGTCCTGAAT", "TGTCTAGAAT", BlockForm.BOTH, 10, 8);
		assertEquals(helper.record.getCigarString(),
				blocklist.toCigar().toString());		
		assertEquals(helper.record.getAttribute("MD"),
				MdElts.listToString(blocklist.toMdEltList()));
	}

}
