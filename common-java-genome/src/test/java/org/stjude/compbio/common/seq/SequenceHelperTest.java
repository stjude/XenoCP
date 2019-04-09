package org.stjude.compbio.common.seq;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

public class SequenceHelperTest {

	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void testEvenLengthRcAndCase() {
		assertEquals("nacgtNACGT",
				SequenceHelper.reverseComplement("ACGTNacgtn"));
		
	}

	@Test
	public void testOddLengthRc() {
		assertEquals("CGT",
				SequenceHelper.reverseComplement("ACG"));
	}

	@Test
	public void testTranslationAllDefaults() {
		String nuc = "TTTTTCTTATTGCTTCTCCTACTGATTATCATAATGGTTGTCGTAGTGTCTTCCTCATCGCCTCCCCCACCGACTACCACAACGGCTGCCGCAGCGTATTACTAATAGCATCACCAACAGAATAACAAAAAGGATGACGAAGAGTGTTGCTGATGGCGTCGCCGACGGAGTAGCAGAAGGGGTGGCGGAGGG";
		String expProt = "FFLLLLLLIIIMVVVVSSSSPPPPTTTTAAAAYY**HHQQNNKKDDEECC*WRRRRSSRRGGGG";
		assertEquals(expProt, SequenceHelper.newSeqTranslator().translate(nuc));
	}

	@Test
	public void testTranslationAllDefaultsAppendExtra() {
		String nuc = "TTTTTCTTATTGCTTCT";
		String expProt = "FFLLLct";
		assertEquals(expProt, SequenceHelper.newSeqTranslator(false, false, true).translate(nuc));
	}

	@Test
	public void testTranslationAllDefaultsStrictExtra() {
		String nuc = "TTTTTCTTATTGCTTCT";
		try {
			SequenceHelper.newSeqTranslator(false, true, true).translate(nuc);
			fail("Should have thrown exception in strict mode");
		} catch(IllegalArgumentException e) {
			// Pass
		}
	}
}
