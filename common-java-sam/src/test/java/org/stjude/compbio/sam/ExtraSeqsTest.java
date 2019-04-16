package org.stjude.compbio.sam;

import static org.junit.Assert.*;

import java.util.Map;

import htsjdk.samtools.fastq.FastqRecord;

import org.junit.Before;
import org.junit.Test;

public class ExtraSeqsTest {

	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void testDecode() {
		Map<String, FastqRecord> decoded =
				ExtraSeqs.decode("foo|CA|CA~bar|GTT|=FH~~~baz|GCA|,/H~~");
		assertEquals(3, decoded.size());
		assertEquals("CA", decoded.get("foo").getReadString());
		assertEquals("CA", decoded.get("foo").getBaseQualityString());
		assertEquals("GTT", decoded.get("bar").getReadString());
		assertEquals("=FH", decoded.get("bar").getBaseQualityString());
		assertEquals("GCA", decoded.get("baz").getReadString());
		assertEquals(",/H", decoded.get("baz").getBaseQualityString());
	}

}
