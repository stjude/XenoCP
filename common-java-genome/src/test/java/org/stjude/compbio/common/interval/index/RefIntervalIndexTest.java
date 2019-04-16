package org.stjude.compbio.common.interval.index;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;
import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.StdRefInterval;

public class RefIntervalIndexTest {
	private RefIntervalIndex<String> unstranded;
	private RefIntervalIndex<String> stranded;
	
	private static final String[] CHRS = { "1", "2" };
	
	private RefInterval ivl(String chr, Strand strand) {
		return new StdRefInterval(chr, strand, 100, 200);
	}
	
	@Before
	public void setUp() throws Exception {
		unstranded = new RefIntervalIndex<String>(false);
		stranded = new RefIntervalIndex<String>(true);
		
		for(String chr: CHRS) {
			for(Strand strand: Strand.values()) {
				RefInterval interval = ivl(chr, strand);
				unstranded.add(interval, interval.toString());
				stranded.add(interval, interval.toString());
			}
		}
	}

	@Test
	public void testGetIndex() {
		RefInterval refStranded;
		
		refStranded = ivl("1", Strand.FORWARD);
		assertEquals(3, unstranded.getIndex(refStranded).size());
		assertEquals(1, stranded.getIndex(refStranded).size());
		
		refStranded = ivl("1", Strand.NA);
		assertEquals(3, unstranded.getIndex(refStranded).size());
		assertEquals(1, stranded.getIndex(refStranded).size());
		
		refStranded = ivl("chr1", Strand.NA);
		assertEquals(3, unstranded.getIndex(refStranded).size());
		assertEquals(1, stranded.getIndex(refStranded).size());
		
		refStranded = ivl("chr2", Strand.NA);
		assertEquals(3, unstranded.getIndex(refStranded).size());
		assertEquals(1, stranded.getIndex(refStranded).size());
		
		refStranded = ivl("3", Strand.NA);
		assertEquals(0, unstranded.getIndex(refStranded).size());
		assertEquals(0, stranded.getIndex(refStranded).size());
	}

	@Test
	public void testQuery() {
		RefInterval query;
		
		query = ivl("1", Strand.FORWARD);
		assertEquals(3, unstranded.query(query).size());
		assertEquals(1, stranded.query(query).size());

		query = ivl("4", Strand.FORWARD);
		assertEquals(0, unstranded.query(query).size());
		assertEquals(0, stranded.query(query).size());
	}

	@Test
	public void testRemove() {
		RefInterval interval = ivl("2", Strand.REVERSE);
		RefInterval i2 = ivl("2", Strand.FORWARD);
		assertTrue(unstranded.remove(interval, interval.toString()));
		assertTrue(stranded.remove(interval, interval.toString()));
		
		assertEquals(5, unstranded.size());
		assertEquals(5, stranded.size());
		
		assertEquals(2, unstranded.query(interval).size());
		assertEquals(0, stranded.query(interval).size());
		assertEquals(2, unstranded.query(i2).size());
		assertEquals(1, stranded.query(i2).size());
	}

}
