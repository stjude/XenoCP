package org.stjude.compbio.common.interval;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

public class StdLongIntervalTest {
	private LongIntervalTestData data;
	
	@Before
	public void setUp() throws Exception {
		data = new LongIntervalTestData();
	}

	@Test
	public void testCopyConstructor() {
		assertNotSame(data.main, new StdLongInterval(data.main));
		assertEquals(data.main, new StdLongInterval(data.main));
	}

	@Test
	public void testGetLeftInt() {
		assertEquals(10, data.main.getLeftInt());
	}

	@Test
	public void testGetRightInt() {
		assertEquals(20, data.main.getRightInt());
	}

	@Test
	public void testLength() {
		assertEquals(10, data.main.length());
	}

	@Test
	public void testGetMidpoint() {
		assertEquals(15, data.main.getMidpoint());
	}

	@Test
	public void testOverlapLength() {
		assertEquals(0, data.main.overlapLength(data.abutL));
		assertEquals(0, data.main.overlapLength(data.abutR));
		assertEquals(5, data.main.overlapLength(data.ovlp5L));
		assertEquals(5, data.main.overlapLength(data.ovlp5R));
	}

	@Test
	public void testShift() {
		assertEquals(data.ovlp5L, data.main.shift(-5));
		assertEquals(data.ovlp5R, data.main.shift(5));
	}

	@Test
	public void testIsEndpoint() {
		assertTrue(data.main.isEndpoint((long)10));
		assertTrue(data.main.isEndpoint((long)20));
		assertFalse(data.main.isEndpoint((long)30));
	}

	@Test
	public void testContainsClosed() {
		assertTrue(data.main.containsClosed((long)10));
		assertTrue(data.main.containsClosed((long)15));
		assertTrue(data.main.containsClosed((long)20));
		assertFalse(data.main.containsClosed((long)25));
	}

	@Test
	public void testContainsOpen() {
		assertFalse(data.main.containsOpen((long)10));
		assertTrue(data.main.containsOpen((long)15));
		assertFalse(data.main.containsOpen((long)20));
		assertFalse(data.main.containsOpen((long)25));
	}

	@Test
	public void testContains() {
		assertTrue(data.main.contains(data.main));
		assertTrue(data.main.contains(data.inside));
		assertFalse(data.inside.contains(data.main));
		assertFalse(data.main.contains(data.disjL));
		assertFalse(data.main.contains(data.abutL));
		assertFalse(data.main.contains(data.ovlp5L));
		assertFalse(data.main.contains(data.disjR));
		assertFalse(data.main.contains(data.abutR));
		assertFalse(data.main.contains(data.ovlp5R));
	}

	@Test
	public void testOverlaps() {
		assertTrue(data.main.overlaps(data.main));
		assertTrue(data.main.overlaps(data.inside));
		assertFalse(data.main.overlaps(data.disjL));
		assertFalse(data.main.overlaps(data.abutL));
		assertTrue(data.main.overlaps(data.ovlp5L));
		assertFalse(data.main.overlaps(data.disjR));
		assertFalse(data.main.overlaps(data.abutR));
		assertTrue(data.main.overlaps(data.ovlp5R));
	}

	@Test
	public void testOverlapsOrAbuts() {
		assertTrue(data.main.overlapsOrAbuts(data.main));
		assertTrue(data.main.overlapsOrAbuts(data.inside));
		assertFalse(data.main.overlapsOrAbuts(data.disjL));
		assertTrue(data.main.overlapsOrAbuts(data.abutL));
		assertTrue(data.main.overlapsOrAbuts(data.ovlp5L));
		assertFalse(data.main.overlapsOrAbuts(data.disjR));
		assertTrue(data.main.overlapsOrAbuts(data.abutR));
		assertTrue(data.main.overlapsOrAbuts(data.ovlp5R));
	}

	@Test
	public void testIntersect() {
		assertEquals(new StdLongInterval(10, 10), data.main.intersect(data.abutL));
		assertEquals(new StdLongInterval(10, 15), data.main.intersect(data.ovlp5L));
		assertEquals(new StdLongInterval(20, 20), data.main.intersect(data.abutR));
		assertEquals(new StdLongInterval(15, 20), data.main.intersect(data.ovlp5R));
	}

	@Test
	public void testUnion() {
		assertEquals(new StdLongInterval(4, 20), data.main.union(data.abutL));
		assertEquals(new StdLongInterval(5, 20), data.main.union(data.ovlp5L));
		assertEquals(new StdLongInterval(10, 26), data.main.union(data.abutR));
		assertEquals(new StdLongInterval(10, 25), data.main.union(data.ovlp5R));
	}

	@Test
	public void testCompareTo() {
		assertEquals(0, data.main.compareTo(data.main));
		assertTrue(data.main.compareTo(data.ovlp5R) < 0);
		assertTrue(data.main.compareTo(data.ovlp5L) > 0);
		assertTrue(data.main.compareTo(data.main.union(data.ovlp5R)) < 0);
	}

}
