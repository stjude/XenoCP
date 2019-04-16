package org.stjude.compbio.sam.block;

import static org.junit.Assert.*;

import java.util.ListIterator;
import java.util.NoSuchElementException;

import org.junit.Before;
import org.junit.Test;
import org.stjude.compbio.sam.block.DirectionalBlockIterator.Direction;

public class DirectionalBlockIteratorTest {
	/**
	 * The test Blocklist
	 */
	private Blocklist blocklist;
	
	/**
	 * Test data object that holds the SAMdata.record
	 */
	private BlockTestHelper helper;
	
	@Before
	public void setUp() throws Exception {
		helper = new BlockTestHelper();
		helper.setStandardRecordInfo();
		
		BlocklistFactory factory = new BlocklistFactory(false, null);
		blocklist = factory.newBlocklist(helper.record);
	}
	

	@Test
	public void testForward() {
		ListIterator<Block> iterator = new DirectionalBlockIterator(
				blocklist, Direction.FROM_LEFT_ANCHOR);
		
		// Starting spans
		assertEquals(0, blocklist.getReadLeft());
		assertEquals(50, blocklist.getReadRight());
		assertEquals(100000, blocklist.getRefLeft());
		assertEquals(100090, blocklist.getRefRight());
		
		// 10H    10M    2I   8M (3A4)
		assertTrue(iterator.hasNext());
		helper.assertBlockCigar("10H", iterator.next());
		assertTrue(iterator.hasNext());
		helper.assertBlockCigar("10M", iterator.next());
		assertTrue(iterator.hasNext());
		helper.assertBlockCigar("2I", iterator.next());
		assertTrue(iterator.hasNext());
		helper.assertBlockCigar("8M", iterator.next());
		
		// 8M -> 3M
		iterator.set(Block.newTemplate(3, "M"));
		assertEquals(0, blocklist.getReadLeft());
		assertEquals(45, blocklist.getReadRight());
		assertEquals(100000, blocklist.getRefLeft());
		assertEquals(100085, blocklist.getRefRight());
		
		// 8M -> 3M -> (nothing)
		iterator.remove();
		assertEquals(0, blocklist.getReadLeft());
		assertEquals(42, blocklist.getReadRight());
		assertEquals(100000, blocklist.getRefLeft());
		assertEquals(100082, blocklist.getRefRight());
		try {
			iterator.remove();
			fail("Remove should not have been allowed after remove.");
		} catch(IllegalStateException e) {
			// Good
		}
		
		// 8M -> 3M -> nothing -> 3= 1X 4=
		iterator.add(Block.newTemplate(3, "="));
		iterator.add(Block.newTemplate(1, "X"));
		iterator.add(Block.newTemplate(4, "EQ"));
		assertEquals(0, blocklist.getReadLeft());
		assertEquals(50, blocklist.getReadRight());
		assertEquals(100000, blocklist.getRefLeft());
		assertEquals(100090, blocklist.getRefRight());
		assertEquals("10H10M2I3=1X4=2D10M50N10M10S",
				blocklist.toCigar().toString());
		
		try {
			iterator.remove();
			fail("Remove should not have been allowed after add.");
		} catch(IllegalStateException e) {
			// Good
		}
		try {
			iterator.set(Block.newTemplate(1, "M"));
			fail("Set should not have been allowed after add.");
		} catch(IllegalStateException e) {
			// Good
		}
		
		assertTrue(iterator.hasNext());
		iterator.next();
		assertTrue(iterator.hasNext());
		iterator.next();
		assertTrue(iterator.hasNext());
		iterator.next();
		assertTrue(iterator.hasNext());
		iterator.next();
		assertTrue(iterator.hasNext());
		iterator.next();
		assertFalse(iterator.hasNext());
		try {
			iterator.next();
			fail("Should have thrown NoSuchElementException on next");
		} catch(NoSuchElementException e) {
			// Good
		}
		
	}

	@Test
	public void testReverse() {
		ListIterator<Block> iterator = new DirectionalBlockIterator(
				blocklist, Direction.FROM_RIGHT_ANCHOR);
		
		// Starting spans
		assertEquals(0, blocklist.getReadLeft());
		assertEquals(50, blocklist.getReadRight());
		assertEquals(100000, blocklist.getRefLeft());
		assertEquals(100090, blocklist.getRefRight());
		
		// 10S 10M
		assertTrue(iterator.hasNext());
		helper.assertBlockCigar("10S", iterator.next());
		assertTrue(iterator.hasNext());
		helper.assertBlockCigar("10M", iterator.next());
		assertTrue(iterator.hasPrevious());
		helper.assertBlockCigar("10M", iterator.previous());
		assertTrue(iterator.hasPrevious());
		helper.assertBlockCigar("10S", iterator.previous());
		assertFalse(iterator.hasPrevious());
		
		// 10S -> 10H
		iterator.set(Block.newTemplate(10, "H"));
		assertEquals(10, blocklist.getReadLeft());
		assertEquals(50, blocklist.getReadRight());
		assertEquals(100000, blocklist.getRefLeft());
		assertEquals(100090, blocklist.getRefRight());
		
		// 10S -> 10H -> (nothing)
		iterator.remove();
		assertEquals(10, blocklist.getReadLeft());
		assertEquals(50, blocklist.getReadRight());
		assertEquals(100000, blocklist.getRefLeft());
		assertEquals(100090, blocklist.getRefRight());
		try {
			iterator.remove();
			fail("Remove should not have been allowed after remove.");
		} catch(IllegalStateException e) {
			// Good
		}
		
		// 10S -> 10H -> (nothing) -> 5S 5M
		iterator.add(Block.newTemplate(5, "S"));
		iterator.add(Block.newTemplate(5, "M"));
		assertEquals(0, blocklist.getReadLeft());
		assertEquals(50, blocklist.getReadRight());
		assertEquals(99995, blocklist.getRefLeft());
		assertEquals(100090, blocklist.getRefRight());
		assertEquals("10H10M2I8M2D10M50N10M5M5S",
				blocklist.toCigar().toString());
		
		try {
			iterator.remove();
			fail("Remove should not have been allowed after add.");
		} catch(IllegalStateException e) {
			// Good
		}
		try {
			iterator.set(Block.newTemplate(1, "M"));
			fail("Set should not have been allowed after add.");
		} catch(IllegalStateException e) {
			// Good
		}
	}

}
