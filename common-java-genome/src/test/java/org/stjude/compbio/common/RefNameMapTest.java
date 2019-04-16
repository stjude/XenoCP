package org.stjude.compbio.common;

import static org.junit.Assert.*;

import org.junit.Before;
import org.junit.Test;

public class RefNameMapTest {
	private RefNameMap<String> map;
	
	@Before
	public void setUp() throws Exception {
		this.map = new RefNameMap<String>();
	}

	@Test
	public void testXY() {
		map.add("1");
		assertEquals(1, map.size());
		map.add("X");
		assertEquals(2, map.size());
		map.add("2");
		assertEquals(3, map.size());
		map.add("Y");
		assertEquals(4, map.size());
		map.add("4");
		assertEquals(5, map.size());
		
		assertEquals("X", map.get("3"));
		assertEquals("X", map.get("5"));
		assertEquals("Y", map.get("6"));
		assertEquals(null, map.get("7"));
	}

}
