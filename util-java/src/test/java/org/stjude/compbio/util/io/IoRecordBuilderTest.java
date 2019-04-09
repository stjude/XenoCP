package org.stjude.compbio.util.io;

import static org.junit.Assert.*;

import java.util.Iterator;
import java.util.SortedMap;
import java.util.TreeMap;

import org.junit.Before;
import org.junit.Test;

public class IoRecordBuilderTest {
	
	private enum Key { X, Y }
	
	private IoFormat<Key> format;
	private IoField<Key> fieldX, fieldY;
	private IoRecord<Key> recA, recB;
	private IoRecordBuilder<Key> builder;

	@Before
	public void setUp() throws Exception {
		format = new IoFormat<Key>();
		format.addField(fieldX = new IoField<Key>("x", 0, Key.X));
		format.addField(fieldY = new IoField<Key>("y", 1, Key.Y));
		
		SortedMap<IoField<Key>, String> data;
		
		data = new TreeMap<IoField<Key>, String>();
		data.put(fieldX, "xA");
		data.put(fieldY, "yA");
		recA = new IoRecord<Key>(format, data);
		
		data = new TreeMap<IoField<Key>, String>();
		data.put(fieldX, "xB");
		data.put(fieldY, "yB");
		recB = new IoRecord<Key>(format, data);
	}
	
	/**
	 * Set this.builder to a new builder initialized using recA
	 */
	private void initBuilderOnA() {
		builder = new IoRecordBuilder<Key>(recA);
	}
	
	/**
	 * Assert that getRecord() on the builder returns a record whose data is
	 * equal to recB
	 */
	private void assertRecB() {
		IoRecord<Key> record = builder.getRecord();
		assertNotNull(record);
		assertNotNull(record.getData());
		assertEquals(recB.getData(), record.getData());
	}

	@Test
	public void testIoRecordBuilderIoFormatOfK() {
		new IoRecordBuilder<Key>(format);
		// If no exception, then we're good
	}

	@Test
	public void testIoRecordBuilderIoRecordOfK() {
		initBuilderOnA();
		// If no exception, then we're good
	}

	@Test
	public void testSetValueForField() {
		initBuilderOnA();
		builder.setValueForField(fieldX, "xB");
		builder.setValueForField(fieldY, "yB");
		assertRecB();
	}

	@Test
	public void testSetValueForName() {
		initBuilderOnA();
		builder.setValueForName("x", "xB");
		builder.setValueForName("y", "yB");
		assertRecB();
	}

	@Test
	public void testSetValueForKey() {
		initBuilderOnA();
		builder.setValueForKey(Key.X, "xB");
		builder.setValueForKey(Key.Y, "yB");
		assertRecB();
	}

	@Test
	public void testSetValueForOrdinal() {
		initBuilderOnA();
		builder.setValueForOrdinal(0, "xB");
		builder.setValueForOrdinal(1, "yB");
		assertRecB();
	}

	@Test
	public void testValues() {
		initBuilderOnA();
		Iterator<String> iterator = builder.values().iterator();
		
		assertTrue(iterator.hasNext());
		assertEquals("xA", iterator.next());
		assertTrue(iterator.hasNext());
		assertEquals("yA", iterator.next());
	}

}
