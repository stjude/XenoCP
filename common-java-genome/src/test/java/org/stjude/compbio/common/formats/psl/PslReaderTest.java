package org.stjude.compbio.common.formats.psl;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;
import org.stjude.compbio.util.io.IoRecord;

public class PslReaderTest {

	public PslReader openReaderOnResource(String rsrcName) throws IOException {
		return new PslReader(getClass().getResourceAsStream(rsrcName));
	}

	@Test
	public void testFiveRecord() throws IOException {
		PslReader reader = openReaderOnResource("five_record.psl");
		
		IoRecord<PslField> record = reader.readIoRecord();
		assertEquals(
				"HWI-M00634:2:000000000-A217N:1:1101:10049:21581/1",
				record.valueForKey(PslField.qName));
		assertEquals("39", record.valueForKey(PslField.matches));
		
		record = reader.readIoRecord();
		assertEquals("50", record.valueForKey(PslField.qSize));
		
		record = reader.readIoRecord();
		assertEquals("193990631,", record.valueForKey(PslField.tStarts));
		
		record = reader.readIoRecord();
		assertEquals("-", record.valueForKey(PslField.strand));
		
		record = reader.readIoRecord();
		assertEquals(
				"HWI-M00634:2:000000000-A217N:1:1101:10674:9560/1",
				record.valueForKey(PslField.qName));
		
		assertNull(reader.readIoRecord());
	}

}
