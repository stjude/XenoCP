package org.stjude.compbio.mapping.standard;

import static org.junit.Assert.*;

import java.io.File;

import org.junit.Before;
import org.junit.Test;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;

public class CleanMappedSamTest {
	
	SAMFileReader reader;
	CleanMappedSam worker;
	
	@Before
	public void setUp() throws Exception {
		reader = new SAMFileReader(
				CleanMappedSamTest.class.getResourceAsStream("test_fix_id.sam"));
		worker = new CleanMappedSam(reader, false, true, true, false, false, false);
		worker.setRefFasta(new File(
				CleanMappedSamTest.class.getResource("test_fix_id_chr1.fa").getPath()));
	}

	private String getCompString(SAMRecord record) {
		return "" + record.getAlignmentStart() + "," + record.getCigarString();
	}
	
	@Test
	public void testFixId() {
		for(SAMRecord record: reader) {
			// Get expected values from comment
			String expected = (String)record.getAttribute("CO");
			System.err.println("Test record " + getCompString(record) + " -> " + expected);
			
			// Clean it
			worker.clean(record);
			
			// Check expected values
			assertEquals(expected, getCompString(record));
		}
	}

}
