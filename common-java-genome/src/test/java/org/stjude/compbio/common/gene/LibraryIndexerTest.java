package org.stjude.compbio.common.gene;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.junit.Before;
import org.junit.Test;
import org.stjude.compbio.common.RefStrand;
import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.formats.refFlat.RefFlatReader;
import org.stjude.compbio.common.interval.index.RefIntervalIndex;
import org.stjude.compbio.util.io.FormatException;

public class LibraryIndexerTest {

	@Before
	public void setUp() throws Exception {
	}

	/**
	 * Creates a test RefIntervalIndex using an input file resource
	 *  
	 * @param stranded whether or not to create a stranded index
	 * @param key what type of object should be used as the interval (tx, 
	 * 				exon, cds)
	 * @param inputFile input filename, in the current package
	 * @return new RefIntervalIndex
	 * @throws IOException
	 * @throws FormatException
	 */
	private RefIntervalIndex<TxHit> createTestIndex(boolean stranded,
			LibraryIndexer.Key key, String inputFile) throws IOException,
			FormatException {
		RefFlatReader reader = new RefFlatReader(
				getClass().getResourceAsStream(inputFile));
		StdGeneHitLibrary library = reader.readAllGeneHits();
		LibraryIndexer indexer =
				new LibraryIndexer(stranded, key);
		RefIntervalIndex<TxHit> newTxHitIndex = indexer.newTxHitIndex(library);
		return newTxHitIndex;
	}

	@Test
	public void testGTF2I() throws IOException, FormatException {
		RefIntervalIndex<TxHit> index = createTestIndex(
				true, LibraryIndexer.Key.TX, "GTF2I_DNAL1_refFlat.txt");
		
		assertEquals(2, index.refNameSet().size());
		
		List<TxHit> hits;
		
		// Check wrong chr 
		hits = index.query(new RefStrand("chr8", Strand.FORWARD), 74150000);
		assertEquals(0, hits.size());
		
		// Check wrong strand 
		hits = index.query(new RefStrand("chr7", Strand.REVERSE), 74150000);
		assertEquals(0, hits.size());
		
		// Check GTF2I 
		hits = index.query(new RefStrand("chr7", Strand.FORWARD), 74150000);
		assertEquals(5, hits.size());
		for(int i = 0; i < 5; i++) {
			assertEquals("Gene name " + i, "GTF2I", hits.get(i).getGeneName());
		}
		
	}

	@Test
	public void testGTF2IFullInput() throws IOException, FormatException {
		RefIntervalIndex<TxHit> index = createTestIndex(
				true, LibraryIndexer.Key.TX, "full_nopar_refFlat.txt");
		
		List<TxHit> hits;
		
		assertEquals(32, index.refNameSet().size());
		
		// Check wrong strand 
		hits = index.query(new RefStrand("chr7", Strand.REVERSE), 74150000);
		assertEquals(0, hits.size());
		
		// Check GTF2I 
		hits = index.query(new RefStrand("chr7", Strand.FORWARD), 74150000);
		assertEquals(5, hits.size());
		for(int i = 0; i < 5; i++) {
			assertEquals("Gene name " + i, "GTF2I", hits.get(i).getGeneName());
		}
		
	}
	

}
