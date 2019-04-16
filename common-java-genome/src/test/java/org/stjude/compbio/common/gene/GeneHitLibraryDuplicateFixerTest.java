package org.stjude.compbio.common.gene;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;

import org.junit.Before;
import org.junit.Test;
import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.formats.refFlat.RefFlatReader;
import org.stjude.compbio.util.io.FormatException;

public class GeneHitLibraryDuplicateFixerTest {

	@Before
	public void setUp() throws Exception {
	}

	private GeneHitLibraryDuplicateFixer getFixer(String inFilename)
	throws IOException, FormatException {
		RefFlatReader reader = new RefFlatReader(
				getClass().getResourceAsStream(inFilename));
		return new GeneHitLibraryDuplicateFixer(reader.readAllGeneHits());
		
	}
	
	@Test
	public void testFixParalogs() throws IOException, FormatException {
		GeneHitLibraryDuplicateFixer fixer = getFixer("paralog_test_in.txt");
		fixer.setParalogEnabled(true);
		GeneHitLibrary fixed = fixer.fix();
		
		List<TxHit> txHits = getSortedTxHitList(fixed);
		assertEquals(2, txHits.size());
		
		Iterator<TxHit> iterator = txHits.iterator();
		TxHit txHit;
		
		txHit = iterator.next();
		assertEquals("BMS1P1_locPar", txHit.getGeneName());
		assertEquals("NR_026566-locPar", txHit.getTxName());
		assertEquals(46737612, txHit.getTxInterval().getLeft().intValue());
		assertEquals(46762892, txHit.getTxInterval().getRight().intValue());
	
		txHit = iterator.next();
		assertEquals("BMS1P5_locPar", txHit.getGeneName());
		assertEquals("NR_003611-locPar", txHit.getTxName());
		assertEquals(48927373, txHit.getTxInterval().getLeft().intValue());
		assertEquals(48952629, txHit.getTxInterval().getRight().intValue());
	
	}
	
	@Test
	public void testOverlappingLoci() throws IOException, FormatException {
		GeneHitLibraryDuplicateFixer fixer = getFixer("loc_ovlp_test_in.txt");
		fixer.setParalogEnabled(true);
		fixer.setLocusSuffixingEnabled(true);
		fixer.setSimpleTxSuffixingEnabled(true);
		GeneHitLibrary fixed = fixer.fix();
		
		List<TxHit> txHits = getSortedTxHitList(fixed);
		assertEquals(2, txHits.size());
		
		Iterator<TxHit> iterator = txHits.iterator();
		TxHit txHit;
		
		txHit = iterator.next();
		assertEquals("DUX4L7", txHit.getGeneName());
		assertEquals("NM_001127387-1", txHit.getTxName());
	
		txHit = iterator.next();
		assertEquals("DUX4L7", txHit.getGeneName());
		assertEquals("NM_001127387-2", txHit.getTxName());
	
	}
	
	@Test
	public void testSPAG11B() throws IOException, FormatException {
		GeneHitLibraryDuplicateFixer fixer = getFixer("SPAG11B_test_in.txt");
		fixer.setParalogEnabled(true);
		fixer.setLocusSuffixingEnabled(true);
		fixer.setSimpleTxSuffixingEnabled(true);
		GeneHitLibrary fixed = fixer.fix();
		
		List<TxHit> txHits = getSortedTxHitList(fixed);
		assertEquals(8, txHits.size());
		
		Iterator<TxHit> iterator = txHits.iterator();
		TxHit txHit;
		
		txHit = iterator.next();
		assertEquals("SPAG11B_locA", txHit.getGeneName());
		assertEquals(Strand.REVERSE, txHit.getTxModel().getStrand());
		assertEquals("SPAG11B_locA", iterator.next().getGeneName());
		assertEquals("SPAG11B_locA", iterator.next().getGeneName());
		assertEquals("SPAG11B_locA", iterator.next().getGeneName());
		assertEquals("SPAG11B_locA", iterator.next().getGeneName());
		assertEquals("SPAG11B_locA", iterator.next().getGeneName());
		assertEquals("SPAG11B_locA", iterator.next().getGeneName());
		
		txHit = iterator.next();
		assertEquals("SPAG11B_locB", txHit.getGeneName());
		assertEquals("NM_058203-locB", txHit.getTxName());
		assertEquals(Strand.FORWARD, txHit.getTxModel().getStrand());
		assertEquals(7705401, txHit.getTxInterval().getLeftInt());
		assertEquals(7707807, txHit.getTxInterval().getRightInt());
	}
	
	@Test
	public void testDUX2() throws IOException, FormatException {
		GeneHitLibraryDuplicateFixer fixer = getFixer("DUX2_test_in.txt");
		fixer.setParalogEnabled(true);
		fixer.setLocusSuffixingEnabled(true);
		fixer.setSimpleTxSuffixingEnabled(true);
		GeneHitLibrary fixed = fixer.fix();
		
		List<TxHit> txHits = getSortedTxHitList(fixed);
		assertEquals(19, txHits.size());
		
		Iterator<TxHit> iterator = txHits.iterator();
		TxHit txHit;
		
		txHit = iterator.next();
		assertEquals("DUX2_locA", txHit.getGeneName());
		assertEquals(Strand.FORWARD, txHit.getTxModel().getStrand());
		assertEquals("NM_012147-locA", txHit.getTxName());
		assertEquals("NM_012147-locB", iterator.next().getTxName());
		assertEquals("NM_012147-locC", iterator.next().getTxName());
		assertEquals("NM_012147-locD", iterator.next().getTxName());
		assertEquals("NM_012147-locE", iterator.next().getTxName());
		assertEquals("NM_012147-locF", iterator.next().getTxName());
		assertEquals("NM_012147-locG", iterator.next().getTxName());
		assertEquals("NM_012147-locH", iterator.next().getTxName());
		assertEquals("NM_012147-locI", iterator.next().getTxName());
		assertEquals("NM_012147-locJ", iterator.next().getTxName());
		assertEquals("NM_012147-locK-1", iterator.next().getTxName());
		assertEquals("NM_012147-locK-2", iterator.next().getTxName());
		assertEquals("NM_012147-locL", iterator.next().getTxName());
		assertEquals("NM_012147-locM", iterator.next().getTxName());
		assertEquals("NM_012147-locN", iterator.next().getTxName());
		assertEquals("NM_012147-locO", iterator.next().getTxName());
		assertEquals("NM_012147-locP", iterator.next().getTxName());
		assertEquals("NM_012147-locQ", iterator.next().getTxName());
		assertEquals("NM_012147-locR", iterator.next().getTxName());
	}

	@Test
	public void testGt26() throws IOException, FormatException {
		GeneHitLibraryDuplicateFixer fixer = getFixer("gt26_test_in.txt");
		fixer.setParalogEnabled(true);
		fixer.setLocusSuffixingEnabled(true);
		fixer.setSimpleTxSuffixingEnabled(true);
		GeneHitLibrary fixed = fixer.fix();
		
		List<TxHit> txHits = getSortedTxHitList(fixed);
		assertEquals(144, txHits.size());
		
		Iterator<TxHit> iterator = txHits.iterator();
		TxHit txHit;
		
		txHit = iterator.next();
		assertEquals("LOC101928804_locAA", txHit.getGeneName());
		
		BufferedReader namesReader = new BufferedReader(new InputStreamReader(
				getClass().getResourceAsStream("gt26_test_out_names.txt")));
		String expName = namesReader.readLine();
		assertEquals(expName, txHit.getTxName());
		while(iterator.hasNext()) {
			txHit = iterator.next();
			expName = namesReader.readLine();
			assertNotNull("Ran out of expected names", expName);
			assertEquals(expName, txHit.getTxName());
		}
		assertNull("Too many expected names", namesReader.readLine());
	}
	
	private List<TxHit> getSortedTxHitList(GeneHitLibrary library) {
		List<TxHit> out = new LinkedList<TxHit>();
		
		List<String> geneNames = new ArrayList<String>(
				library.getAllGeneNames());
		Collections.sort(geneNames);
		
		for(String geneName: geneNames) {
			GeneHitGroup geneGrp = library.getGeneHitGroup(geneName);
			SortedMap<String, TxHitGroup> txHitGroupMap =
					new TreeMap<String, TxHitGroup>(geneGrp.getTxHitGroupMap());
			
			for(String txName: txHitGroupMap.keySet()) {
				TxHitGroup txGrp = txHitGroupMap.get(txName);
				List<TxHit> txHits = new ArrayList<TxHit>(txGrp.getAllTxHits());
				
				Collections.sort(txHits, new Comparator<TxHit>() {
					public int compare(TxHit o1, TxHit o2) {
						return o1.getTxInterval().compareTo(o2.getTxInterval());
					}});
				out.addAll(txHits);
			}
		}
		return out;
	}

}
