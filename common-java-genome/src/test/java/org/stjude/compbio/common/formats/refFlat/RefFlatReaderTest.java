package org.stjude.compbio.common.formats.refFlat;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.Collection;

import org.junit.Test;
import org.stjude.compbio.common.gene.GeneModel;
import org.stjude.compbio.common.gene.Spliceform;
import org.stjude.compbio.common.gene.StdGeneModelLibrary;
import org.stjude.compbio.common.gene.TxModel;
import org.stjude.compbio.common.interval.SubRefInterval;
import org.stjude.compbio.util.io.FormatException;

public class RefFlatReaderTest {
	@Test
	public void testReadAllGeneModels() throws IOException, FormatException {
		// Read in to gene model library
		RefFlatReader reader = get5RecordReader();
		StdGeneModelLibrary library = reader.readAllGeneModels();
		Collection<GeneModel> geneModels = library.getGeneModels();
		assertEquals(2, geneModels.size());
		
		GeneModel gene;
		Collection<Spliceform> spliceforms;
		Spliceform spliceform;
		TxModel txModel;
		
		// Check LOC100996291
		gene = library.getGeneModel("LOC100996291");
		assertNotNull(gene);
		spliceforms = gene.getSpliceforms();
		assertEquals(1, spliceforms.size());
		spliceform = spliceforms.iterator().next();
		assertEquals("NR_073178", spliceform.getName());
		txModel = spliceform.getTxModel();
		assertNull(txModel.getCdsInterval());
		
		// Check GDA
		gene = library.getGeneModel("GDA");
		assertNotNull(gene);
		assertEquals(4, gene.getSpliceforms().size());
		assertEquals("GDA@chr9:74729510-74867140(+)",
				gene.getInterval().toString());
		spliceform = gene.getSpliceform("NM_004293");
		assertEquals("GDA", spliceform.getGeneName());
		txModel = spliceform.getTxModel();
		assertEquals("chr9:74764266-74867140(+)",
				txModel.getTxInterval().toString());
		assertEquals("chr9:74764475-74863258(+)",
				txModel.getCdsInterval().toString());
		assertEquals(new Integer(5),
				((SubRefInterval)txModel.getExon(5)).getOrdinal());
		assertEquals(new Integer(0), txModel.getCdsStartExon());
		assertEquals(new Integer(13), txModel.getCdsEndExon());
		
	}

	public RefFlatReader get5RecordReader() {
		return new RefFlatReader(
				RefFlatReaderTest.class.getResourceAsStream(
						"refFlat-5records.txt"));
	}

}
