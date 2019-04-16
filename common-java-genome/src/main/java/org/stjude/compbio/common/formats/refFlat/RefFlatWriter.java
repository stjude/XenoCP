package org.stjude.compbio.common.formats.refFlat;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.stjude.compbio.common.gene.GeneHitGroup;
import org.stjude.compbio.common.gene.GeneHitLibrary;
import org.stjude.compbio.common.gene.TxHit;
import org.stjude.compbio.common.gene.TxHitGroup;
import org.stjude.compbio.common.gene.TxModel;
import org.stjude.compbio.common.interval.LongInterval;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.StdLongInterval;
import org.stjude.compbio.util.DelimitedListBuilder;

/**
 * Class that can write refFlat-formatted files
 */
public class RefFlatWriter {
	/**
	 * Writer used to write the file
	 */
	private PrintWriter writer;
	/**
	 * Whether to close the underlying writer on close
	 */
	private boolean closeWriter = true;

	public RefFlatWriter(PrintWriter writer) {
		this.writer = writer;
	}
	
	public RefFlatWriter(File file) throws IOException {
		this(new PrintWriter(new FileWriter(file)));
	}
	
	/**
	 * Constructs a writer to System.out
	 */
	public RefFlatWriter() {
		this(new PrintWriter(System.out));
		this.closeWriter = false;
	}
	
	/**
	 * Writes an entire GeneHitLibrary
	 */
	public void write(GeneHitLibrary library) {
		for(GeneHitGroup geneGrp: library.getAllGeneHitGroups()) {
			write(geneGrp);
		}
	}

	/**
	 * Writes all records for the GeneHitGroup
	 */
	public void write(GeneHitGroup geneGrp) {
		for(TxHitGroup txGrp: geneGrp.getAllTxHitGroups()) {
			write(txGrp);
		}
	}

	/**
	 * Writes all records for the TxHitGroup
	 */
	public void write(TxHitGroup txGrp) {
		for(TxHit txHit: txGrp.getAllTxHits()) {
			write(txHit);
		}
	}

	/**
	 * Writes the TxHit
	 */
	public void write(TxHit txHit) {
		// Much is based off the TxModel
		TxModel txModel = txHit.getTxModel();
		
		// Build up left/right lists
		DelimitedListBuilder lefts = new DelimitedListBuilder(",");
		DelimitedListBuilder rights = new DelimitedListBuilder(",");
		for(LongInterval exon: txModel.getExons()) {
			lefts.append(exon.getLeft());
			rights.append(exon.getRight());
		}
		
		// Get intervals and handle null cds interval
		RefInterval txInterval = txModel.getTxInterval();
		LongInterval cdsInterval = txModel.getCdsInterval();
		if(cdsInterval == null) {
			cdsInterval = new StdLongInterval(
					txInterval.getRight(), txInterval.getRight());
		}
		
		// Write
		writer.println(DelimitedListBuilder.implode("\t", new String[] {
				txHit.getGeneName(),
				txHit.getTxName(),
				txHit.getRefName(),
				txHit.getStrand().toString(),
				String.valueOf(txInterval.getLeft()),
				String.valueOf(txInterval.getRight()),
				String.valueOf(cdsInterval.getLeft()),
				String.valueOf(cdsInterval.getRight()),
				String.valueOf(txModel.getExons().size()),
				lefts.toString(),
				rights.toString()
		}));
	}
	
	/**
	 * Closes the writer
	 */
	public void close() {
		if(closeWriter) writer.close(); else writer.flush();
	}
}
