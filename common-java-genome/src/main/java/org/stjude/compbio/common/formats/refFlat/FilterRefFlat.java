package org.stjude.compbio.common.formats.refFlat;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.Options;
import org.stjude.compbio.common.GenomeMainHelper;
import org.stjude.compbio.common.Strand;
import org.stjude.compbio.common.gene.GeneHitLibrary;
import org.stjude.compbio.common.gene.LibraryIndexer;
import org.stjude.compbio.common.gene.TxHit;
import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.StdRefInterval;
import org.stjude.compbio.common.interval.index.RefIntervalIndex;
import org.stjude.compbio.util.MainHelper;
import org.stjude.compbio.util.io.FormatException;

public class FilterRefFlat {

	/**
	 * Input library
	 */
	private GeneHitLibrary library;
	
	/**
	 * Intervals to query by
	 */
	private List<RefInterval> intervals;
	
	public FilterRefFlat(GeneHitLibrary library) {
		this.library = library;
	}

	public List<RefInterval> getIntervals() {
		return intervals;
	}

	public void setIntervals(List<RefInterval> intervals) {
		this.intervals = intervals;
	}

	/**
	 * Sets query intervals from an intervals spec of the format
	 * chr:left-right,chr:left-right...
	 */
	public void setIntervals(String intervalsSpec) {
		Pattern pattern =
				Pattern.compile("^([A-Za-z0-9_]+):([0-9]+)(-([0-9]+))?$");
		
		String[] intervalSpecs = intervalsSpec.split(",");
		this.intervals = new ArrayList<RefInterval>(intervalSpecs.length); 
		for(String spec: intervalSpecs) {
			Matcher matcher = pattern.matcher(intervalsSpec);
			if(!matcher.matches()) {
				throw new IllegalArgumentException(
					"Could not use intervals " + intervalsSpec + " because "
					+ "the following interval could not be parsed: " + spec);
			}
			
			String refName = matcher.group(1);
			long left = Long.parseLong(matcher.group(2));
			long right = (matcher.groupCount() > 3)
					? Long.parseLong(matcher.group(4)) : left;
					
			intervals.add(new StdRefInterval(refName, Strand.NA, left, right));
		}
	}
	
	/**
	 * Performs the filtering, writing output to the given writer.
	 */
	public void filter(RefFlatWriter writer) {
		// TODO when there are more filters, we can't assume that there are intervals

		// Output to write
		Set<TxHit> output = new HashSet<TxHit>();
		
		// Build unstranded tx-level index
		LibraryIndexer indexer =
				new LibraryIndexer(false, LibraryIndexer.Key.TX);
		RefIntervalIndex<TxHit> index = indexer.newTxHitIndex(library);
		
		// Query each interval, adding to the output set
		for(RefInterval interval: intervals) {
			List<TxHit> txHits = index.query(interval);
			output.addAll(txHits);
		}
		
		// Write output
		for(TxHit txHit: output) {		
			writer.write(txHit);
		}
	}

	/**
	 * Helper for the main method
	 */
	private static GenomeMainHelper main =
			new GenomeMainHelper(getMainOptions());

	private static void usage(int exitcode) {
		main.usage(
				"java " + FilterRefFlat.class.getName() + " [OPTION]...\n\n" +
				"Filters a refFlat file, writing output in refFlat format.\n\n"
				);
		System.exit(exitcode);
	}

	public static Options getMainOptions() {
		Options options = new Options();
		options.addOption("o", true, "output file if not stdout");
		options.addOption("I", true, 
				"query using interval(s) chr:left-right,chr:left-right...");
		return options;
	}

	public static void main(String[] args) throws IOException {
		main.addGeneAnnoOption();
		if(args.length == 0) usage(0);
		main.setArgs(args);

		// Parse command line
		main.setStrict(false);
		File geneAnnoFile = main.getGeneAnnoFile();
		File outputFile = main.getOutputFile("o");
		String intervalsSpec = main.getArgString("I");
		
		// Read the refFlat file into memory
		GeneHitLibrary library;
		try {
			if(geneAnnoFile != null) {
				System.err.println(
						"Reading gene library from refFlat file " + geneAnnoFile);
				library = new RefFlatReader(geneAnnoFile).readAllGeneHits();
			}
			else {
				System.err.println("Reading gene library from stdin");
				library = new RefFlatReader(System.in).readAllGeneHits();
			}
		} catch (FormatException e) {
			throw new MainHelper.CommandLineException(e);
		}
		
		// Open output
		RefFlatWriter writer;
		if(outputFile != null) {
			try {
				writer = new RefFlatWriter(outputFile);
			}
			catch(IOException e) {
				throw new MainHelper.CommandLineException(e);
			}
		}
		else {
			writer = new RefFlatWriter();
		}
		
		// Create the worker
		FilterRefFlat worker = new FilterRefFlat(library);
		worker.setIntervals(intervalsSpec);

		try {
			worker.filter(writer);
		} finally {
			// Close writer
			writer.close();
		}
	}


}
