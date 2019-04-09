package org.stjude.compbio.common.gene;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.Set;
import java.util.TreeMap;

import org.stjude.compbio.common.interval.RefInterval;
import org.stjude.compbio.common.interval.index.RefIntervalIndex;
import org.stjude.compbio.util.HashMultimap;
import org.stjude.compbio.util.LinkedHashSetQueue;
import org.stjude.compbio.util.Multimap;

/**
 * Fixes problems with duplicate entries in a GeneHitLibrary and creates a new
 * modified version.
 */
public class GeneHitLibraryDuplicateFixer {
	/**
	 * Custom filter
	 */
	public static interface Filter {
		/**
		 * Returns true if the TxHit passes the filter and should be copied
		 */
		boolean isAcceptable(TxHit txHit);
	}
	
	/**
	 * Extra filters
	 */
	private Filter filter;
	/**
	 * GeneHitLibrary being tweaked; it will not be modified
	 */
	private GeneHitLibrary src;
	/**
	 * Whether or not to do paralog processing
	 */
	private boolean paralogEnabled;
	/**
	 * Whether or not to do locus-based suffixing
	 */
	private boolean locusSuffixingEnabled;
	/**
	 * Whether or not to do simple tx suffixing
	 */
	private boolean simpleTxSuffixingEnabled;
	
	/**
	 * Index on the input GeneHitLibrary
	 */
	private RefIntervalIndex<TxHit> srcTxIndex;
	
	/**
	 * Output library that is being built
	 */
	private StdGeneHitLibrary dest;
	/**
	 * GeneHitGroups remaining to be processed
	 */
	private LinkedHashSetQueue<GeneHitGroup> queue;
	/**
	 * Current pile map
	 */
	private SortedMap<RefInterval, Pile> piles;
	
	/**
	 * Constructs a GeneHitLibraryDuplicateFixer with all options
	 * 
	 * @param src source library
	 * @param paralogEnabled whether to perform paralog processing
	 * @param locusSuffixingEnabled whether to do locus-based suffixing
	 * @param simpleTxSuffixingEnabled whether to do simple numeric
	 *      suffixing of duplicate tx names
	 */
	public GeneHitLibraryDuplicateFixer(GeneHitLibrary src,
			boolean paralogEnabled, boolean locusSuffixingEnabled,
			boolean simpleTxSuffixingEnabled) {
		this.filter = null;
		this.src = src;
		this.paralogEnabled = paralogEnabled;
		this.locusSuffixingEnabled = locusSuffixingEnabled;
		this.simpleTxSuffixingEnabled = simpleTxSuffixingEnabled;
		this.srcTxIndex = null;
	}
	
	/**
	 * Constructs a fixer with all suffixing disabled.
	 */
	public GeneHitLibraryDuplicateFixer(GeneHitLibrary src) {
		this(src, false, false, false);
	}
	
	/**
	 * Gets the current filter
	 */
	public Filter getFilter() {
		return filter;
	}

	/**
	 * Sets the filter to use
	 * @param filter new filter, or null for no filtering
	 */
	public void setFilter(Filter filter) {
		this.filter = filter;
	}

	/**
	 * Returns true if a record passes filters
	 */
	private boolean isAcceptable(TxHit txHit) {
		return (filter == null) ? true : filter.isAcceptable(txHit);
	}
	
	public boolean isParalogEnabled() {
		return paralogEnabled;
	}

	public void setParalogEnabled(boolean paralogEnabled) {
		this.paralogEnabled = paralogEnabled;
	}

	public boolean isLocusSuffixingEnabled() {
		return locusSuffixingEnabled;
	}

	public void setLocusSuffixingEnabled(boolean locusSuffixingEnabled) {
		this.locusSuffixingEnabled = locusSuffixingEnabled;
	}
	
	public boolean isSimpleTxSuffixingEnabled() {
		return simpleTxSuffixingEnabled;
	}

	public void setSimpleTxSuffixingEnabled(boolean simpleTxSuffixingEnabled) {
		this.simpleTxSuffixingEnabled = simpleTxSuffixingEnabled;
	}

	public GeneHitLibrary fix() {
		// Index the library, if paralogs are enabled
		if(paralogEnabled) {
			LibraryIndexer indexer = new LibraryIndexer(
					true, LibraryIndexer.Key.TX);
			srcTxIndex = indexer.newTxHitIndex(src);
		}
		
		// Create the empty output library
		dest = new StdGeneHitLibrary();
		
		// Load all genes into a queue for processing
		queue = new LinkedHashSetQueue<GeneHitGroup>(
				src.getAllGeneHitGroups());
		
		// Iterate over it until complete
		while(!queue.isEmpty()) {
			GeneHitGroup geneGrp = queue.poll();
			
			// Check for the paralog case first, since this supercedes all else
			if(paralogEnabled) {
				// If it's a paralog, then fix and move on to next gene
				if(checkAndFixParalog(geneGrp)) continue;
			}
			
			// Perform suffixing as necessary
			// Locus suffixing will do simple as well if appropriate
			if(locusSuffixingEnabled) {
				fixByLocusSuffixing(geneGrp);
			}
			// Simple suffixing will only suffix if enabled
			else {
				fixByOptionalSimpleSuffixing(geneGrp);
			}
		}
		return dest;
	}

	/**
	 * Checks if a GeneHitGroup is a paralog, and if so, fixes it.
	 * 
	 * The paralog case means that there are n genes in n loci with the same
	 * set of transcripts for each gene.
	 */
	private boolean checkAndFixParalog(GeneHitGroup geneGrp) {
		// Make piles without filtering, as filtering before we do the detection
		// could lead to false negatives.
		piles = makePiles(geneGrp, false);
		
		// If there is only one pile, then this is not a paralog
		if(piles.size() < 2) return false;
		
		// Pick a single interval to use to find the other genes
		RefInterval qryInterval = piles.keySet().iterator().next();
		
		// Check other genes
		
		// We'll keep track of which genes have been checked, by gene name, and
		// also a sorted map of the other piles.  This map is sorted by the key,
		// which is the pile's RefInterval, so that there is a stable,
		// predictable sort order across multiple runs.  To initialize, put in
		// the gene that we started off with.
		Set<String> checked = new HashSet<String>();
		checked.add(geneGrp.getName());
		SortedMap<String, Map<RefInterval, Pile>> matches =
				new TreeMap<String, Map<RefInterval, Pile>>();
		matches.put(geneGrp.getName(), piles);
		
		// Now, look for others
		List<TxHit> txHits = srcTxIndex.query(qryInterval);
		for(TxHit txHit: txHits) {
			// Skip if we've done this gene already
			String geneName = txHit.getGeneName();
			if(checked.contains(geneName)) continue;
			checked.add(geneName);
			
			// Make piles for the other gene, and look for match with the
			// original.  Pile is set up so that two piles are equal if they
			// contain the same models, so we can use pile equality to do the
			// test
			GeneHitGroup otherGeneGrp = txHit.getTxHitGroup().getGeneHitGroup();
			Map<RefInterval, Pile> otherPiles = makePiles(otherGeneGrp, false);
			if(piles.equals(otherPiles)) {
				matches.put(geneName, otherPiles);
			}
		}
		
		// Finally, see if the number of matches is the same as the number of
		// piles, and if so, we have the paralog case
		if(matches.size() != piles.size()) return false;
		
		// Choose one pile for each gene in a predictable way, and add to dest
		// Note: we are using SortedMaps, so keySet() will be in stable order
		List<String> geneNames =
				new ArrayList<String>(matches.keySet());
		List<RefInterval> pileIntervals =
				new ArrayList<RefInterval>(piles.keySet());
		// Since each gene exists at each location, you can imagine it as a
		// square table, and we are iterating along the diagonal.
		for(int i = 0; i < geneNames.size(); i++) {
			Pile pile = matches.get(geneNames.get(i)).get(pileIntervals.get(i));
			for(TxHit txHit: pile.txHits.values()) {
				// Remove from queue (yes, this will try to remove the same one
				// many times, but it won't hurt)
				queue.remove(txHit.getTxHitGroup().getGeneHitGroup());
				// Copy to dest
				checkAndWrite(txHit, txHit.getGeneName() + "_locPar",
						txHit.getTxName() + "-locPar");
			}
		}
		
		return true;
	}
	
	/**
	 * Used by fixByLocusSuffixing
	 */
	private static class SuffixCtr {
		private final int degree;
		private int ctr = 0;
		
		/**
		 * Returns a new SuffixCtr if numPiles > 1, otherwise null
		 * @param numPiles number of piles
		 */
		public static SuffixCtr newIfApplicable(int numPiles) {
			return (numPiles > 1) ? new SuffixCtr(numPiles) : null;
		}
		
		/**
		 * @param numPiles number of piles; must be > 1
		 */
		public SuffixCtr(int numPiles) {
			this.degree = (int)
					Math.floor((Math.log(numPiles - 1) / Math.log(26))) + 1;
		}
		
		/**
		 * Increments the counter
		 */
		public void increment() {
			++ctr;
		}
		
		/**
		 * Returns the suffix for the current counter
		 * @return uppercase alphabetic suffix
		 */
		public String getSuffix() {
			int remaining = ctr;
			char chars[] = new char[degree];
			for(int i = 0; i < degree; i++) {
				chars[degree - 1 - i] = (char)('A' + remaining % 26);
				remaining = remaining / 26;
			}
			return new String(chars);
		}
	}

	/**
	 * Fixes duplicate gene and tx names by grouping by locus.
	 */
	private void fixByLocusSuffixing(GeneHitGroup geneGrp) {
		// Make piles, using the filter
		piles = makePiles(geneGrp, true);
		
		SuffixCtr suffixCtr = SuffixCtr.newIfApplicable(piles.size());
		for(Pile pile: piles.values()) {
			Buffer buffer = new Buffer();
			for(TxHit txHit: pile.txHits.values()) {
				String geneName = txHit.getGeneName();
				String txName = txHit.getTxName();
				if(suffixCtr != null) {
					geneName += "_loc" + suffixCtr.getSuffix();
					txName += "-loc" + suffixCtr.getSuffix();
				}
				buffer.add(txHit, geneName, txName);
			}
			fixBySimpleSuffixing(buffer);
			if(suffixCtr != null) suffixCtr.increment();
		}
	}

	/**
	 * Buffer for simple suffixing
	 */
	private class Buffer {
		private class Record {
			public String geneName;
			public TxHit txHit;
			public Record(String geneName, TxHit txHit) {
				this.geneName = geneName;
				this.txHit = txHit;
			}
		}
		public Multimap<String, Record, List<Record>> map =
				HashMultimap.newListInstance();
		public void add(TxHit txHit, String geneName, String txName) {
			map.add(txName, new Record(geneName, txHit));
		}
	}
	
	/**
	 * Fixes tx names only by appending 1, 2, 3, ...
	 */
	private void fixBySimpleSuffixing(Buffer buffer) {
		for(String txName: buffer.map.keySet()) {
			List<Buffer.Record> records = buffer.map.get(txName);
			
			// Write, with suffix if necessary
			boolean doSuffix = simpleTxSuffixingEnabled && (records.size() > 1);
			int ctr = 0;
			for(Buffer.Record record: records) {
				String useTxName = txName;
				if(doSuffix) useTxName += "-" + (++ctr);
				checkAndWrite(record.txHit, record.geneName, useTxName);
			}
		}
	}

	/**
	 * Filters and writes everything for a GeneHitGroup, with simple suffixing
	 * if appropriate
	 * 
	 * @param geneGrp GeneHitGroup to write
	 */
	private void fixByOptionalSimpleSuffixing(GeneHitGroup geneGrp) {
		// Filtering is done here, suffixing is done by delegate method
		for(TxHitGroup txGrp: geneGrp.getAllTxHitGroups()) {
			// Perform filtering and add to buffer
			Buffer buffer = new Buffer();
			for(TxHit txHit: txGrp.getAllTxHits()) {
				if(isAcceptable(txHit)) {
					buffer.add(txHit, txHit.getGeneName(), txHit.getTxName());
				}
			}
			// Delegate
			fixBySimpleSuffixing(buffer);
		}
	}

	/**
	 * Checks the filter, and writes the tx hit if it passes
	 * 
	 * @param txHit TxHit to write
	 * @param geneName gene name to use
	 * @param txName tx name to use
	 */
	private void checkAndWrite(TxHit txHit, String geneName, String txName) {
		if(!isAcceptable(txHit)) return;
		dest.addTxHit(geneName, txName, txHit.getTxModel());
	}

	private SortedMap<RefInterval, Pile> makePiles(
			GeneHitGroup geneGrp, boolean applyFilter) {
		// Make index
		RefIntervalIndex<Pile> idx = new RefIntervalIndex<Pile>(true);
		// Loop over hits
		for(TxHitGroup txGrp: geneGrp.getAllTxHitGroups()) {
			for(TxHit txHit: txGrp.getAllTxHits()) {
				// Filter
				if(applyFilter && !isAcceptable(txHit)) continue;
				// Look for overlapping piles, and create/add/combine as necessary
				List<Pile> hits = idx.query(txHit.getTxInterval());
				// New pile
				if(hits.isEmpty()) {
					idx.add(txHit.getTxInterval(), new Pile(txHit));
				}
				// Merge/Add to pile
				else {
					Pile merged = null;
					for(Pile pile: hits) {
						// Remove existing entry
						idx.removeAllEntries(idx.getDataEntries(pile));
						// On first iteration, choose this as THE pile, and
						// add the new entry
						if(merged == null) {
							merged = pile;
							merged.add(txHit);
						}
						// Merge subsequent piles into the first
						else {
							merged.addAll(pile);
						}
					}
					// Add to the index
					idx.add(merged.interval, merged);						
				}
			}
		}
		
		// Return piles
		SortedMap<RefInterval, Pile> piles = new TreeMap<RefInterval, Pile>();
		for(Pile pile: idx.dataSet()) {
			piles.put(pile.interval, pile);
		}
		return piles;
	}
	
	private class Pile {
		public RefInterval interval;
		public Map<TxModel, TxHit> txHits;
		
		public Pile(TxHit txHit) {
			this.interval = txHit.getTxInterval();
			this.txHits = new HashMap<TxModel, TxHit>();
			this.txHits.put(txHit.getTxModel(), txHit);
		}

		public void add(TxHit txHit) {
			interval = interval.union(txHit.getTxInterval());
			txHits.put(txHit.getTxModel(), txHit);
		}
		
		public void addAll(Pile pile) {
			for(TxHit txHit: pile.txHits.values()) add(txHit);
		}
		
		public int hashCode() {
			return interval.hashCode();
		}
		
		public boolean equals(Object obj) {
			if(obj == this) return true;
			if(obj == null) return false;
			if(!(obj instanceof Pile)) return false;
			
			Pile other = (Pile)obj;
			return interval.equals(other.interval) &&
					txHits.keySet().equals(other.txHits.keySet());
		}
	}
}
