package org.stjude.compbio.common.seq;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A translation table.
 * 
 * Each entry contains an amino acid short and long name, and a flag indicating
 * whether or not it may be initiating.
 * 
 * Keys are uppercase codons with no ambiguity codes.
 */
public class TranslTable {
	/**
	 * Name of the standard transl table
	 */
	private static final String STANDARD_TRANSL_TABLE_NAME = "Standard";
	/**
	 * Pattern for parsing tables in the NCBI "t" format
	 */
	private static final Pattern NCBI_T_ENTRY = Pattern.compile(
			"([AaCcGgTt]{3}) ([A-Za-z*]) ([A-Za-z]{3})( i)?");
	
	/**
	 * Cache of named transl tables
	 */
	private static final Map<String, TranslTable> namedCache =
			new HashMap<String, TranslTable>();

	public static class Entry {
		/**
		 * Three-character unambiguous uppercase codon
		 */
		public final String codon;
		/**
		 * Single-character amino acid code (* for terminating)
		 */
		public final char aa;
		/**
		 * Three-character amino acid abbreviation (Ter for terminating)
		 */
		public final String aaAbbrev;
		/**
		 * Whether or not it can be initiating
		 */
		public final boolean initiating;
		
		private Entry(String codon, char aa, String aaAbbrev, boolean initiating) {
			super();
			this.codon = codon;
			this.aa = aa;
			this.aaAbbrev = aaAbbrev;
			this.initiating = initiating;
		}
	}
	
	/**
	 * The map defining the table
	 */
	private Map<String, Entry> map;
	
	/**
	 * Constructs a TranslTable with an empty map 
	 */
	private TranslTable() {
		this.map = new HashMap<String, Entry>();
	}

	/**
	 * Returns the standard TranslTable (i.e. the one named "Standard")
	 * @throws IOException 
	 */
	public static TranslTable getStandardTranslTable() throws IOException {
		return getNamedTranslTable(STANDARD_TRANSL_TABLE_NAME);
	}
	
	/**
	 * Gets an entry for a codon
	 */
	public Entry getEntry(String codon) {
		return map.get(codon);
	}
	
	/**
	 * Translates a codon
	 */
	public char translate(String codon) {
		return map.get(codon).aa;
	}
	
	/**
	 * Reads a TranslTable from an available bundled resource file using
	 * the name.
	 * 
	 * The resource file is named TranslTableNAME.txt, where NAME is the
	 * specified name.
	 * @throws IOException 
	 */
	public static TranslTable getNamedTranslTable(String name)
	throws IOException {
		// Check cache
		if(namedCache.containsKey(name)) return namedCache.get(name);
		
		// Read it
		BufferedReader reader = new BufferedReader(new InputStreamReader(
				TranslTable.class.getResourceAsStream(getResourceName(name))));
		TranslTable tt = readTranslTable(reader);
		
		// Put it in the cache before returning it
		namedCache.put(name, tt);
		return tt;
	}

	/**
	 * Builds a TranslTable by reading from an input stream.
	 * 
	 * The expected format is the one used by NCBI that includes the three-
	 * letter amino acid abbreviations.  See 
	 * http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1
	 * 
	 * @param reader BufferedReader open on input data
	 * @return new TranslTable
	 * @throws IOException
	 */
	public static TranslTable readTranslTable(BufferedReader reader)
	throws IOException {
		// Construct empty table
		TranslTable tt = new TranslTable();
		
		// Loop through the entries in the table
		String line;
		while((line = reader.readLine()) != null) {
			Matcher matcher = NCBI_T_ENTRY.matcher(line);
			int start = 0;
			// For each one, parse and add a new entry to the map
			while(matcher.find(start)) {
				String codon = matcher.group(1).toUpperCase();
				tt.map.put(codon, new Entry(
						codon,
						matcher.group(2).charAt(0),
						matcher.group(3),
						" i".equals(matcher.group(4))
						));
				start = matcher.end();
			}
		}
		
		// Return the new table
		return tt;
	}

	/**
	 * Gets the resource name for a transl table name
	 * @param name transl table name
	 * @return resource name without package
	 */
	private static String getResourceName(String name) {
		return "TranslTable" + name + ".txt";
	}
	
	
}
