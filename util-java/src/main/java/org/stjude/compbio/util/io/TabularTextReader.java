package org.stjude.compbio.util.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.regex.Pattern;

/**
 * A reader for tabular text that returns IoRecords.
 * 
 * Note that this extends LineReader so you can also peek and read raw lines.
 * 
 * You specify in the construction whether or not there should be a header.  If
 * there is a header, then it is read during construction and is immediately
 * available from getHeader().
 * 
 * @param <K> type used as an alternate key for referring to fields
 */
public class TabularTextReader<K> extends AbstractLineReader
implements Iterable<IoRecord<K>> {
	/**
	 * Header options
	 */
	public static class Header {
		/**
		 * Whether the header is allowed
		 */
		private boolean allowed;
		/**
		 * Whether the header is required (cannot be true if allowed is false)
		 */
		private boolean required;
		/**
		 * Regex pattern that identifies a valid header; null to indicate that
		 * all are valid
		 */
		private Pattern pattern;
		/**
		 * Number of lines that are assumed to be ignored header
		 */
		private int numSkipped;
		
		/**
		 * Constructor using all fields
		 */
		private Header(boolean allowed, boolean required, Pattern pattern,
				int numSkipped) {
			super();
			this.allowed = allowed;
			this.required = required;
			this.pattern = pattern;
			this.numSkipped = numSkipped;
		}
		/**
		 * Constructor that compiles regex
		 */
		private Header(boolean allowed, boolean required, String regex,
				int numSkipped) {
			this(allowed, required,
					regex == null ? null : Pattern.compile(regex), numSkipped);
		}
		
		/**
		 * Whether the given header line is valid
		 */
		public boolean isValid(String header) {
			if(header == null) return false;
			return pattern == null || pattern.matcher(header).find();
		}
		
		/**
		 * Returns a header options object for a required header
		 * 
		 * @param regex regex that matches (via find()) a valid header, or null
		 * 			if any header is valid
		 */
		public static Header required(String regex) {
			return new Header(true, true, regex, 0);
		}
		
		/**
		 * Returns a header options object for an optional header.
		 * 
		 * You may not specify a null regex.
		 * 
		 * @param regex regex that matches (via find()) a valid header
		 */
		public static Header optional(String regex) {
			if(regex == null) {
				throw new IllegalArgumentException(
						"No regex given for optional header.");
			}
			return new Header(true, false, regex, 0);
		}
		
		/**
		 * Returns a header options object for a forbidden header.
		 */
		public static Header forbidden() {
			return new Header(false, false, (Pattern)null, 0);
		}
		
		/**
		 * Returns a header options object for an ignored header of n lines.
		 */
		public static Header skipped(int numLines) {
			return new Header(false, false, (Pattern)null, numLines);
		}
	}
	
	/**
	 * Field delimiter (regex)
	 */
	private String delim;
	/**
	 * IoFormat
	 */
	private IoFormat<K> ioFormat;
	/**
	 * Header line
	 */
	private String header;
	
	/**
	 * Opens a keyless TabularTextReader on a File (a common use case)
	 * @throws IOException 
	 */
	public static TabularTextReader<Object> open(
			File file, String delim, Header headerOptions) throws IOException {
		return new TabularTextReader<Object>(file, delim, headerOptions);
	}
	
	public TabularTextReader(BufferedReader reader,
			String delim, Header headerOptions) throws IOException {
		super(reader);
		init(delim, headerOptions);
	}
	public TabularTextReader(File file,
			String delim, Header headerOptions) throws IOException {
		super(file);
		init(delim, headerOptions);
	}
	public TabularTextReader(InputStream stream,
			String delim, Header headerOptions) throws IOException {
		super(stream);
		init(delim, headerOptions);
	}
	
	/**
	 * Performs initialization, which includes header processing
	 * 
	 * @param delim delimiter
	 * @param headerOptions header options
	 * @param comparator optional (may be null) custom field comparator
	 * @throws IOException on any read problem
	 * @throws FormatException if the header doesn't pass specs
	 */
	private void init(String delim, Header headerOptions)
	throws IOException, FormatException {
		this.delim = delim;
		this.ioFormat = new IoFormat<K>();
		this.header = null;
		
		// If there are lines to skip, then do so
		for(int i = 0; i < headerOptions.numSkipped; i++) readLine();
		
		// If there is no header allowed, then we are ready for data
		if(!headerOptions.allowed) return;
		
		// Check validity of the header, and set up the IoFormat if there is a
		// valid header available
		boolean isHeaderValid = headerOptions.isValid(peekLine());
		if(isHeaderValid) {
			header = readLine();
			addFieldsFromHeader(header);
		}
		else if(headerOptions.required) {
			throw new FormatException("Header is not valid.");
		}
	}
	/**
	 * Uses a header line to add fields
	 * @param header header line
	 */
	private void addFieldsFromHeader(String header) {
		String[] parts = header.split(delim);
		for(int i = 0; i < parts.length; i++) {
			ioFormat.addField(new IoField<K>(parts[i], i));
		}
	}
	
	/**
	 * Gets the delimiter
	 */
	public String getDelim() {
		return delim;
	}
	
	/**
	 * Gets the IoFormat
	 */
	public IoFormat<K> getIoFormat() {
		return ioFormat;
	}
	
	/**
	 * Gets the header line, or null if there was none
	 */
	public String getHeader() {
		return header;
	}
	
	/**
	 * Assigns keys to fields based on name.
	 * 
	 * For entries in the map whose names are not present in the format, it may
	 * optionally add a new field.  Whether or not it does depends on the
	 * addNew parameter.  If you specify addNew = false, these will be ignored.
	 * 
	 * @param map map from key to field name
	 * @param addNew true to add fields for new names
	 */
	public void addKeys(Map<K, String> map, boolean addNew) {
		for(Entry<K, String> entry: map.entrySet()) {
			String name = entry.getValue();
			// Update if possible
			if(ioFormat.getNameMap().containsKey(name)) {
				ioFormat.addKey(
						ioFormat.getNameMap().get(name), entry.getKey());
			}
			// Otherwise add if addNew
			else if(addNew) {
				IoField<K> field = new IoField<K>(name, null);
				field.setKey(entry.getKey());
				ioFormat.addField(field);
			}
		}
	}
	
	/**
	 * Assigns keys to fields based on position.
	 * 
	 * For entries in the list beyond the last field in the format, it may
	 * optionally add a new field.  Whether or not it does depends on the
	 * addNew parameter.  If you specify addNew = false, these will be ignored.
	 * 
	 * @param keys list of keys in order (some elements may be null)
	 * @param addNew true to add fields for new names
	 */
	public void addKeys(List<K> keys, boolean addNew) {
		int i = 0;
		for(K key: keys) {
			// Update if possible
			if(ioFormat.getOrdinalMap().containsKey(i)) {
				ioFormat.addKey(ioFormat.getOrdinalMap().get(i), key);
			}
			// Otherwise add
			else if(addNew) {
				IoField<K> field = new IoField<K>(null, i);
				field.setKey(key);
				ioFormat.addField(field);
			}
			// Increment to next field
			++i;
		}
	}
	
	/**
	 * Reads a record.
	 * 
	 * @return IoRecord for next line or null if input is exhausted
	 * @throws IOException 
	 */
	public IoRecord<K> readIoRecord() throws IOException {
		// Read line and bail if input is exhausted
		String line = readLine();
		if(line == null) return null;
		
		// Parse into record data
		SortedMap<IoField<K>, String> data =
				new TreeMap<IoField<K>, String>();
		String[] parts = line.split(delim);
		for(int i = 0; i < parts.length; i++) {
			IoField<K> field = getOrAddField(i);
			data.put(field, parts[i]);
		}
		
		// Return record
		return new IoRecord<K>(ioFormat, data);
	}
	
	/**
	 * Gets the field for an ordinal, creating one if none exists
	 * @param ordinal the ordinal
	 * @return field
	 */
	private IoField<K> getOrAddField(int ordinal) {
		if(ioFormat.getOrdinalMap().containsKey(ordinal)) {
			return ioFormat.getOrdinalMap().get(ordinal);
		}
		else {
			IoField<K> field = new IoField<K>(null, ordinal);
			ioFormat.addField(field);
			return field;
		}
	}
	
	public Iterator<IoRecord<K>> iterator() {
		return new AbstractIterator<IoRecord<K>>() {
			public IoRecord<K> next() {
				try {
					return readIoRecord();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
			}
		};
	}

}
