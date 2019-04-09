package org.stjude.compbio.util;

import java.util.Arrays;

/**
 * Simple convenience class for creating comma-separated list strings
 */
public class DelimitedListBuilder {
	/**
	 * Delimiter
	 */
	private final String delim;
	/**
	 * StringBuilder used to build the string
	 */
	private final StringBuilder builder;
	/**
	 * Whether or not the list is started (if true, then a delimiter will
	 * precede next item)
	 */
	private boolean listStarted;
	/**
	 * Number of fields
	 */
	private int numFields;
	/**
	 * String value of null (default "")
	 */
	private String nullStr = "";

	/**
	 * Constructs list builder with empty list and "," delimiter.
	 */
	public DelimitedListBuilder() {
		this(",");
	}

	/**
	 * Constructs list builder with empty list and specified delimiter.
	 * @param delim delimiter to use
	 */
	public DelimitedListBuilder(String delim) {
		this(delim, new StringBuilder(), false);
	}

	/**
	 * Constructs list builder with specified string builder and delimiter.
	 * @param delim delimiter to use
	 * @param builder StringBuilder to use for building the list
	 * @param listStarted true if there are list items in builder already
	 */
	public DelimitedListBuilder(String delim, StringBuilder builder,
			boolean listStarted) {
		this.delim = delim;
		this.builder = builder;
		this.listStarted = listStarted;
	}

	/**
	 * Constructs list builder with specified list and delimiter
	 * @param delim delimiter to use
	 * @param init partial list as a starting point
	 * @param listStarted true if there are list items in builder already
	 */
	public DelimitedListBuilder(String delim, String init,
			boolean listStarted) {
		this(delim, new StringBuilder(init), listStarted);
	}

	/**
	 * Gets the string representation of null to use.
	 * 
	 * When you pass null to one of the append methods, this is the string that
	 * is actually added
	 */
	public String getNullStr() {
		return nullStr;
	}

	/**
	 * Sets the string representation of null to use.
	 * 
	 * When you pass null to one of the append methods, this is the string that
	 * is actually added (default is "")
	 */
	public void setNullStr(String nullStr) {
		this.nullStr = nullStr;
	}
	
	public void append(String item) {
		if(item == null) item = nullStr;
		if(listStarted) builder.append(delim); else listStarted = true;
		builder.append(item);
		++numFields;
	}

	public void append(Object item) {
		append(item == null ? null : item.toString());
	}

	public void appendAll(Iterable<?> items) {
		for(Object item: items) append(item);
	}

	/**
	 * Appends item until there are a total of numFields fields
	 * 
	 * @param numFields target number of fields
	 * @param item item to append
	 */
	public void padTo(int numFields, Object item) {
		while(this.numFields < numFields) append(item);
	}
	
	public void clear() {
		builder.delete(0, builder.length());
		listStarted = false;
		numFields = 0;
	}

	public static String implode(String delim, Iterable<?> items) {
		DelimitedListBuilder dlb = new DelimitedListBuilder(delim);
		dlb.appendAll(items);
		return dlb.toString();
	}

	public static String implode(String delim, Object[] items) {
		return implode(delim, Arrays.asList(items));
	}
	
	/**
	 * Returns the string representation of the list
	 */
	@Override
	public String toString() {
		return builder.toString();
	}
}
