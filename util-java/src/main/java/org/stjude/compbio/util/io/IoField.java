package org.stjude.compbio.util.io;

import org.stjude.compbio.util.DelimitedListBuilder;

/**
 * A field in an external file/resource that is being read/written.
 * 
 * This would correspond to a column in a tabular text file for example.  It
 * does not actually contain data; it is descriptive.
 * 
 * The name is what the field is called outside of the application, and the
 * key is an application-specific identification.  This allows you to have
 * multiple file formats that contain the same data with different names/
 * headings.
 * 
 * Ordinal typically corresponds to column number.
 * 
 * @param <K> the key type
 */
public class IoField<K> implements Comparable<IoField<K>> {
	/**
	 * The external name/header (may be null)
	 */
	private String name;
	/**
	 * The internal key (may be null)
	 */
	private K key;
	/**
	 * Ordinal, which may or may not be used
	 */
	private Integer ordinal;
	
	public IoField(String name, Integer ordinal) {
		this(name, ordinal, null);
	}
	
	public IoField(String name, Integer ordinal, K key) {
		this.name = name;
		this.ordinal = ordinal;
		this.key = key;
	}
	
	/**
	 * Copy constructor
	 */
	public IoField(IoField<K> src) {
		this(src.name, src.ordinal);
		this.key = src.key;
	}
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public K getKey() {
		return key;
	}
	public void setKey(K key) {
		this.key = key;
	}
	public Integer getOrdinal() {
		return ordinal;
	}
	public void setOrdinal(Integer ordinal) {
		this.ordinal = ordinal;
	}
	
	@Override
	public String toString() {
		DelimitedListBuilder dlb = new DelimitedListBuilder("/");
		dlb.append(ordinal == null ? "?" : ordinal);
		dlb.append(name == null ? "<noname>" : name);
		dlb.append(key == null ? "<nokey>" : key);
		return dlb.toString();
	}
	
	public int compareTo(IoField<K> o) {
		return ordinalCompare(this, o);
	}

	/**
	 * Performs ordinal-based comparison of two fields
	 * 
	 * This is designed to be a part of a custom Comparator.
	 */
	public static <K> int ordinalCompare(IoField<K> o1, IoField<K> o2) {
		return (o1.ordinal == null ? Integer.MAX_VALUE : o1.ordinal) -
				(o2.ordinal == null ? Integer.MAX_VALUE : o2.ordinal);
	}
	
	/**
	 * Performs key-based comparison of two fields, when the keys are
	 * comparable.
	 * 
	 * This is designed to be a part of a custom Comparator.
	 */
	public static <K extends Comparable<K>> int keyCompare(
			IoField<K> o1, IoField<K> o2) {
		return o1.getKey().compareTo(o2.getKey());
	}
}
