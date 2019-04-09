package org.stjude.compbio.util.io;

import java.util.Collection;
import java.util.Collections;
import java.util.SortedMap;

/**
 * A record in an external file/resource that is being read/written.
 * 
 * This structure holds values for a single record and provides access by
 * IoField.
 * 
 * All values are Strings
 * 
 * @see IoField
 * 
 * @param <K> the key type
 */
public class IoRecord<K> {
	/**
	 * The IoFormat, for mapping names, keys, and ordinals to fields
	 */
	private final IoFormat<K> format;
	/**
	 * The data
	 */
	private SortedMap<IoField<K>, String> data;

	public IoRecord(IoFormat<K> format, SortedMap<IoField<K>, String> data) {
		this.format = format;
		this.data = data;
	}
	
	IoFormat<K> getFormat() {
		return format;
	}
	
	/**
	 * Unmodifiable map view of the underlying data
	 */
	public SortedMap<IoField<K>, String> getData() {
		return Collections.unmodifiableSortedMap(data);
	}

	/**
	 * Retrieves a value using the field
	 */
	public String valueForField(IoField<K> field) {
		if(field == null) return null;
		return data.get(field);
	}
	
	/**
	 * Retrieves a value using the name
	 */
	public String valueForName(String name) {
		return valueForField(format.getNameMap().get(name));
	}
	
	/**
	 * Retrieves a value using the key
	 */
	public String valueForKey(K key) {
		return valueForField(format.getKeyMap().get(key));
	}
	
	/**
	 * Retrieves a value using the ordinal
	 */
	public String valueForOrdinal(Integer ordinal) {
		return valueForField(format.getOrdinalMap().get(ordinal));
	}
	
	@Override
	public String toString() {
		return data.toString();
	}

	/**
	 * Returns an unmodifiable Collection view of the values in the record.
	 * 
	 * The returned Collection's iterator will iterate over the values in the
	 * order of the record's fields.
	 */
	public Collection<String> values() {
		return Collections.unmodifiableCollection(data.values());
	}
}
