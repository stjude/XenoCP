package org.stjude.compbio.util.io;

import java.util.Collection;
import java.util.Collections;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * Builder for IoRecords
 * 
 * @see IoField
 * 
 * @param <K> the key type
 */
public class IoRecordBuilder<K> {
	/**
	 * The IoFormat, for mapping names, keys, and ordinals to fields
	 */
	private final IoFormat<K> format;
	/**
	 * The data
	 */
	private SortedMap<IoField<K>, String> data;

	public IoRecordBuilder(IoFormat<K> format) {
		this.format = format;
		this.data = new TreeMap<IoField<K>, String>();
	}
	
	public IoRecordBuilder(IoRecord<K> source) {
		this(source.getFormat());
		init(source);
	}
	
	/**
	 * Removes all data from the record
	 */
	public void clear() {
		data.clear();
	}
	
	/**
	 * Removes all data and copies in values from a source record.
	 */
	public void init(IoRecord<K> source) {
		data.clear();
		data.putAll(source.getData());
	}
	
	/**
	 * Sets a value using the field
	 */
	public void setValueForField(IoField<K> field, String value) {
		if(field == null) return;
		if(!format.getFields().contains(field)) {
			throw new IllegalArgumentException(
					"Field " + field + " is not in format.");
		}
		data.put(field, value);
	}
	
	/**
	 * Sets a value using the name
	 */
	public void setValueForName(String name, String value) {
		setValueForField(format.getNameMap().get(name), value);
	}
	
	/**
	 * Sets a value using the key
	 */
	public void setValueForKey(K key, String value) {
		setValueForField(format.getKeyMap().get(key), value);
	}
	
	/**
	 * Sets a value using the ordinal
	 */
	public void setValueForOrdinal(Integer ordinal, String value) {
		setValueForField(format.getOrdinalMap().get(ordinal), value);
	}

	
	public IoRecord<K> getRecord() {
		return new IoRecord<K>(format, data);
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
