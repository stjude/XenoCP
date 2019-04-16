package org.stjude.compbio.util.io;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * Basic description of the format of an external file/resource that is being
 * read/written.
 * 
 * All IoFields are assigned an ordinal when they come in if they do not already
 * have one, so all fields coming out of an IoFormat will have a non-null
 * ordinal.
 * 
 * @see IoField
 *
 * @param <K> the key type
 */
public class IoFormat<K> {
	/**
	 * All fields, sorted by ordinal
	 */
	private SortedSet<IoField<K>> fields;
	/**
	 * Map from names to fields
	 */
	private Map<String, IoField<K>> nameMap;
	/**
	 * Map from keys to fields
	 */
	private Map<K, IoField<K>> keyMap;
	/**
	 * Map from ordinals to fields
	 */
	private SortedMap<Integer, IoField<K>> ordinalMap;
	
	/**
	 * Constructs an empty IoFormat
	 */
	public IoFormat() {
		this.fields = new TreeSet<IoField<K>>(new Comparator<IoField<K>>() {
			public int compare(IoField<K> o1, IoField<K> o2) {
				return IoField.ordinalCompare(o1, o2);
			}});
		this.nameMap = newMap();
		this.keyMap = newMap();
		this.ordinalMap = new TreeMap<Integer, IoField<K>>();
	}

	/**
	 * Generic method for creating a new HashMap to IoField<K>
	 */
	private <T> Map<T, IoField<K>> newMap() {
		return new HashMap<T, IoField<K>>();
	}

	/**
	 * Gets the full list of fields
	 */
	public SortedSet<IoField<K>> getFields() {
		return Collections.unmodifiableSortedSet(fields);
	}

	/**
	 * Gets map from names to fields
	 */
	public Map<String, IoField<K>> getNameMap() {
		return Collections.unmodifiableMap(nameMap);
	}

	/**
	 * Gets map from keys to fields
	 */
	public Map<K, IoField<K>> getKeyMap() {
		return Collections.unmodifiableMap(keyMap);
	}

	/**
	 * Gets map from ordinals to fields
	 */
	public Map<Integer, IoField<K>> getOrdinalMap() {
		return Collections.unmodifiableMap(ordinalMap);
	}
	
	/**
	 * Adds a field.
	 * 
	 * If there is no ordinal on the field, then it is assigned automatically.
	 * This updates the field object that is passed in.
	 */
	public void addField(IoField<K> field) {
		// Set ordinal if it is not already
		if(field.getOrdinal() == null) {
			field.setOrdinal(getNextOrdinal());
		}
		// Add the field
		fields.add(field);
		addFieldToMap(nameMap, field, field.getName());
		addFieldToMap(keyMap, field, field.getKey());
		addFieldToMap(ordinalMap, field, field.getOrdinal());
	}
	
	/**
	 * Adds a field to the end using name only.
	 */
	public void addFieldByName(String name) {
		addField(new IoField<K>(name, null));
	}
	
	/**
	 * Adds a key setting to an existing key-less field.
	 * 
	 * Use this method rather than setting it directly, as this allows the
	 * indexes to be updated.
	 * 
	 * @param field the field to update
	 * @param key the key to add
	 */
	public void addKey(IoField<K> field, K key) {
		field.setKey(key);
		keyMap.put(key, field);
	}

	/**
	 * Gets the next ordinal to assign
	 */
	private int getNextOrdinal() {
		if(fields.isEmpty()) {
			return 0;
		}
		else {
			return fields.last().getOrdinal() + 1;
		}
	}

	/**
	 * Adds an entry for a field to a map if the appropriate map key is
	 * present.
	 * 
	 * @param map map to add to
	 * @param field the field to add
	 * @param mapKey the map key value
	 */
	private <T> void addFieldToMap(
			Map<T, IoField<K>> map, IoField<K> field, T mapKey) {
		if(mapKey == null) return;
		map.put(mapKey, field);
	}
	
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for(IoField<K> field: fields) {
			sb.append(field.getOrdinal() + ". " + field + "\n");
		}
		return sb.toString();
	}
}
