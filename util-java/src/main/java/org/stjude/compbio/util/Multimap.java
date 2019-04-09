package org.stjude.compbio.util;

import java.util.Collection;
import java.util.Map;

/**
 * A Map whose values are collections with methods for accessing the members of
 * the collections.
 * 
 * You in effect specify whether this maps to lists or sets, etc.
 * 
 * @param <K> key type
 * @param <V> collections' value type
 * @param <C> type of the collection 
 */
public interface Multimap<K, V, C extends Collection<V>> extends Map<K, C> {
	/**
	 * Adds a value to the collection with the given key.
	 * 
	 * @param key map key
	 * @param value value to add to collection
	 * @return true if added, false otherwise
	 */
	boolean add(K key, V value);
	
	/**
	 * Removes a value from the collection with the given key.
	 * 
	 * @param key map key
	 * @param value value to remove from collection
	 * @return true if removed, false otherwise
	 */
	boolean removeSimple(K key, V value);
	
	/**
	 * Removes a value from the collection with the given key, and removes the
	 * collection if it is then empty.
	 * 
	 * If the removal fails but the collection was empty to begin with, then the
	 * collection is NOT removed.
	 * 
	 * @param key map key
	 * @param value value to remove from collection
	 * @return true if removed, false otherwise
	 */
	boolean removeCascade(K key, V value);
	
	/**
	 * Returns true if a key's collection contains the value.
	 * 
	 * @param key map key
	 * @param value value to test
	 * @return true if the collection contains the value
	 */
	boolean contains(K key, V value);
	
	/**
	 * Re-implementation of get that creates and returns an empty collection
	 * instead of null
	 * 
	 * @param key map key
	 * @return collection for that key, lazily initialized to empty collection
	 */
	C getLazyInit(K key);
	
	/**
	 * Re-implementation of remove that creates and returns an empty collection
	 * instead of null
	 * 
	 * @param key map key
	 * @return collection for that key, lazily initialized to empty collection
	 */
	C removeLazyInit(K key);

	/**
	 * Adds all mappings in another Multimap into this one.
	 * 
	 * This is the equivalent of calling add(K, V) for every key-value pair
	 * expressed by the other map
	 * 
	 * @param other the map from which you are adding values
	 */
	void addAll(Multimap<K, V, C> other);
}
