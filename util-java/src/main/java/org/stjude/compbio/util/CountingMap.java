package org.stjude.compbio.util;

import java.util.LinkedHashMap;
import java.util.Map;

/**
 * HashMap implementation that is used for counting objects.
 * 
 * Use the increment method to increment the count for a key.  It is safe to
 * use even if the key has not been added to the map yet.
 * 
 * The getCount method returns 0 instead of null if there is no entry.
 *
 * @param <T> type of object being counted
 */
public class CountingMap<T> extends LinkedHashMap<T, Integer> {
	private static final long serialVersionUID = 1L;

	/**
	 * Increments the count for object key
	 * @param key the key identifying which count to increment
	 * @return new value
	 */
	public int increment(T key) {
		return increment(key, 1);
	}

	/**
	 * Increments the count for object key
	 * @param key the key identifying which count to increment
	 * @param amount amount to increment
	 * @return new value
	 */
	public int increment(T key, int amount) {
		int val = getCount(key) + amount;
		put(key, val);
		return val;
	}

	/**
	 * Increments and adds counts from another Map
	 * @param added the map whose counts to add
	 */
	public void addAll(Map<T, Integer> addend) {
		for(Map.Entry<T, Integer> entry: addend.entrySet()) {
			increment(entry.getKey(), entry.getValue());
		}
	}
	
	/**
	 * Gets the count for a key, which will return 0 if the key is not in the 
	 * map.
	 * 
	 * @param key key for which you are getting count
	 * @return count associated with that key or 0 if there is none
	 */
	public int getCount(Object key) {
		return containsKey(key) ? get(key) : 0;
	}

	/**
	 * Sets the value to val for all keys; this is often used to intialize
	 * several counts to 0.
	 * 
	 * @param keys keys to set to the value
	 * @param val value to set (often 0)
	 */
	public void setAll(Iterable<T> keys, int val) {
		for(T key: keys) put(key, val);
	}

	/**
	 * Sets the value to val for all keys; this is often used to intialize
	 * several counts to 0.
	 * 
	 * @param keys keys to set to the value
	 * @param val value to set (often 0)
	 */
	public void setAll(T[] keys, int val) {
		for(T key: keys) put(key, val);
	}
	
	/**
	 * Create a counting map with a set of initialized counters.
	 * 
	 * This is useful in cases where some counters may never increment, but you
	 * want to use keySet to get the full set.
	 * 
	 * @param keys keys to initialize
	 */
	public static <T> CountingMap<T> newInitialized(int val, T... keys) {
		CountingMap<T> out = new CountingMap<T>();
		out.setAll(keys, val);
		return out;
	}
	
	/**
	 * Formats an entry in the map for printing
	 */
	public static String formatEntry(Map.Entry<?, Integer> entry) {
		return String.format(
				"%10d %s",
				entry.getValue(),
				entry.getKey().toString());
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		for(Map.Entry<T, Integer> entry: entrySet()) {
			sb.append(formatEntry(entry) + "\n");
		}
		return sb.toString();
	}
}
