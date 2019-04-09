package org.stjude.compbio.util;

import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

/**
 * Map implementation whose keys are regex patterns and for which hits are based
 * on key regex matches, not key object equality.
 * 
 * When searching, search objects are converted to string and then matched
 * against the patterns.  Note that even Patterns are searched this way, so you
 * must use containsPattern() and getUsingPattern() instead.  Since a string may 
 * hit more than one pattern, the first pattern put into the map is the one that
 * is used.  
 * 
 * Note that this is a significant divergence from the way maps usually work.
 * 
 * Lookup time is linear with the size of the map.  So, it is slow for large
 * numbers of entries.
 * 
 * Null values are not permitted.
 *
 * @param <V> value type
 */
public class RegexMap<V> implements Map<Pattern, V> {
	/**
	 * The underlying map
	 */
	LinkedHashMap<Pattern, V> map;
	
	/**
	 * Last searched-for object
	 */
	Object lastSearchedObject = null;
	/**
	 * Last returned value
	 */
	V lastReturnedValue = null;
	
	public RegexMap() {
		map = new LinkedHashMap<Pattern, V>();
	}
	
	public int size() {
		return map.size();
	}

	public boolean isEmpty() {
		return map.isEmpty();
	}

	public boolean containsKey(Object key) {
		return get(key) != null;
	}

	public boolean containsValue(Object value) {
		if(value == null) return false;
		for(V v: map.values()) if(v.equals(value)) return true;
		return false;
	}

	public V get(Object key) {
		for(Pattern pattern: map.keySet()) {
			if(pattern.matcher(key.toString()).find()) return map.get(pattern);
		}
		return null;
	}

	public V put(Pattern key, V value) {
		return map.put(key, value);
	}

	public V remove(Object key) {
		return map.remove(key);
	}

	public void putAll(Map<? extends Pattern, ? extends V> m) {
		map.putAll(m);
	}

	public void clear() {
		map.clear();
	}

	public Set<Pattern> keySet() {
		return map.keySet();
	}

	public Collection<V> values() {
		return map.values();
	}

	public Set<java.util.Map.Entry<Pattern, V>> entrySet() {
		return map.entrySet();
	}

}
