package org.stjude.compbio.common;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.stjude.compbio.util.HashMultimap;
import org.stjude.compbio.util.Multimap;

/**
 * Map implementation whose keys are reference names; allows fuzzy key matching
 * to account for reference naming convention variations.
 *
 * @param <V> map value type
 */
public class RefNameMap<V> implements Map<String, V> {
	/**
	 * Smallest acceptable lenght of integral substring to be considered
	 */
	public static final int MIN_INTEGRAL_SUBSTRING_LEN = 4;
	/**
	 * Regex pattern used to find maximal integral substrings
	 */
	private static final Pattern MIS_PATTERN = Pattern.compile("[0-9]+");
	
	/**
	 * Returns the alternate form of a ref name, formed by adding or removing
	 * the chr prefix.
	 * 
	 * @param refName reference name
	 * @return alternate form
	 */
	public static String altForm(String refName) {
		if(refName.startsWith("chr")) {
			return refName.substring(3);
		}
		else {
			return "chr" + refName;
		}
	}

	/**
	 * Does a quick check to see if a string would be considered to match a key
	 * if that key were in a RefNameMap.
	 */
	public static boolean matches(String key, String refName) {
		// Check equality first, for improved performance
		if(key.equals(refName)) return true;
		
		// Build a single-entry map and use containsKey to do the check.  This
		// is inefficient, but ensures that the logic is identical.
		RefNameMap<Boolean> map = new RefNameMap<Boolean>();
		map.put(key, true);
		return map.containsKey(refName);
	}
	
	/**
	 * Underlying map; there are no constraints on this map
	 */
	private Map<String, V> map;
	/**
	 * Mapping from maximal integral substrings to keys
	 */
	private Multimap<String, String, List<String>> misKeys;
	
	/**
	 * Constructs a RefNameMap with an empty HashMap underlying it
	 */
	public RefNameMap() {
		this(new HashMap<String, V>());
	}
	/**
	 * Constructs a RefNameMap from an existing map.
	 * 
	 * The map is copied, so updates to the map will not affect the RefNameMap
	 * and vice versa.
	 * 
	 * @param map initial map
	 */
	public RefNameMap(Map<String, V> map) {
		this.map = copyMap(map);
		this.misKeys = buildMisKeys(map);
	}
	
	/**
	 * Attempts to copy a map using the same class as the original.
	 * 
	 * @param src original map
	 * @return copy
	 */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	private Map<String, V> copyMap(Map<String, V> src) {
		try {
			Constructor<? extends Map> ctor =
					src.getClass().getConstructor(Map.class);
			return ctor.newInstance(src);
		} catch(Exception e) {
			return new HashMap<String, V>(src);
		}
	}
	
	/**
	 * Builds a mapping from maximal integral subkeys to the map keys that have
	 * that key.
	 * 
	 * @param map map from which to read keys
	 * @return multimap from maximal integral subkey to keys
	 */
	private Multimap<String, String, List<String>> buildMisKeys(
			Map<String, V> map) {
		Multimap<String, String, List<String>> out =
			HashMultimap.newListInstance();
		for(String key: map.keySet()) {
			List<String> miss = getMiss(key);
			for(String mis: miss) out.add(mis, key);
		}
		return out;
	}
	
	/**
	 * Gets the maximal integral substrings of a string above the minimum
	 * 
	 * @param key string
	 * @return list of maximal integral substrings
	 */
	private List<String> getMiss(String key) {
		List<String> out = new ArrayList<String>();
		Matcher matcher = MIS_PATTERN.matcher(key);
		int start = 0;
		while(matcher.find(start)) {
			String is = matcher.group();
			if(is.length() >= MIN_INTEGRAL_SUBSTRING_LEN) out.add(is);
			start = matcher.end();
		}
		return out;
	}
	
	/**
	 * Gets the key that should be used for a key lookup; this is where the
	 * fuzzy logic lies.
	 * 
	 * @param key key to look up (will be converted to String using toString()))
	 * @return key to use, or null if no match in the map
	 */
	public String getMappedKey(Object key) {
		return getMappedKey(key.toString());
	}

	/**
	 * Gets the key that should be used for a key lookup; this is where the
	 * fuzzy logic lies.
	 * 
	 * Rules:
	 * 1. If the map contains the key as-is, then use it
	 * 2. Try adding or removing chr prefix
	 * 3. Try switching M/MT
	 * 4. Try longest integer substring match
	 * 
	 * @param key key to look up
	 * @return key to use, or null if no match in the map
	 */
	private String getMappedKeySimple(String key) {
		// Try key itself
		if(map.containsKey(key)) return key;
		
		// Try alternate form (w/ or w/out chr- prefix), and drop chr prefix
		// for following tests
		if(key.startsWith("chr")) {
			key = key.substring(3);
			if(map.containsKey(key)) return key;
		}
		else {
			String altFormKey = "chr" + key;
			if(map.containsKey(altFormKey)) return altFormKey;
		}
		
		// Now, key lacks chr prefix
		
		// Try M/MT
		String altMKey = null;
		if(key.equals("M")) {
			altMKey = "MT";
		}
		else if(key.equals("MT")) {
			altMKey = "M";
		}
		if(altMKey != null) {
			if(map.containsKey(altMKey)) return altMKey;
			altMKey = altForm(altMKey);
			if(map.containsKey(altMKey)) return altMKey;
		}
		
		return null;
	}

	/**
	 * Gets the key that should be used for a key lookup; this is where the
	 * fuzzy logic lies.
	 * 
	 * Rules:
	 * 1. If the map contains the key as-is, then use it
	 * 2. Try adding or removing chr prefix
	 * 3. Try switching M/MT
	 * 4. Try longest integer substring match
	 * 
	 * @param key key to look up
	 * @return key to use, or null if no match in the map
	 */
	private String getMappedKey(String key) {
		String mappedKey = getMappedKeySimple(key);
		if(mappedKey != null) return mappedKey;
		
		// Try X=N+1, Y=N+2
		// Requires numeric key
		int num = 0;
		try {
			num = Integer.parseInt(key);
		} catch(NumberFormatException e) {}
		if(num > 0) {
			String xyKey;
			// Try X (num should be 1 more than max)
			if(getMappedKeySimple(Integer.toString(num - 1)) != null) {
				if((xyKey = getMappedKey("X")) != null) return xyKey;
			}
			// Try Y (num should be 2 more than max)
			else if(getMappedKeySimple(Integer.toString(num - 2)) != null) {
				if((xyKey = getMappedKey("Y")) != null) return xyKey;
			}
		}
		
		// Finally, try longest maximal substring matching
		List<String> miss = getMiss(key);
		Set<String> candidates = new HashSet<String>();
		for(String mis: miss) {
			List<String> misCandidates = misKeys.get(mis);
			if(misCandidates != null) {
				candidates.addAll(misCandidates);
				if(candidates.size() > 1) break;
			}
		}
		if(candidates.size() == 1) return candidates.iterator().next();
		
		// Nothing worked; return null
		return null;
	}
	
	/**
	 * Adds a mapping from a ref name to itself; useful for the special case
	 * where the map is used to canonicalize names.
	 * 
	 * @param refName the ref name to add
	 */
	public void add(V refName) {
		put(refName.toString(), refName);
	}

	public int size() {
		return map.size();
	}

	public boolean isEmpty() {
		return map.isEmpty();
	}

	public boolean containsKey(Object key) {
		return getMappedKey(key) != null;
	}

	public boolean containsValue(Object value) {
		return map.containsValue(value);
	}

	public V get(Object key) {
		key = getMappedKey(key);
		if(key == null) return null;
		return map.get(key);
	}

	public V put(String key, V value) {
		if(!map.containsKey(key)) {
			for(String mis: getMiss(key)) misKeys.add(mis, key);
		}
		return map.put(key, value);
	}

	public V remove(Object key) {
		String mappedKey = getMappedKey(key);
		if(mappedKey == null) {
			return map.remove(key);
		}
		else {
			for(String mis: getMiss(mappedKey)) misKeys.removeSimple(mis, mappedKey);
			return map.remove(mappedKey);
		}
	}

	public void putAll(Map<? extends String, ? extends V> m) {
		map.putAll(m);
	}

	public void clear() {
		map.clear();
	}

	/**
	 * Returns the set of actual keys which does not include all mapped values.
	 * 
	 * Key matching is fuzzy in the map, but the set is just a regular set, so
	 * there may be values k such that keySet().contains(k) == false, but
	 * containsKey(k) == true.
	 */
	public Set<String> keySet() {
		return map.keySet();
	}

	public Collection<V> values() {
		return map.values();
	}

	public Set<java.util.Map.Entry<String, V>> entrySet() {
		return map.entrySet();
	}

}
