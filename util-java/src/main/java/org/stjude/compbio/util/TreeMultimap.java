package org.stjude.compbio.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * TreeMap-based Multimap implementation
 *
 * @param <K> key type
 * @param <V> value type
 * @param <C> type of collection that contains the values
 */
public class TreeMultimap<K, V, C extends Collection<V>> extends TreeMap<K, C>
implements SortedMultimap<K, V, C> {
	private static final long serialVersionUID = 1L;

	/**
	 * Makes a new empty collection
	 * 
	 * @param <T> collection type
	 */
	public static interface Factory<T> {
		T newInstance();
	}
	
	/**
	 * The factory to make empty collections
	 */
	private final Factory<C> factory;

	private TreeMultimap(Factory<C> factory) {
		super();
		this.factory = factory;
	}

	/**
	 * Copy constructor
	 * 
	 * For non-copy construction use static factory methods
	 */
	private TreeMultimap(TreeMultimap<K, V, C> src) {
		super();
		this.factory = src.factory;
		this.addAll(src);
	}
	
	public boolean add(K key, V value) {
		return getLazyInit(key).add(value);
	}

	public boolean removeSimple(K key, V value) {
		return removeOptionalCascade(key, value, false);
	}

	public boolean removeCascade(K key, V value) {
		return removeOptionalCascade(key, value, true);
	}

	private boolean removeOptionalCascade(K key, V value, boolean cascade) {
		C coll = get(key);
		if(coll == null) return false;
		if(coll.remove(value)) {
			if(cascade && coll.isEmpty()) remove(key);
			return true;
		}
		else {
			return false;
		}
	}

	public boolean contains(K key, V value) {
		C coll = get(key);
		if(coll == null) return false;
		return coll.contains(value);
	}

	public C getLazyInit(K key) {
		C coll = get(key);
		if(coll == null) {
			put(key, coll = factory.newInstance());
		}
		return coll;
	}

	public C removeLazyInit(K key) {
		C coll = remove(key);
		if(coll == null) {
			coll = factory.newInstance();
		}
		return coll;
	}

	public void addAll(Multimap<K, V, C> other) {
		for(K key: other.keySet()) {
			getLazyInit(key).addAll(other.get(key));
		}
	}

	/**
	 * Constructs an empty HashMultimap using HashSet collection.
	 */
	public static <K, V> TreeMultimap<K, V, Set<V>> newHashSetInstance() {
		return new TreeMultimap<K, V, Set<V>>(new Factory<Set<V>>() {
			public Set<V> newInstance() {
				return new HashSet<V>();
			}
		});
	}

	/**
	 * Constructs an empty HashMultimap using TreeSet collection.
	 */
	public static <K, V> TreeMultimap<K, V, TreeSet<V>>
			newTreeSetInstance() {
		return new TreeMultimap<K, V, TreeSet<V>>(
				new Factory<TreeSet<V>>() {
					public TreeSet<V> newInstance() {
						return new TreeSet<V>();
					}
		});
	}

	/**
	 * Constructs an empty HashMultimap using a List collection.
	 */
	public static <K, V> TreeMultimap<K, V, List<V>> newListInstance() {
		return new TreeMultimap<K, V, List<V>>(
				new Factory<List<V>>() {
					public List<V> newInstance() {
						return new ArrayList<V>();
					}
				});
	}

	/**
	 * Constructs an empty HashMultimap using ArrayList collection.
	 */
	public static <K, V> TreeMultimap<K, V, ArrayList<V>>
	newArrayListInstance() {
		return new TreeMultimap<K, V, ArrayList<V>>(
				new Factory<ArrayList<V>>() {
					public ArrayList<V> newInstance() {
						return new ArrayList<V>();
					}
				});
	}

	/**
	 * Constructs an empty HashMultimap using LinkedList collection.
	 */
	public static <K, V> TreeMultimap<K, V, LinkedList<V>>
	newLinkedListInstance() {
		return new TreeMultimap<K, V, LinkedList<V>>(
				new Factory<LinkedList<V>>() {
					public LinkedList<V> newInstance() {
						return new LinkedList<V>();
					}
				});
	}
	
	/**
	 * Constructs a HashMultimap which is a copy of an existing one.
	 * 
	 * The Map and the C collections are all copied; that is as deep as the
	 * copy goes.
	 */
	public static <K, V, C extends Collection<V>> TreeMultimap<K, V, C>
	newCopyInstance(TreeMultimap<K, V, C> src) {
		return new TreeMultimap<K, V, C>(src);
	}
}
