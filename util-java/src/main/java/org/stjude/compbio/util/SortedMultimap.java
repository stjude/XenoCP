package org.stjude.compbio.util;

import java.util.Collection;
import java.util.SortedMap;

/**
 * A multimap where the map implementation is a SortedMap
 *
 * @param <K> key type
 * @param <V> collections' value type
 * @param <C> type of the collection 
 */
public interface SortedMultimap<K, V, C extends Collection<V>>
extends Multimap<K, V, C>, SortedMap<K, C> {

}
