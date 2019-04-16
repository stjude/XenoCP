package org.stjude.compbio.util;

import java.util.AbstractCollection;
import java.util.Collection;
import java.util.Iterator;

/**
 * Class that provides an unmodifiable Collection view of a Multimap's values.
 * 
 * This is different than doing multimap.values() because values() will
 * return a collection of the collections that are the values of the outer map.
 * 
 * @param <V> the value type
 */
public class MultimapValuesView<V> extends AbstractCollection<V> {
	private final Multimap<?, V, ? extends Collection<V>> multimap;

	public MultimapValuesView(Multimap<?, V, ?> multimap) {
		this.multimap = multimap;
	}

	@Override
	public int size() {
		int out = 0;
		for(Collection<V> coll: multimap.values()) out += coll.size();
		return out;
	}

	@Override
	public boolean isEmpty() {
		return multimap.isEmpty();
	}

	private class IteratorImpl implements Iterator<V> {
		private Iterator<? extends Collection<V>> colls;
		private Iterator<V> vals;
		private V next;
		
		public IteratorImpl() {
			this.colls = multimap.values().iterator();
			this.vals = null;
			this.next = null;
			advance();
		}
		
		private void advance() {
			// Loop until the vals iterator has a next value
			while(vals == null || !vals.hasNext()) {
				// If there are no more colls, then we're done
				if(!colls.hasNext()) {
					next = null;
					return;
				}
				// Otherwise, load the next coll's iterator
				vals = colls.next().iterator();
			}
			// Now, load next with the next vlaue from vals
			next = vals.next();
		}
		
		public boolean hasNext() {
			return next != null;
		}

		public V next() {
			V out = next;
			advance();
			return out;
		}

		public void remove() {
			throw new UnsupportedOperationException();			
		}
	}

	@Override
	public Iterator<V> iterator() {
		return new IteratorImpl();
	}
}
