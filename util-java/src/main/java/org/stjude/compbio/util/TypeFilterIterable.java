package org.stjude.compbio.util;

import java.util.Iterator;

/**
 * Wraps an Iterable and provides access to only those elements that are of the
 * specified type or a subtype.
 * 
 * Note: the "right" way to do this is not to use this class but rather the 
 * visitor pattern.  However, the visitor pattern can be harder to understand
 * for those not familiar with it, and can add some complexity, so this class is
 * provided.
 *
 * @param <T> the type to return (subtypes are also returned)
 */
public class TypeFilterIterable<T> implements Iterable<T> {
	/**
	 * The underlying iterable being filtered
	 */
	private Iterable<? super T> iterable;
	/**
	 * Class object for the type
	 */
	private Class<T> type;
	
	/**
	 * Constructs a filtering iterable over the underlying iterable.
	 * 
	 * @param iterable underlying iterable to filter
	 * @param type class object for the type
	 */
	public TypeFilterIterable(Iterable<? super T> iterable, Class<T> type) {
		this.iterable = iterable;
		this.type = type;
	}

	public Iterator<T> iterator() {
		return new TypeFilterIterator();
	}

	private class TypeFilterIterator implements Iterator<T> {
		private Iterator<? super T> iterator;
		private T buffer;
		
		public TypeFilterIterator() {
			iterator = iterable.iterator();
			loadBuffer();
		}
		
		@SuppressWarnings("unchecked")
		private void loadBuffer() {
			// Clear first
			buffer = null;
			// If buffer is clear and the underlying iterator has records, then
			// load
			while(buffer == null && iterator.hasNext()) {
				// Load next
				Object next = iterator.next();
				// If type is bad, then re-clear buffer
				if(type.isInstance(next)) buffer = (T)next;
			}
		}
		
		public boolean hasNext() {
			return buffer != null;
		}

		public T next() {
			T out = buffer;
			loadBuffer();
			return out;
		}

		public void remove() {
			throw new UnsupportedOperationException(
					"Cannot remove through TypeFilterIterable wrapper");
		}
		
	}
}
