package org.stjude.compbio.util;

import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

/**
 * Iterator over any list in reverse order
 * 
 * This will support remove if the underlying ListIterator supports remove.
 * 
 * @param <E> element type
 */
public class ListReverseIterator<E> implements Iterator<E> {
	/**
	 * Underlying (forward) ListIterator
	 */
	private ListIterator<E> listIterator;
	
	/**
	 * Constructs a reverse iterator from a list
	 */
	public ListReverseIterator(List<E> list) {
		this(list.listIterator(list.size()));
	}
	
	/**
	 * Constructs a reverse iterator from a pre-positioned normal ListIterator.
	 * 
	 * Note: this will begin at the current position of the listIterator.
	 */
	public ListReverseIterator(ListIterator<E> listIterator) {
		this.listIterator = listIterator;
	}

	public boolean hasNext() {
		return listIterator.hasPrevious();
	}

	public E next() {
		return listIterator.previous();
	}

	public void remove() {
		listIterator.remove();
	}
}
