package org.stjude.compbio.util;

import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.NoSuchElementException;
import java.util.Queue;

/**
 * Subclass of LinkedHashSet that also implements the Queue interface
 *
 * @param <E> type of element in the set/queue
 */
public class LinkedHashSetQueue<E> extends LinkedHashSet<E> 
implements Queue<E> {
	private static final long serialVersionUID = -3596725693740301423L;
	
	/**
	 * The head element, if known
	 */
	private E head = null;

	public LinkedHashSetQueue() {
		super();
	}

	public LinkedHashSetQueue(Collection<? extends E> c) {
		super(c);
	}

	public LinkedHashSetQueue(int initialCapacity, float loadFactor) {
		super(initialCapacity, loadFactor);
	}

	public LinkedHashSetQueue(int initialCapacity) {
		super(initialCapacity);
	}

	public boolean offer(E e) {
		try {
			return add(e);
		} catch(IllegalStateException ex) {
			return false;
		}
	}

	public E remove() {
		E out = poll();
		if(out == null) throw new NoSuchElementException();
		return out;
	}

	public E poll() {
		// Get the iterator
		Iterator<E> iterator = iterator();
		// If there's a next element, then load it
		if(iterator.hasNext()) {
			head = iterator.next();
		}
		// Otherwise, reset head to null, and return null--we're done
		else {
			head = null;
			return null;
		}
		// There's a head, so that will be the output
		E out = head;
		// Remove it
		iterator.remove();
		// Try to load another one, since we already have the iterator; if that
		// fails, we have to set head to null; this is an optimization for
		// peek()
		head = iterator.hasNext() ? iterator.next() : null;
		// Return the predetermined output
		return out;
	}

	public E element() {
		E out = peek();
		if(out == null) throw new NoSuchElementException();
		return out;
	}

	public E peek() {
		if(head == null && !isEmpty()) {
			head = iterator().next();
		}
		return head;
	}

}
