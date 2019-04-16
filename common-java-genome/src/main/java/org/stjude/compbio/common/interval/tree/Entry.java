package org.stjude.compbio.common.interval.tree;

import org.stjude.compbio.common.interval.Interval;

/**
 * IntervalTree entry
 *
 * @param <T> coordinate type
 * @param <I> interval type
 * @param <D> data type
 */
public class Entry<T extends Comparable<T>, I extends Interval<T>, D> {
	/**
	 * Interval; may not be null
	 */
	private final I interval;
	/**
	 * Data; may be null
	 */
	private final D data;
	
	public Entry(I interval, D data) {
		if(interval == null) {
			throw new NullPointerException("Interval may not be null.");
		}
		this.interval = interval;
		this.data = data;
	}

	public I getInterval() {
		return interval;
	}

	public D getData() {
		return data;
	}

	@Override
	public int hashCode() {
		return ((data == null) ? 0 : data.hashCode()) +
				37 * interval.hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		if(this == obj) return true;
		if(obj == null) return false;
		if(getClass() != obj.getClass()) return false;
		
		Entry<?, ?, ?> other = (Entry<?, ?, ?>) obj;
		if(data == null) {
			if(other.data != null) return false;
		}
		else if (!data.equals(other.data)) {
			return false;
		}
		return interval.equals(other.interval);
	}
}
