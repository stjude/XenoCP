package org.stjude.compbio.common.interval.index;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import org.stjude.compbio.common.interval.Interval;
import org.stjude.compbio.common.interval.tree.Entry;

/**
 * Factory class for creating unmodifiable wrappers of interval indexes.
 */
public class IntervalIndexes {
	public static <T extends Comparable<T>, I extends Interval<T>, D>
	IntervalIndex<T, I, D> unmodifiable(IntervalIndex<T, I, D> src) {
		return new NonLinearImpl<T, I, D>(src);
	}
	
	public static <T extends Comparable<T>, I extends Interval<T>, D>
	LinearIntervalIndex<T, I, D> unmodifiable(LinearIntervalIndex<T, I, D> src) {
		return new LinearImpl<T, I, D>(src);
	}
	
	private static class NonLinearImpl
	<T extends Comparable<T>, I extends Interval<T>, D>
	implements IntervalIndex<T, I, D> {
		private IntervalIndex<T, I, D> src;
		public NonLinearImpl(IntervalIndex<T, I, D> src) {
			this.src = src;
		}
		public List<D> query(I interval) {
			return src.query(interval);
		}
		public Set<D> queryDistinct(I interval) {
			return src.queryDistinct(interval);
		}
		public void add(I interval, D data) {
			throw new UnsupportedOperationException("Index is unmodifiable.");
		}
		public boolean remove(I interval, D data) {
			throw new UnsupportedOperationException("Index is unmodifiable.");
		}
		public List<Entry<T, I, D>> getDataEntries(D data) {
			return src.getDataEntries(data);
		}
		public List<Entry<T, I, D>> queryEntries(I interval) {
			return src.queryEntries(interval);
		}
		public void addEntry(Entry<T, I, D> entry) {
			throw new UnsupportedOperationException("Index is unmodifiable.");
		}
		public boolean removeEntry(Entry<T, I, D> entry) {
			throw new UnsupportedOperationException("Index is unmodifiable.");
		}
		public List<Entry<T, I, D>> removeAllEntries(
				Collection<Entry<T, I, D>> entries) {
			throw new UnsupportedOperationException("Index is unmodifiable.");
		}
		public int size() {
			return src.size();
		}
		public Set<D> dataSet() {
			return Collections.unmodifiableSet(src.dataSet());
		}
	}

	private static class LinearImpl
	<T extends Comparable<T>, I extends Interval<T>, D>
	extends NonLinearImpl<T, I, D>
	implements LinearIntervalIndex<T, I, D> {
		private LinearIntervalIndex<T, I, D> src;
		public LinearImpl(LinearIntervalIndex<T, I, D> src) {
			super(src);
			this.src = src;
		}
		public List<D> query(T point) {
			return src.query(point);
		}
		public Set<D> queryDistinct(T point) {
			return src.queryDistinct(point);
		}
		public List<Entry<T, I, D>> queryEntries(T point) {
			return src.queryEntries(point);
		}		
		public Interval<T> getScope() {
			return src.getScope();
		}
	}
}
