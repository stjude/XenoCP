package org.stjude.compbio.common.interval.index;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;
import org.stjude.compbio.common.interval.LongInterval;
import org.stjude.compbio.common.interval.LongIntervalTestData;
import org.stjude.compbio.common.interval.StdLongInterval;
import org.stjude.compbio.common.interval.tree.Entry;

public class LongIntervalIndexTest {
	private LongIntervalIndex<String> index;
	private LongIntervalTestData low;
	private LongIntervalTestData high;
	private List<LongInterval> intervals;
	private List<LongInterval> intervalsTwice;
	private List<LongInterval> allIntervals;
	
	@Before
	public void setUp() throws Exception {
		// Instantiate the test index
		index = new LongIntervalIndex<String>();
		// Get test data
		intervals = new LinkedList<LongInterval>();
		low = new LongIntervalTestData();
		intervals.addAll(low.toList());
		high = new LongIntervalTestData(100);
		intervals.addAll(high.toList());
		
		// Build intervalsTwice
		intervalsTwice = new LinkedList<LongInterval>(intervals);
		intervalsTwice.addAll(intervals);
		
		// Put all of that in the index, as well as a second copy shifted 200,
		// but with data according to the unshifted one (so there will be
		// multiple entries for each data object)
		allIntervals = new LinkedList<LongInterval>(intervals);
		for(LongInterval interval: intervals) {
			index.add(interval, interval.toString());
			index.add(interval.shift(200), interval.toString());
			allIntervals.add(interval.shift(200));
		}
	}
	
	private void assertQueryResults(List<String> actual,
			LongInterval... expected) {
		// Easiest way to compare is to sort both
		Collections.sort(actual);
		
		// Convert expected to strings and sort
		List<String> expData = new ArrayList<String>(expected.length);
		for(LongInterval interval: expected) {
			expData.add(interval.toString());
		}
		Collections.sort(expData);
		
		assertEquals(expData, actual);
	}

	private void assertQueryResultEntries(
			List<Entry<Long, LongInterval, String>> actual,
			LongInterval... expected) {
		// Make sure all entries are matched interval/data, and convert to an
		// interval-only list
		List<LongInterval> actIvls = new ArrayList<LongInterval>(actual.size());
		for(Entry<Long, LongInterval, String> entry: actual) {
			// Get the interval and add to interval list
			LongInterval interval = entry.getInterval();
			actIvls.add(interval);
			// For the consistency check, you need to shift the interval first
			// if it's the +200 one
			if(interval.getLeft() >= 200) {
				interval = interval.shift(-200);
			}
			assertEquals(entry.getData(), interval.toString());
		}

		// Now sort and compare the lists
		List<LongInterval> expIvls = Arrays.asList(expected);
		Collections.sort(expIvls);
		Collections.sort(actIvls);
		assertEquals(expIvls, actIvls);
	}

	@Test
	public void testQueryT() {
		assertQueryResults(index.query(13),
				low.main, low.inside, low.ovlp5L);
		assertQueryResults(index.query(50)
				);
		assertQueryResults(index.query(105),
				high.abutL, high.disjL, high.ovlp5L);
		assertQueryResults(index.query(230),
				low.disjR);
	}

	@Test
	public void testQueryI() {
		assertQueryResults(index.query(low.main),
				low.main, low.inside, low.ovlp5L, low.ovlp5R);
		assertQueryResults(index.query(new StdLongInterval(50, 60))
				);
		assertQueryResults(index.query(new StdLongInterval(125, 125)),
				high.abutR);
		assertQueryResults(index.query(new StdLongInterval(0, 400)),
				intervalsTwice.toArray(new LongInterval[] {}));
	}
	
	@Test
	public void testQueryDistinctI() {
		Set<String> results = index.queryDistinct(new StdLongInterval(0, 400));
		assertQueryResults(new ArrayList<String>(results),
				intervals.toArray(new LongInterval[] {}));
	}

	@Test
	public void testRemove() {
		index.remove(low.main, low.main.toString());
		assertQueryResults(index.query(low.main),
				low.inside, low.ovlp5L, low.ovlp5R);
	}

	@Test
	public void testGetDataEntries() {
		assertQueryResultEntries(index.getDataEntries(low.main.toString()),
				low.main, low.main.shift(200));
	}

	@Test
	public void testQueryEntriesT() {
		assertQueryResultEntries(index.queryEntries(13),
				low.main, low.inside, low.ovlp5L);
		assertQueryResultEntries(index.queryEntries(50)
				);
		assertQueryResultEntries(index.queryEntries(105),
				high.abutL, high.disjL, high.ovlp5L);
		assertQueryResultEntries(index.queryEntries(230),
				low.disjR.shift(200));
	}

	@Test
	public void testQueryEntriesI() {
		assertQueryResultEntries(index.queryEntries(low.main),
				low.main, low.inside, low.ovlp5L, low.ovlp5R);
		assertQueryResultEntries(index.queryEntries(new StdLongInterval(50, 60))
				);
		assertQueryResultEntries(index.queryEntries(new StdLongInterval(125, 125)),
				high.abutR);
		assertQueryResultEntries(index.queryEntries(new StdLongInterval(0, 400)),
				allIntervals.toArray(new LongInterval[] {}));
	}

	@Test
	public void testRemoveEntry() {
		index.removeEntry(new Entry<Long, LongInterval, String>(
				low.main, low.main.toString()));
		assertQueryResults(index.query(low.main),
				low.inside, low.ovlp5L, low.ovlp5R);
	}

	@Test
	public void testRemoveAllEntries() {
		int origSize = index.size();
		index.removeAllEntries(index.getDataEntries(low.main.toString()));
		assertEquals(origSize - 2, index.size());
		assertQueryResults(index.query(low.main),
				low.inside, low.ovlp5L, low.ovlp5R);
		assertQueryResults(index.query(low.main.shift(200)),
				low.inside, low.ovlp5L, low.ovlp5R);
	}

}
