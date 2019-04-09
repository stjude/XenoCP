package org.stjude.compbio.common.interval.index;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Set;

import org.stjude.compbio.common.interval.LongInterval;
import org.stjude.compbio.common.interval.StdLongInterval;

public class Benchmarker {
	private static class UcscRecord {
		public final LongInterval tx;
		public UcscRecord(LongInterval tx) {
			super();
			this.tx = tx;
		}
	}
	private static class LocIndex {
		LongIntervalIndex<UcscRecord> index;

		public LocIndex(File refFlat) throws IOException {
			this.index = new LongIntervalIndex<UcscRecord>();
			BufferedReader reader = new BufferedReader(new FileReader(refFlat));
			String line;
			while((line = reader.readLine()) != null) {
				String[] parts = line.split("\t");
				if(parts[2].equals("chr1")) {
					LongInterval interval = new StdLongInterval(
							Long.parseLong(parts[4]), Long.parseLong(parts[5]));
					index.add(interval, new UcscRecord(interval));
				}
			}
			reader.close();
		}
		
		public Set<UcscRecord> get(String chr, long pos) {
			return index.queryDistinct(pos);
		}
		
		public Set<UcscRecord> get(String chr, long left, long right) {
			return index.queryDistinct(new StdLongInterval(left, right));
		}
		
	}
	static long t0;
	static long[] times;
	static int timeIndex;
	private static void reset() {
		t0 = System.nanoTime();
		times = new long[3];
		timeIndex = 0;
	}
	private static void tick() {
		times[timeIndex++] = System.nanoTime() - t0;
		t0 = System.nanoTime();
	}
	private static void report(String label) {
		tick();
		Arrays.sort(times);
		System.out.printf("%12.9f  %s\n", ((double)times[1]) / 1000000000L, label);
	}
	public static void main(String[] args) throws IOException {
		File refFlat = new File(args[0]);
		reset();
		LocIndex locIndex = new LocIndex(refFlat);
		tick();
		locIndex = new LocIndex(refFlat);
		tick();
		locIndex = new LocIndex(refFlat);
		report("initial_build");

		reset();
		locIndex.get("chr1", 14600);
		tick();
		locIndex.get("chr1", 14600);
		tick();
		locIndex.get("chr1", 14600);
		report("chr1:14600");

		reset();
		locIndex.get("chr1", 119575000);
		tick();
		locIndex.get("chr1", 119575000);
		tick();
		locIndex.get("chr1", 119575000);
		report("chr1:119575000");

		reset();
		locIndex.get("chr1", 240497500);
		tick();
		locIndex.get("chr1", 240497500);
		tick();
		locIndex.get("chr1", 240497500);
		report("chr1:240497500");

		reset();
		locIndex.get("chr1", 118570000, 120576000);
		tick();
		locIndex.get("chr1", 118570000, 120576000);
		tick();
		locIndex.get("chr1", 118570000, 120576000);
		report("chr1:118570000-120576000");
		
		reset();
		locIndex.get("chr1", 1, 240497500);
		tick();
		locIndex.get("chr1", 1, 240497500);
		tick();
		records = locIndex.get("chr1", 1, 240497500);
		report("chr1:1-240497500");
		
		reset();
		testBuildNQuery();
		tick();
		testBuildNQuery();
		tick();
		testBuildNQuery();
		report("chr1:1-240497500");
		
	}

	static Set<UcscRecord> records;

	private static void testBuildNQuery() {
		LongIntervalIndex<UcscRecord> index = new LongIntervalIndex<UcscRecord>();
		for(UcscRecord record: records) {
			index.add(record.tx, record);
			index.query(record.tx);
		}
	}

}
