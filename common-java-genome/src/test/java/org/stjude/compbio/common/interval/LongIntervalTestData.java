package org.stjude.compbio.common.interval;

import java.util.Arrays;
import java.util.List;

/**
 * Some test LongIntervals
 * 
 *                     1                   2                   3
 * 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0
 *                       |---- inside ---|
 *         |-- abutL --|------- main ------|-- abutR --|
 * |- disjL -|------ ovlp5L -----|------ ovlp5R -----|- disjR -|
 */
public class LongIntervalTestData {
	public final LongInterval main;
	public final LongInterval inside;
	public final LongInterval abutL;
	public final LongInterval abutR;
	public final LongInterval ovlp5L;
	public final LongInterval ovlp5R;
	public final LongInterval disjL;
	public final LongInterval disjR;
	
	public LongIntervalTestData() {
		this(0);
	}
	
	/**
	 * Constructor that takes the normal locations and shifts them all
	 * @param shift shift amount
	 */
	public LongIntervalTestData(long shift) {
		this.main = new StdLongInterval(10, 20).shift(shift);
		this.inside = new StdLongInterval(11, 19).shift(shift);
		this.abutL = new StdLongInterval(4, 10).shift(shift);
		this.abutR = new StdLongInterval(20, 26).shift(shift);
		this.ovlp5L = new StdLongInterval(5, 15).shift(shift);
		this.ovlp5R = new StdLongInterval(15, 25).shift(shift);
		this.disjL = new StdLongInterval(0, 5).shift(shift);
		this.disjR = new StdLongInterval(25, 30).shift(shift);		
	}
	
	/**
	 * Returns all as a list
	 */
	public List<LongInterval> toList() {
		return Arrays.asList(new LongInterval[] {
				main,
				inside,
				abutL,
				abutR,
				ovlp5L,
				ovlp5R,
				disjL,
				disjR
		});
	}
}
