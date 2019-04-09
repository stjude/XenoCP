package org.stjude.compbio.sam.block;

import org.stjude.compbio.common.interval.RefInterval;

/**
 * Provider of reference bases
 */
public interface RefBasesFetcher {
	/**
	 * Provides reference bases for an interval
	 * 
	 * @param interval the ref interval being queried
	 * @return bases, or null if they cannot all be fetched for any reason
	 */
	byte[] getRefBases(RefInterval interval);
}
