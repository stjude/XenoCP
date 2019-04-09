package org.stjude.compbio.sam;

import htsjdk.samtools.SAMRecord;

/**
 * Information used to identify a read: read name and pairing flags
 * 
 * The pairing flags are the flags for the SAMRecord masked by the three
 * bits that pertain to mate pairing: 0x01, 0x20, 0x40.  The flag mask is
 * stored in the MASK constant.
 * 
 * You may test read ids for equality to each other or for equality to a
 * SAMRecord.  A read ID is equal to the SAMRecord with the same name and
 * pairing flags.
 */
public class ReadId {
	/**
	 * The flag mask used
	 */
	public static final int MASK = 0x00C1;
	/**
	 * Flag indicating the read is paired
	 */
	public static final int PAIRED_MASK = 0x0001;
	/**
	 * Masked flags for first in pair
	 */
	public static final int FIRST = 0x0041;
	/**
	 * Masked flags for second in pair
	 */
	public static final int SECOND = 0x0081;
	/**
	 * Flags to XOR in order to get complementary read id
	 */
	public static final int MATE_XOR = 0x00C0;
	
	/**
	 * Read name
	 */
	private String name;
	/**
	 * Pairing-related flags (flags masked by MASK)
	 */
	private int pairingFlags;

	/**
	 * Constructs a ReadId from name and pairing flags
	 * @param name read name
	 * @param pairingFlags flags that are related to pairing
	 */
	public ReadId(String name, int pairingFlags) {
		this.name = name;
		this.pairingFlags = pairingFlags;
	}
	
	/**
	 * Parses a ReadId from a SAMRecord
	 * @param record
	 */
	public ReadId(SAMRecord record) {
		this(record.getReadName(), record.getFlags() & MASK);
	}

	/**
	 * Parses a read name into a ReadId, guessing at mate pairing.
	 * 
	 * It looks for a two-character suffix of [non alphanumeric][12].  If it
	 * finds it, then it assumes paired with the mate in the pair determined by
	 * the last character, otherwise unpaired.
	 * 
	 * @return ReadId
	 */
	public static ReadId parse(String name) {
		if(name.matches(".*[^0-9A-Za-z][12]")) {
			int delimPos = name.length() - 2;
			if(name.substring(delimPos + 1).equals("1")) {
				return new ReadId(name.substring(0, delimPos), FIRST);
			}
			else {
				return new ReadId(name.substring(0, delimPos), SECOND);
			}
		}
		else {
			return new ReadId(name, 0);
		}
	}
	
	/**
	 * Flexible factory method that uses flags if they indicate paired and
	 * parses name otherwise
	 * 
	 * @param record SAM record with pairing indicated in flags OR name
	 * @return new ReadId instance
	 */
	public static ReadId newInstanceFlex(SAMRecord record) {
		if((record.getFlags() & MASK) == 0) {
			return parse(record.getReadName());
		}
		else {
			return new ReadId(record);
		}
	}
	
	public String getName() {
		return name;
	}

	public int getPairingFlags() {
		return pairingFlags;
	}
	
	public boolean isPaired() {
		return (pairingFlags & PAIRED_MASK) > 0;
	}
	
	public boolean isFirst() {
		return pairingFlags == FIRST;
	}

	public boolean isSecond() {
		return pairingFlags == SECOND;
	}
	
	/**
	 * Returns a new ReadId that represents this one's mate
	 * 
	 * @return new ReadId with same name but different flags
	 */
	public ReadId getMate() {
		return new ReadId(name, pairingFlags ^ MATE_XOR);
	}
	
	/**
	 * Checks if two records have equal ReadIds
	 */
	public static boolean equals(SAMRecord r1, SAMRecord r2) {
		return new ReadId(r1).equals(new ReadId(r2));
	}

	@Override
	public boolean equals(Object obj) {
		// Easy checks
		if(obj == null) return false;
		if(obj == this) return true;
		
		// Compare as ReadId or SAMRecord if possible
		ReadId readId;
		if(obj instanceof ReadId) {
			readId = (ReadId)obj;
		}
		else if(obj instanceof SAMRecord) {
			readId = new ReadId((SAMRecord)obj);
		}
		// If it's neither, then just return false
		else {
			return false;			
		}
		
		// Do the comparison
		return readId.name.equals(name) && readId.pairingFlags == pairingFlags;
	}

	@Override
	public int hashCode() {
		return name.hashCode() + 37 * pairingFlags;
	}

	@Override
	public String toString() {
		return name + "(" + pairingFlags + ")";
	}
	
	
}
