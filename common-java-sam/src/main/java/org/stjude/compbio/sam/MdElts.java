package org.stjude.compbio.sam;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;

/**
 * Class with static methods for working with MdElts (similar to Arrays or
 * Collections)
 * 
 */
public class MdElts {
	/**
	 * Parses an MD string into a list of MdElts
	 * 
	 * @param record SAMRecord from which to retrieve MD annotation
	 * @return list of MdElts
	 */
	public static List<MdElt> parse(SAMRecord record) {
		Object mdAttr = record.getAttribute(SAMTag.MD.name());
		if(mdAttr == null) return null;
		return parse((String)mdAttr);
	}

	/**
	 * Parses an MD string into a list of MdElts
	 * 
	 * @param mdStr the MD tag string
	 * @return list of MdElts
	 */
	public static List<MdElt> parse(String mdStr) {
		// Output list to build
		List<MdElt> out = new LinkedList<MdElt>();
		
		while(!mdStr.isEmpty()) {
			// Determine the type of the next element based on the first char
			MdElt.Type type;
			char firstChar = mdStr.charAt(0);
			if(firstChar == '0') {
				type = null;
			}
			else if(firstChar == '^') {
				type = MdElt.Type.DELETION;
			}
			else if('1' <= firstChar && firstChar <= '9') {
				type = MdElt.Type.MATCH;
			}
			else if('A' <= firstChar && firstChar <= 'Z') {
				type = MdElt.Type.MISMATCH;
			}
			else if('a' <= firstChar && firstChar <= 'z') {
				type = MdElt.Type.MISMATCH;
			}
			else {
				throw new IllegalArgumentException("Could not parse " + mdStr + 
						" at char " + firstChar);
			}
			
			// In the case of a 0-length match, just skip the '0'
			if(type == null) {
				mdStr = mdStr.substring(1);
				continue;
			}
			
			// Do the regexp match
			Matcher matcher = type.pattern.matcher(mdStr);
			if(!matcher.find()) throw new RuntimeException("Unexpected ex");
			String piece = matcher.group();
			
			// Cut the matched piece off the string for the next iteration
			mdStr = mdStr.substring(piece.length());
			
			// Add the element
			switch(type) {
			case MATCH:
				out.add(new MdElt(type, Integer.parseInt(piece), null));
				break;
			case DELETION:
				// Cut off the ^, and...
				piece = piece.substring(1);
				// ...fall through
			case MISMATCH:
				out.add(new MdElt(type, piece.length(), piece.getBytes()));
			}
		}
		
		return out;
	}
	
	/**
	 * Returns a sublist of an MdElt list based on an offset in reference bases
	 * 
	 * @param source source list
	 * @param start starting offset in number of reference bases
	 * @return new sublist
	 */
	public static List<MdElt> sublist(List<MdElt> source, int start) {
		// Init
		List<MdElt> out = new ArrayList<MdElt>(source.size());
		Iterator<MdElt> iterator = source.iterator();
		MdElt elt = null;
		
		// Skip over as many whole elements in source as possible
		while(start > 0 && iterator.hasNext()) {
			elt = iterator.next();
			// Check length cases
			// If the element is completely before start position, then just
			// adjust start position and continue
			if(elt.length <= start) {
				start -= elt.length;
			}
			// Otherwise, if the element spans start position, then split it up,
			// and add the split part after start to the output.  Set start to 0
			// to indicate that we have made it to the start
			else {
				out.add(elt.subMdElt(start));
				start = 0;
			}
		}
		
		// Copy remaining elements whole
		while(iterator.hasNext()) out.add(iterator.next());
		
		return out;
	}
	
	/**
	 * Converts a List<MdElt> to an MD string
	 * 
	 * @param list the list of MdElts
	 * @return the string representation
	 * @throws IllegalArgumentException if two MATCH elts are next to each other
	 */
	public static String listToString(List<MdElt> list) {
		MdElt.Type lastType = null;
		StringBuilder sb = new StringBuilder();
		for(MdElt elt: list) {
			// Check for two MATCH in a row
			if(lastType == MdElt.Type.MATCH &&
					elt.type == MdElt.Type.MATCH) {
				throw new IllegalArgumentException(
						"MdElt list had two MATCH in a row");
			}
			// Check for del-mismatch (need to insert 0 match between)
			else if(lastType == MdElt.Type.DELETION &&
					elt.type == MdElt.Type.MISMATCH) {
				sb.append("0");
			}
			
			// Append to string
			sb.append(elt.toString());
			
			// Remember type
			lastType = elt.type;
		}
		return sb.toString();
	}
	
	/**
	 * Builds reference sequence from read bases and MdElt list
	 * 
	 * @param mdElts md elt list
	 * @param readBases read bases
	 * @return reference bases (CIGAR skips get 0 positions in output)
	 */
	public static byte[] getPackedRefBases(List<MdElt> mdElts,
			byte[] readBases) {
		StringBuilder refBases = new StringBuilder();
		int readIndex = 0;
		for(MdElt mdElt: mdElts) {
			if(mdElt.bases != null) {
				refBases.append(new String(mdElt.bases));
			}
			else {
				refBases.append(new String(Arrays.copyOfRange(
						readBases, readIndex, readIndex + mdElt.length)));
			}
			if(mdElt.type != MdElt.Type.DELETION) readIndex += mdElt.length;
		}
		return refBases.toString().getBytes();
	}
	
	/**
	 * Returns a new ListBuilder which may be used to build up a List<MdElt>
	 */
	public static ListBuilder newListBuilder() {
		return new ListBuilder();
	}
	
	/**
	 * Class with functionality for building a List of MdElts, similar to a
	 * StringBuilder
	 *
	 */
	public static class ListBuilder {
		/**
		 * List of MdElts as built up so far, minus the last element
		 */
		private List<MdElt> committed;
		/**
		 * Pending element being built
		 */
		private MdElt pending;
		
		public ListBuilder() {
			this.committed = new LinkedList<MdElt>();
			this.pending = null;
		}
		
		/**
		 * Adds a base of match
		 */
		public void addMatch() {
			addMatch(1);
		}
		
		/**
		 * Adds some more bases of match
		 * @param num number of matching bases
		 */
		public void addMatch(int num) {
			// If there's a pending element, then process it
			if(pending != null) {
				// If it's match, then just add to the amount, and we are done
				if(pending.type == MdElt.Type.MATCH) {
					pending = new MdElt(
							pending.type, pending.length + num, null);
					return;
				}
				// Otherwise, commit it
				committed.add(pending);
			}
			
			// Add a new match element as pending
			pending = new MdElt(MdElt.Type.MATCH, num, null);
		}
		
		/**
		 * Adds a base of mismatch
		 * @param base the reference base
		 */
		public void addMismatch(byte base) {
			addMismatch(new byte[] { base });
		}
		
		/**
		 * Adds some more mismatch bases
		 * @param bases the reference bases that are mismatched
		 */
		public void addMismatch(byte[] bases) {
			addNonMatch(bases, MdElt.Type.MISMATCH, true);
		}
		
		/**
		 * Adds some more deleted bases
		 * @param bases the reference bases that are deleted
		 */
		public void addDeleted(byte[] bases) {
			addNonMatch(bases, MdElt.Type.DELETION, false);
		}
		
		/**
		 * Adds mismatch or deletion bases
		 * @param bases bases to add
		 * @param type type
		 * @param append whether to append to existing element if possible
		 */
		private void addNonMatch(byte[] bases, MdElt.Type type,
				boolean append) {
			// Get the length to be added
			int addLen = bases.length;
			
			// If there's a pending element, then process it
			if(pending != null) {
				// If types match and append is on, then append and return
				if(append && pending.type == type) {
					int oldLen = pending.length;
					int newLen = oldLen + addLen;
					byte[] newBases = new byte[newLen];
					System.arraycopy(pending.bases, 0, newBases, 0, oldLen);
					System.arraycopy(bases, 0, newBases, oldLen, addLen);
					pending = new MdElt(type, newLen, newBases);
					return;
				}
				// Otherwise, commit it
				committed.add(pending);
			}
			
			// Add a new element as pending
			pending = new MdElt(type, addLen, bases.clone());
		}
		
		/**
		 * Returns the MdElt list as built up so far
		 */
		public List<MdElt> getList() {
			List<MdElt> out = new ArrayList<MdElt>(committed.size() + 1);
			out.addAll(committed);
			if(pending != null) out.add(pending);
			return out;
		}
	}
}
