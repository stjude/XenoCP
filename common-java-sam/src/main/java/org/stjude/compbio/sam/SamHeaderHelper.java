package org.stjude.compbio.sam;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Set;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class SamHeaderHelper {
	/**
	 * Reads a header from a file and closes it
	 * @param file the file to read
	 */
	public static SAMFileHeader ingestHeader(File file) {
		SamReader reader = SamReaderFactory
			.makeDefault()
			.validationStringency(ValidationStringency.SILENT)
			.open(file);
		SAMFileHeader fileHeader = reader.getFileHeader();
		try {
			reader.close();
		}
		catch(IOException e) {
			// No-op
		}
		return fileHeader;
	}

	/**
	 * Clones a SAMFileHeader without cloning the SAMSequenceDictionary; the
	 * seq dict instance is shared by original and clone.
	 * 
	 * This saves memory when the sequence dictionary is large, but keep in
	 * mind that the original and clone will share the same instance!
	 * 
	 * @param fileHeader the header to clone
	 * @return the clone
	 */
	public static SAMFileHeader sharedSeqDictHeaderClone(SAMFileHeader fileHeader) {
		// Temporarily take out the sequence dictionary so that it will not be
		// cloned
		SAMSequenceDictionary seqDict = fileHeader.getSequenceDictionary();
		fileHeader.setSequenceDictionary(new SAMSequenceDictionary());
		
		// Clone the header without the seq dict
		SAMFileHeader clone = fileHeader.clone();
		
		// Put the sequence dictionary back in the original, and in the new
		// instance--note that we have made the clone a little shallower by
		// re-using the instance
		fileHeader.setSequenceDictionary(seqDict);
		clone.setSequenceDictionary(seqDict);
		
		return clone;
	}

	/**
	 * Gets program group "tips" (those that are not depended upon by another)
	 */
	public static Collection<SAMProgramRecord> getProgramRecordTips(
			SAMFileHeader header) {
		// Get depended-upon IDs and ID->PG mapping
		Set<String> nonTipIds = new HashSet<String>();
		Map<String, SAMProgramRecord> pgs =
				new HashMap<String, SAMProgramRecord>();
		for(SAMProgramRecord pg: header.getProgramRecords()) {
			pgs.put(pg.getId(), pg);
			String prevId = pg.getPreviousProgramGroupId();
			if(prevId != null) nonTipIds.add(prevId);
		}
		// Remove depended-upon PGs from mapping
		for(String id: nonTipIds) pgs.remove(id);
		return pgs.values();
	}
	
	/**
	 * Adds a new program group to each existing PG that is not depended upon by
	 * another.
	 * 
	 * The name is taken as the given class's name, and the version is taken
	 * from version.properties in the same package as the class.
	 * 
	 * @param header the header to modify
	 * @param idPart part to add to the end of the ID
	 * @param cls the main class
	 * @return old to new id mapping
	 * @throws RuntimeException if the version cannot be read from properties
	 */
	public static Map<String, String> addProgramRecordToTips(
			SAMFileHeader header, String idPart, Class<?> cls) {
		Properties p = new Properties();
		try {
			p.load(cls.getResourceAsStream("version.properties"));
		} catch (IOException e) {
			throw new RuntimeException("Could not load version", e);
		}
		SAMProgramRecord template = newProgramRecord(
				idPart, cls.getName(), p.getProperty("version"), null);
		return addProgramRecordToTips(header, template);
	}
	
	/**
	 * Adds a new program group to each existing PG that is not depended upon by
	 * another.
	 * 
	 * @param header the header to modify
	 * @param template the template record to add
	 * @return old to new id mapping
	 */
	public static Map<String, String> addProgramRecordToTips(
			SAMFileHeader header, SAMProgramRecord template) {
		Map<String, String> out = new HashMap<String, String>();
		for(SAMProgramRecord pg: getProgramRecordTips(header)) {
			// Get prev/new IDs
			String prevId = pg.getId();
			String newId = prevId + "." + template.getId();
			// Create the group
			SAMProgramRecord newPg = new SAMProgramRecord(
					newId, template);
			newPg.setPreviousProgramGroupId(prevId);
			// Add it
			header.addProgramRecord(newPg);
			out.put(prevId, newId);
		}
		return out;
	}
	
	/**
	 * Constructs a program group.
	 * 
	 * @param id PGID
	 * @param progName program name (may be null)
	 * @param progVers program version (may be null)
	 * @param cli command line (may be null)
	 * @return new SAMProgramRecord
	 */
	public static SAMProgramRecord newProgramRecord(String id,
			String progName, String progVers, String cli) {
		SAMProgramRecord pg = new SAMProgramRecord(id);
		if(progName != null) pg.setProgramName(progName);
		if(progVers != null) pg.setProgramVersion(progVers);
		if(cli != null) pg.setCommandLine(cli);
		return pg;
	}

	/**
	 * Adds a suffix to all PGIDs
	 * 
	 * @param header header to update
	 * @param suff suffix to add
	 */
	public static void suffixAllProgramRecordIds(SAMFileHeader header,
			String suff) {
		List<SAMProgramRecord> pgs = new LinkedList<SAMProgramRecord>();
		for(SAMProgramRecord pg: header.getProgramRecords()) {
			pg = new SAMProgramRecord(pg.getId() + suff, pg);
			String prevPgid = pg.getPreviousProgramGroupId();
			if(prevPgid != null) pg.setPreviousProgramGroupId(prevPgid + suff);
			pgs.add(pg);
		}
		header.setProgramRecords(pgs);
	}
	
	/**
	 * Updates a records PGID according to a PGID map from a previous
	 * add-to-tips.
	 * 
	 * @param record the record to update
	 * @param pgidMap the map from old to new PGIDs
	 */
	public static void updatePgid(SAMRecord record,
			Map<String, String> pgidMap) {
		// Get the PGID and make sure there's a mapping
		Object pgid = record.getAttribute(SAMTag.PG.name());
		if(!pgidMap.containsKey(pgid)) return;
		
		// Update
		record.setAttribute(SAMTag.PG.name(), pgidMap.get(pgid));
	}
}
