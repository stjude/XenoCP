package org.stjude.compbio.util;

import java.io.File;
import java.io.IOException;

public class FsHelper {
	/**
	 * Wrapper around File.createTempFile that registers the file for deletion
	 * upon VM exit
	 * 
	 * @see File#createTempFile(String, String)
	 * @return new file that will be removed on exit
	 */
	public static File createAutoRemovedTempFile(String prefix, String suffix)
	throws IOException {
		File file = File.createTempFile(prefix, suffix);
		file.deleteOnExit();
		return file;
	}

	/**
	 * Creates a temporary directory that will be automatically removed when
	 * the VM exits.
	 * 
	 * Note that this method is not safe against race conditions:  it creates a
	 * temporary file in a safe manner, and then removes it and creates a
	 * directory with the same name.  It's possible that between the delete and
	 * mkdir that something else will take that directory.  If that happens,
	 * then an IOException is thrown.
	 * 
	 * @see File#createTempFile(String, String)
	 * @return new dir (will not be removed on exit)
	 */
	public static File createTempDir(String prefix, String suffix)
	throws IOException {
		File file = File.createTempFile(prefix, suffix);
		file.delete();
		if(!file.mkdir()) {
			throw new IOException("Could not create directory " + file);
		}
		return file;
	}
	
	/**
	 * Delete a file or remove a directory and throw an IOException on failure
	 * 
	 * @param file file/dir to delete
	 * @throws IOException if delete() returns false
	 */
	public static void intolerantDelete(File file) throws IOException {
		if(!file.delete()) {
			throw new IOException("Could not delete " + file);
		}
	}
	
	/**
	 * Deletes a file or recursively deletes a directory.
	 * 
	 * @param file file or directory to delete.
	 * @throws IOException on any problem deleting
	 */
	public static void recursiveDelete(File file) throws IOException {
		if(file.isFile()) {
			intolerantDelete(file);
		}
		else {
			for(File child: file.listFiles()) {
				recursiveDelete(child);
			}
		}
	}
}
