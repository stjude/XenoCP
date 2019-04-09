package org.stjude.compbio.util.io;

import java.io.IOException;

/**
 * Exception indicating a formatting problem
 * 
 * This extends IOException because it is pretty much exclusively used in an
 * IO context.  This is somewhat intuitive because getting data that doesn't
 * look like it should is similar to not being able to get the data at all.
 * Also, it makes exception handling easier because most methods that can
 * throw a FormatException can also throw an IOException, so in cases where
 * exceptions are rethrown, you only need to throw IOException.
 */
public class FormatException extends IOException {
	private static final long serialVersionUID = 5958258187398329581L;

	public FormatException() {
		super();
	}

	public FormatException(String message, Throwable cause) {
		super(message, cause);
	}

	public FormatException(String message) {
		super(message);
	}

	public FormatException(Throwable cause) {
		super(cause);
	}

}
