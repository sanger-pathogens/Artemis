package uk.ac.sanger.artemis.io;

import java.util.concurrent.atomic.AtomicLong;

/**
 * An ID generator for creating line group 
 * unique IDs. We could use a UUID for this but
 * that's twice the size of a Long!
 * <br/>
 * Thread Safe.
 * 
 * @author kp11
 *
 */
public class LineGroupIdGenerator
{
	/** A counter variable */
	private static AtomicLong idCounter = new AtomicLong();

	/**
	 * Increment the counter and return the new value.
	 * @return Long
	 */
	public static Long createID()
	{
	    return idCounter.getAndIncrement();
	}
}
