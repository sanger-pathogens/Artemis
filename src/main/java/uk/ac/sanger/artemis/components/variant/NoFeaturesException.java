package uk.ac.sanger.artemis.components.variant;

/**
 * An exception to be thrown when
 * a contig that has no features
 * is requested to be displayed.
 * 
 * @author kp11
 *
 */
public class NoFeaturesException extends RuntimeException
{
	private static final long serialVersionUID = 1L;

	/**
	 * Constructor.
	 * @param message - reason for exception
	 */
	public NoFeaturesException(String message)
	{
		super(message);
	}	
}
