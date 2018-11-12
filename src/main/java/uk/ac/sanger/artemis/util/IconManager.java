package uk.ac.sanger.artemis.util;

import java.awt.Image;
import java.awt.MediaTracker;
import java.awt.Taskbar;
import java.awt.image.BufferedImage;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;

import uk.ac.sanger.artemis.components.Splash;

/**
 * Provides utility methods for setting 
 * an application icon.
 * 
 * @author kp11
 *
 */
public class IconManager
{
	public static final String ARTEMIS_NAME = "Artemis";
	public static final String ACT_NAME = "ACT";
	public static final String DNAPLOTTER_NAME = "DNAPlotter";
	public static final String BAMVIEW_NAME = "BamView";
	
	/*
	 * Private constructor.
	 */
	private IconManager()
	{
		// Do nothing
	}
	
	/**
	 * Set an application icon.
	 * 
	 * @param app String - application name/identifier
	 */
	public static void setApplicationIcon(String app)
	{
		// Load a window icon image
	    //
	    if (Taskbar.isTaskbarSupported())
	    {
		    BufferedImage image = null;
			try
			{
				image = ImageIO.read(Splash.class.getClassLoader().getResource(getApplicationIcon(app)));
				
				if (image != null)
				{
					Taskbar.getTaskbar().setIconImage(image);
				}
			} 
			catch (IOException e1)
			{
				// TODO Log4j does not seem to be set up for DnaPlotter/BamView.
				//logger4j.error("ERROR: Unable to locate application icon");
			}
	    }
		
	}
	
	/**
	 * Set a dock icon.
	 * 
	 * @param frame JFrame - parent frame
	 * @param app String - application name/identifier
	 */
	public static void setDockIcon(JFrame frame, String app) 
	{
		final ClassLoader cl = frame.getClass().getClassLoader();
		ImageIcon icon = new ImageIcon(cl.getResource(getDockIcon(app)));

	    if(icon != null) 
	    {
	      final Image icon_image = icon.getImage();
	      MediaTracker tracker = new MediaTracker(frame);
	      tracker.addImage(icon_image, 0);

	      try
	      {
	        tracker.waitForAll();
	        frame.setIconImage(icon_image);
	      }
	      catch(InterruptedException e) 
	      {
	    	// TODO Log4j does not seem to be set up for DnaPlotter/BamView.
	    	//logger4j.error("ERROR: Unable to locate dock icon");
	      }
	    }
	}
	
	/**
	 * Utility method to determine the application icon file to load.
	 * 
	 * @param app
	 *            String - application name
	 * @return String - icon file name to load or null if not known.
	 */
	public static String getApplicationIcon(String app)
	{
		if (ARTEMIS_NAME.equalsIgnoreCase(app))
		{
			return "images/artemis_icon.png";
		} 
		else if (ACT_NAME.equalsIgnoreCase(app))
		{
			return "images/act_icon.png";
		}
		else if (DNAPLOTTER_NAME.equalsIgnoreCase(app))
		{
			return "images/dnaplotter_icon.png";
		}
		else if (BAMVIEW_NAME.equalsIgnoreCase(app))
		{
			return "images/bamview_icon.png";
		}
		else
		{
			throw new IllegalArgumentException("Application name '" + app + "' not recognised");
		}
	}
	
	/**
	 * Utility method to determine the dock icon file to load.
	 * 
	 * @param app
	 *            String - application name
	 * @return String - icon file name to load or null if not known.
	 */
	public static String getDockIcon(String app)
	{
		if (ARTEMIS_NAME.equalsIgnoreCase(app))
		{
			return "images/icon.gif";
		} 
		else if (ACT_NAME.equalsIgnoreCase(app))
		{
			return "images/act.gif";
		}
		else
		{
			throw new IllegalArgumentException("Application name '" + app + "' not recognised");
		}
	}
}
