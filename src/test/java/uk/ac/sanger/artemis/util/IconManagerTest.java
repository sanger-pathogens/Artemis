package uk.ac.sanger.artemis.util;

import static org.junit.Assert.*;

import java.awt.GraphicsEnvironment;
import java.awt.Image;

import javax.swing.JFrame;

import org.junit.Test;

import junit.framework.AssertionFailedError;

/**
 * Unit test for the IconManager class.
 * 
 * @author kp11
 *
 */
public class IconManagerTest
{

	@Test
	public void testSetApplicationIcon()
	{
		
		// Test would only work for Mac!
		
	}
	
	@Test
	public void testSetDockIcon()
	{
		// Ignore test if in headless mode
		if(GraphicsEnvironment.isHeadless())
		      return;
		
		// Create fake app window...
		class IconTestFrame extends JFrame 
		{ 
			private static final long serialVersionUID = 1L;

			public IconTestFrame (String arg)
			{
				super(arg);
			}
		};
		
		
		// Pretend we are on non-Mac system
		
		IconTestFrame win1 = new IconTestFrame("test frame 1");
		
		System.clearProperty("mrj.version");
		System.setProperty("os.name", "unix");
		IconManager.setDockIcon(win1, IconManager.ARTEMIS_NAME);
		
		assertNotNull("Dock image has been set", win1.getIconImage());
		
		
		// Pretend we are on a Mac - dock icon not currently set
		
		IconTestFrame win2 = new IconTestFrame("test frame 2");
		
		System.setProperty("mrj.version", "mrj");
		System.setProperty("os.name", "");
		IconManager.setDockIcon(win2, IconManager.ARTEMIS_NAME);

		assertNull("Dock image not set", win2.getIconImage());
		
		System.clearProperty("mrj.version");
		System.setProperty("os.name", "mac os x");
		IconManager.setDockIcon(win2, IconManager.ARTEMIS_NAME);

		assertNull("Dock image not set", win2.getIconImage());
		
	}
	
	@Test
	public void testGetApplicationIcon()
	{
		assertTrue(IconManager.getApplicationIcon(IconManager.ARTEMIS_NAME).contains("artemis"));
		assertTrue(IconManager.getApplicationIcon(IconManager.ACT_NAME).contains("act"));
		assertTrue(IconManager.getApplicationIcon(IconManager.DNAPLOTTER_NAME).contains("dnaplotter"));
		assertTrue(IconManager.getApplicationIcon(IconManager.BAMVIEW_NAME).contains("bamview"));
		
		try 
		{
			IconManager.getApplicationIcon("bob").contains("artemis");
			throw new AssertionFailedError("Expected IllegalArgumentException for invalid option");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}
	
	@Test
	public void testGetDockIcon()
	{
		assertTrue(IconManager.getDockIcon(IconManager.ARTEMIS_NAME).contains("icon"));
		assertTrue(IconManager.getDockIcon(IconManager.ACT_NAME).contains("act"));
		assertTrue(IconManager.getDockIcon(IconManager.DNAPLOTTER_NAME).contains("dnaplotter"));
		assertTrue(IconManager.getDockIcon(IconManager.BAMVIEW_NAME).contains("bamview"));
		
		try 
		{
			IconManager.getApplicationIcon("bob").contains("artemis");
			throw new AssertionFailedError("Expected IllegalArgumentException for invalid option");
		}
		catch (IllegalArgumentException e)
		{
			// expected
		}
	}
}
