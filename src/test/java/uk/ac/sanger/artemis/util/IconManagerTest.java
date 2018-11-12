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
		
		IconTestFrame win = new IconTestFrame("test frame");
		
		IconManager.setDockIcon(win, IconManager.ARTEMIS_NAME);
		
		Image icon = win.getIconImage();
		
		assertNotNull("Dock image has been set", icon);
		
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
