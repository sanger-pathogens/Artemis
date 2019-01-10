package uk.ac.sanger.artemis.components.alignment;

import java.awt.Desktop;
import java.awt.desktop.AboutEvent;
import java.awt.desktop.AboutHandler;
import java.awt.desktop.OpenFilesEvent;
import java.awt.desktop.OpenFilesHandler;
import java.awt.desktop.QuitEvent;
import java.awt.desktop.QuitHandler;
import java.awt.desktop.QuitResponse;
import java.io.File;
import java.util.List;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

 class BamFrame extends JFrame
  {
    private static final long serialVersionUID = 1L;
    private String bamFile = null;

    BamFrame()
    {
      super();
      if(isMac())
        initMac();
    }
    
    private void initMac()
    {
          try 
          { 
        	// Special changes for Java 9...
            Desktop desktop = Desktop.getDesktop();
            desktop.setAboutHandler(new AboutHandler() 
            {
            	@Override
    			public void handleAbout(AboutEvent e)
    			{
    				about();
    			}
            });
            
            desktop.setQuitHandler(new QuitHandler()
            {
    			@Override
    			public void handleQuitRequestWith(QuitEvent e, QuitResponse response)
    			{
    				exitApp();
    			}
            });
            
            desktop.setOpenFileHandler(new OpenFilesHandler()
            {
    			@Override
    			public void openFiles(OpenFilesEvent e)
    			{
    				List<File> files = e.getFiles();
    				
    				if (files != null && files.size() > 0) 
    				{
    					try 
    					{
    						loadFile(files.get(0).getCanonicalPath());
    					} 
    					catch (Exception ex)
    					{
    						ex.printStackTrace();
    					}
    				}
    			}
    			
    			
            });
            
          } 
          catch (Exception e)
          {
        	  e.printStackTrace();
          }
    }

    /**
     * Display about dialog.
     */
    protected void about()
    {
  		final JOptionPane pane = new JOptionPane("BamView\nthis is free software and is distributed"
  				+ "\nunder the terms of the GNU General Public License.", JOptionPane.INFORMATION_MESSAGE);
  		final JDialog d = pane.createDialog((JFrame) null, "About");
  		d.setLocation(10, 10);
  		d.setVisible(true);
    }

    protected void loadFile(final String bamFile)
    {
      this.bamFile = bamFile;
    }

    protected void exitApp()
    {
      System.exit(0);
    }

    protected String getBamFile()
    {
      return bamFile;
    }
    
    protected static boolean isMac()
    {
      return System.getProperty("mrj.version") != null ||
             System.getProperty("os.name").toLowerCase().indexOf("mac") >= 0;
    }
  }
