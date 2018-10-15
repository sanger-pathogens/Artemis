package uk.ac.sanger.artemis.components.alignment;

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
        // Generate and register the OSXAdapter, passing it a hash of all the
        // methods we wish to
        // use as delegates for various com.apple.eawt.ApplicationListener
        // methods
        Class<?> bamClass = Class.forName("uk.ac.sanger.artemis.components.alignment.BamFrame");
        BamOSXAdapter.setQuitHandler(this, bamClass.getDeclaredMethod(
            "exitApp", (Class[]) null));
        BamOSXAdapter.setAboutHandler(this, bamClass.getDeclaredMethod("about",
            (Class[]) null));
        // OSXAdapter.setPreferencesHandler(this,
        // splashClass.getDeclaredMethod("preferences", (Class[])null));
        BamOSXAdapter.setFileHandler(this, bamClass.getDeclaredMethod(
            "loadFile", new Class[]
            { String.class }));
      }
      catch (Exception e)
      {
        e.printStackTrace();
      }
    }

    protected void about()
    {
      JOptionPane.showMessageDialog(this,
          "BamView\nthis is free software and is distributed"
              + "\nunder the terms of the GNU General Public License.",
          "About", JOptionPane.INFORMATION_MESSAGE);
    }

    protected void loadFile(final String bamFile)
    {
      this.bamFile = bamFile;
    }

    protected void exitApp()
    {
      int status = JOptionPane.showConfirmDialog(this, "Exit?",
          "BamView", JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
      if(status == JOptionPane.CANCEL_OPTION)
        return;
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
