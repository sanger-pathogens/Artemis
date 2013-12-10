/* Splash.java
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2000,2002  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.EntrySourceVector;
import uk.ac.sanger.artemis.Logger;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.util.InputStreamProgressEvent;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.management.ManagementFactory;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.MediaTracker;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;

import java.util.Enumeration;
import java.util.Properties;

/**
 *  Base class that creates a generic "Splash Screen"
 *  @author Kim Rutherford
 **/

abstract public class Splash extends JFrame
{
  private static final long serialVersionUID = 1L;

  /** Do any necessary cleanup then exit. */
  abstract protected void exit();

  /** A label for status and error messages. */
  final private JLabel status_line = new JLabel("");

  /** The program name that was passed to the constructor. */
  private String program_name;

  /** program version passed to the constructor. */
  private String program_version;

  /** JComponent to draw the main splash screen into */
  private JComponent helix_canvas;

  private JMenuBar menu_bar;
  protected JMenu file_menu;
  protected JMenu options_menu;

  private JCheckBoxMenuItem geneCode[];
  private String geneticCode;

  private static boolean save_wd_properties = false;
  public static boolean save_display_name = false;
  public static boolean save_systematic_names = false;
  /**  The Artemis LogViewer. */
  private final static LogViewer logger = new LogViewer();
  private static String optionsLogString[];
  public static org.apache.log4j.Logger logger4j = 
         org.apache.log4j.Logger.getLogger(Splash.class);
  
  public Splash(final String program_title,
                final String program_name)
  {
    super(program_title);
    initLogger();
    
    logger4j.info(System.getProperty("java.version"));
    logger4j.info(System.getProperty("java.vendor"));
    logger4j.info(System.getProperty("java.home"));
    logger4j.info(System.getProperty("os.name"));
    logger4j.info("Max. Heap Memory / Mb: "+ManagementFactory.getMemoryMXBean().getHeapMemoryUsage().getMax()/1000000);
    logger4j.info("Starting application: "+program_name);
    
    if(optionsLogString != null)
    {
      for(int i=0; i<optionsLogString.length; i++)
        logger4j.info(optionsLogString[i]);
    }
  }
  
  /**
   *  Create a new JFrame for a Splash screen.
   *  @param program_name The full name of the program.
   *  @param program_title The name to use in the title.
   *  @param program_version The version string.
   **/
  public Splash(final String program_name,
                final String program_title,
                final String program_version) 
  {
    this(program_title+" "+program_version, program_name);
    
    this.program_name    = program_name;
    this.program_version = program_version;
    
    final ClassLoader cl = this.getClass().getClassLoader();
    try
    {
      String line;
      InputStream in = cl.getResourceAsStream("etc/versions");
      BufferedReader reader = new BufferedReader(new InputStreamReader(in));
      while((line = reader.readLine()) != null)
      {
        if(line.startsWith(program_title))
        {
          this.program_version = line.substring( program_title.length() ).trim();
          setTitle(program_title+" "+this.program_version);
        }
      }
    }
    catch (Exception ex)
    {
      logger4j.debug(ex.getMessage());
    }

    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        exitApp();
      }
    });

    //final javax.swing.LookAndFeel look_and_feel =
    //  javax.swing.UIManager.getLookAndFeel();

    final javax.swing.plaf.FontUIResource font_ui_resource =
      Options.getOptions().getFontUIResource();

    final Enumeration<Object> keys = UIManager.getDefaults().keys();
    while(keys.hasMoreElements()) 
    {
      Object key = keys.nextElement();
      Object value = UIManager.get(key);
      if(value instanceof javax.swing.plaf.FontUIResource) 
        UIManager.put(key, font_ui_resource);
    }

    getContentPane().setLayout(new BorderLayout());
    makeAllMenus();
    helix_canvas = makeHelixCanvas();
    status_line.setFont(Options.getOptions().getFont());

    final int font_height = 
        this.getFontMetrics(status_line.getFont()).getHeight()+10;

    status_line.setMinimumSize(new Dimension(100, font_height));
    status_line.setPreferredSize(new Dimension(100, font_height));
    status_line.setBorder(BorderFactory.createCompoundBorder(
        BorderFactory.createRaisedBevelBorder(), BorderFactory.createLoweredBevelBorder()));

    getContentPane().add(helix_canvas, "Center");
    getContentPane().add(status_line, "South");

    //ClassLoader cl = this.getClass().getClassLoader();
    ImageIcon icon = new ImageIcon(cl.getResource("images/icon.gif"));

    if(icon != null) 
    {
      final Image icon_image = icon.getImage();
      MediaTracker tracker = new MediaTracker(this);
      tracker.addImage(icon_image, 0);

      try
      {
        tracker.waitForAll();
        setIconImage(icon_image);
      }
      catch(InterruptedException e) 
      {
        // ignore and continue
      }
    }
    
    pack();
    setSize(460, 250);

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation(new Point((screen.width - getSize().width) / 2,
                          (screen.height - getSize().height) / 2));
    registerForMacOSXEvents();
  }

  public static void initLogger()
  {
    logger.log("");
    logger.getFileViewer().setHideOnClose(true);
    final InputStream options_input_stream =
      Splash.class.getResourceAsStream("/etc/log4j.properties");

    if(options_input_stream != null)
    {
      Properties logProperties = new Properties();
      try
      {
        logProperties.load(options_input_stream);
        org.apache.log4j.PropertyConfigurator.configure(logProperties);
      }
      catch(FileNotFoundException e)
      {
      }
      catch(IOException e)
      {
      }
    }
  }

  public static void appendToLog(String s)
  {
    if(optionsLogString == null)
      optionsLogString = new String[]{s};
    else
    {
      String logStringTmp[] = new String[optionsLogString.length];
      for(int i=0; i<logStringTmp.length; i++)
        logStringTmp[i] = optionsLogString[i];
      optionsLogString = new String[logStringTmp.length+1];
      
      for(int i=0; i<logStringTmp.length; i++)
        optionsLogString[i] = logStringTmp[i];
      
      optionsLogString[optionsLogString.length-1] = s;
    }
  }
  
  /**
   * Generic registration with the Mac OS X application menu
   * Checks the platform, then attempts to register with the Apple EAWT
   */
  private void registerForMacOSXEvents()
  {
    if(isWindows())
      setWorkingDirectory();
    if(isMac()) 
    {
      setWorkingDirectory();
      try 
      {
        // Generate and register the OSXAdapter, passing it a hash of all the methods we wish to
        // use as delegates for various com.apple.eawt.ApplicationListener methods
        Class<?> splashClass = Class.forName("uk.ac.sanger.artemis.components.Splash");
        OSXAdapter.setQuitHandler(this, 
            splashClass.getDeclaredMethod("exitApp", (Class[])null));
        OSXAdapter.setAboutHandler(this, 
            splashClass.getDeclaredMethod("about", (Class[])null));
        //OSXAdapter.setPreferencesHandler(this, 
        //  splashClass.getDeclaredMethod("preferences", (Class[])null));
        OSXAdapter.setFileHandler(this, 
            splashClass.getDeclaredMethod("loadFile", new Class[] { String.class }));
      } 
      catch (Exception e)
      {
        logger4j.error("Error while loading the OSXAdapter:");
        logger4j.error(e.getMessage());
      }
    }
    logger4j.info("Working directory: "+System.getProperty("user.dir"));
  }
  
  private boolean isMac() 
  {
    return System.getProperty("mrj.version") != null ||
           System.getProperty("os.name").toLowerCase().indexOf("mac") >= 0;
  }
  
  private boolean isWindows() 
  {  
    String os = System.getProperty("os.name").toLowerCase();
    logger4j.info("os.name "+os);
    return (os.indexOf("win") >= 0);
  }
 
  protected void about()
  {
    ClassLoader cl = this.getClass().getClassLoader();
    ImageIcon icon = new ImageIcon(cl.getResource("images/icon.gif"));

    JOptionPane.showMessageDialog(this,
            getTitle()+ "\nthis is free software and is distributed"+
                        "\nunder the terms of the GNU General Public License.",
            "About", JOptionPane.INFORMATION_MESSAGE,
            icon);
  }
  
  /**
   * Web start properties need to begin with "javaws." or "jnlp.
   */
  protected static void processJnlpAttributes()
  {
    // JNLP properties
    final Properties properties = System.getProperties();
    for(String key : properties.stringPropertyNames())
    {
      if( key.equals("jnlp.black_belt_mode") ||
          key.equals("jnlp.chado") ||
          key.equals("jnlp.offset") ||
          key.equals("jnlp.artemis.environment") ||
          key.equals("jnlp.sanger_options") ||
          key.equals("jnlp.read_only") || 
          key.startsWith("jnlp.bam") ||
          key.startsWith("jnlp.userplot") ||
          key.startsWith("jnlp.loguserplot") ||
          key.startsWith("jnlp.show_") )
        System.setProperty(key.substring(5), System.getProperty(key));
    }
  }
  
  protected void loadFile(final String fileName)
  {
    if(this instanceof ArtemisMain)
      ((ArtemisMain)this).readArgsAndOptions(new String[]{ fileName }, this);
  }


  /**
   *  Return a JComponent object that will display a helix and a short
   *  copyright notice.
   **/
  private JComponent makeHelixCanvas() 
  {
    return new JPanel() 
    {
      /** */
      private static final long serialVersionUID = 1L;

      public void update(Graphics g) 
      {
        paint(g);
      }

      /**
       *  Draws the splash screen text.
       **/
      public int textPaint(final Graphics graphics)
      {
        FontMetrics fm = this.getFontMetrics(graphics.getFont());
        final int font_height = fm.getHeight() + 3;
        
        Graphics2D g = (Graphics2D) graphics;
        g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING,
                           RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        
        g.setColor(Color.black);
        final int left_margin = 150;
        final int yPos = helix_height+5;

        g.drawString(program_name,
                     left_margin, yPos+font_height);
        g.drawString(program_version,
                     left_margin, yPos+(font_height * 2));
        g.drawString(geneticCode,
                     left_margin, yPos+(font_height * 3));

        g.drawString("Copyright 1998 - 2014",
                     left_margin, yPos+(font_height * 9 / 2));
        g.drawString("Genome Research Limited",
                     left_margin, yPos+(font_height * 11 / 2));

        return font_height;
      }

      public void paint(final Graphics g) 
      {
        /*final boolean simple_splash_screen =
          Options.getOptions().getPropertyTruthValue("simple_splash_screen");*/

        g.setColor(Color.white);

        g.fillRect(0, 0, this.getSize().width, this.getSize().height);

        if(helix == null) 
        {
          ClassLoader cl = this.getClass().getClassLoader();
          ImageIcon helix_icon = new ImageIcon(cl.getResource("images/PSUlogo.gif"));
          helix = helix_icon.getImage();

          tracker = new MediaTracker(this);
          tracker.addImage(helix, 0);

          try 
          {
            tracker.waitForAll();
            helix_height = helix.getHeight(this);
            //helix_width = helix.getWidth(this);
          }
          catch(InterruptedException e) 
          {
            return;
          }
        }

        g.drawImage(helix,
                    0, 0, this);
                    //helix_height, this);

        textPaint(g);
      }

      private MediaTracker tracker = null;

      /** Sanger DNA logo.  This is set in paint() */
      private Image helix = null;


      /** height of the Sanger DNA logo.  This is set in paint(). */
      private int helix_height;

    };
  }

  /**
   *  Return the reference of the Label used as a status line.
   **/
  public JLabel getStatusLabel() 
  {
    return status_line;
  }

  /**
   *  The possible sources for reading Entry objects.
   **/
  public EntrySourceVector getEntrySources(final JFrame frame) 
  {
    return Utilities.getEntrySources(frame, stream_progress_listener);
  }

  /**
   *  Return an InputStreamProgressListener which updates the error label with
   *  the current number of chars read while reading
   **/
  public InputStreamProgressListener getInputStreamProgressListener() 
  {
    return stream_progress_listener;
  }

  /**
   *  Force the options files to be re-read and the EntryEdit components to be
   *  redisplayed.
   **/
  private void resetOptions() 
  {
    Options.getOptions().reset();
  }

  /**
   *  Make all the menus and menu items for the main window.  Also sets up
   *  suitable ActionListener objects for each item.
   */
  private void makeAllMenus() 
  {
    menu_bar = new JMenuBar();
    file_menu = new JMenu("File");
    file_menu.setMnemonic(KeyEvent.VK_F);

    options_menu = new JMenu("Options");
    options_menu.setMnemonic(KeyEvent.VK_O);

    menu_bar.add(file_menu);
    menu_bar.add(options_menu);

    setJMenuBar(menu_bar);

    ActionListener menu_listener = new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        resetOptions();
      }
    };
    makeMenuItem(options_menu, "Re-read Options", menu_listener);

    final JCheckBoxMenuItem enable_direct_edit_item =
      new JCheckBoxMenuItem("Enable Direct Editing");
    enable_direct_edit_item.setState(Options.getOptions().canDirectEdit());
    enable_direct_edit_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        final boolean item_state = enable_direct_edit_item.getState();
        Options.getOptions().setDirectEdit(item_state);
      }
    });
    options_menu.add(enable_direct_edit_item);

    options_menu.addSeparator();
    options_menu.add(new JLabel(" --- Genetic Codes Tables ---"));

    makeGeneticCodeMenu(options_menu);
    options_menu.addSeparator();

    final JCheckBoxMenuItem j2ssh_option = new JCheckBoxMenuItem(
                                         "Send Searches via SSH");

    j2ssh_option.setState((System.getProperty("j2ssh") != null));
    j2ssh_option.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        final boolean item_state = j2ssh_option.getState();
        if(item_state) 
          System.setProperty("j2ssh", "");
        else
          System.setProperty("j2ssh", "false");
      }
    });
    options_menu.add(j2ssh_option);
    options_menu.addSeparator();
    

    final JCheckBoxMenuItem autohide_option = new JCheckBoxMenuItem(
                                         "Auto hide scrollbar");
    autohide_option.setState((System.getProperty("autohide") != null));

    autohide_option.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent event)
      {
        final boolean item_state = autohide_option.getState();
        if(item_state)
          System.setProperty("autohide", "");
        else
          System.setProperty("autohide", "false");
      }
    });
    options_menu.add(autohide_option);
    options_menu.addSeparator();

    
    final JCheckBoxMenuItem highlight_active_entry_item =
      new JCheckBoxMenuItem("Highlight Active Entry");
    final boolean highlight_active_entry_state =
      Options.getOptions().highlightActiveEntryFlag();
    highlight_active_entry_item.setState(highlight_active_entry_state);
    highlight_active_entry_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        final boolean item_state = highlight_active_entry_item.getState();
        Options.getOptions().setHighlightActiveEntryFlag(item_state);
      }
    });
    options_menu.add(highlight_active_entry_item);

    if(Options.getOptions().getProperty("black_belt_mode") != null) 
    {
      final JCheckBoxMenuItem black_belt_mode_item =
        new JCheckBoxMenuItem("Black Belt Mode");

      black_belt_mode_item.setState(Options.isBlackBeltMode());
      black_belt_mode_item.addItemListener(new ItemListener() 
      {
        public void itemStateChanged(ItemEvent event) 
        {
          final boolean item_state = black_belt_mode_item.getState();
          if(item_state) 
            Options.getOptions().put("black_belt_mode", "true");
          else 
            Options.getOptions().put("black_belt_mode", "false");
        }
      });
      options_menu.add(black_belt_mode_item);
    }

    options_menu.addSeparator();
    menu_listener = new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        showLog();
      }
    };
    makeMenuItem(options_menu, "Show Log Window", menu_listener);

    menu_listener = new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        Options.getOptions().setProperty("artemis.user.dir.prompt","yes"); 
        setWorkingDirectory();
      }
    };
    makeMenuItem(options_menu, "Set Working Directory...", menu_listener);

    JMenu lafMenu = new JMenu("Look and Feel");
    final LookAndFeelInfo[] lafInfo = UIManager.getInstalledLookAndFeels();
    ButtonGroup group = new ButtonGroup();
    for(int i=0; i<lafInfo.length; i++)
    {
      JCheckBoxMenuItem laf = new JCheckBoxMenuItem(lafInfo[i].getName());
      group.add(laf);
      lafMenu.add(laf);
      final LookAndFeelInfo thisLAF = lafInfo[i];
      
      if(UIManager.getLookAndFeel().getClass().getName().equals(thisLAF.getClassName()))
        laf.setSelected(true);
      
      laf.addItemListener(new ItemListener() 
      {
        public void itemStateChanged(ItemEvent event) 
        {
          if(event.getStateChange() == ItemEvent.DESELECTED)
            return;
          
          Frame[] frames = JFrame.getFrames();
          try
          {
            UIManager.setLookAndFeel(thisLAF.getClassName());
          }
          catch(Exception e1)
          {
            e1.printStackTrace();
            return;
          }

          for(int i = 0; i < frames.length; i++)
          {
            SwingUtilities.updateComponentTreeUI(frames[i]);
            frames[i].repaint();
          }
          logger4j.debug("Set look and feel to: " + thisLAF.getClassName());
        }
      });
    }
    options_menu.add(lafMenu);
    
    // list JFrames that are open
    final JMenu framesMenu = new JMenu("Windows");
    menu_bar.add(framesMenu);
    framesMenu.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent event)
      {
        framesMenu.removeAll();
        Frame[] frames = JFrame.getFrames();
        for(int i=0;i<frames.length;i++)
        {
          if( !(frames[i] instanceof JFrame) ||
              !frames[i].isVisible())
            continue;
          
          final JFrame thisFrame = (JFrame)frames[i];
          JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem(thisFrame.getTitle(), false);
          if(thisFrame.isActive())
            menuItem.setSelected(true);
          framesMenu.add(menuItem);
          
          menuItem.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              thisFrame.toFront();
              thisFrame.requestFocus();
            }
          });
        }
      }
    });
  }

  /**
   *
   * Set the working directory, used by the file manager.
   *
   */
  public static void setWorkingDirectory()
  {
    final JTextField wdir = new JTextField(System.getProperty("user.dir")+"   ");

    if(Options.getOptions() != null &&
       Options.getOptions().getProperty("artemis.user.dir") != null)
      wdir.setText( Options.getOptions().getProperty("artemis.user.dir") );

    if( Options.getOptions().getProperty("artemis.user.dir.prompt") != null &&
       !Options.getOptions().getPropertyTruthValue("artemis.user.dir.prompt") )
      return;
    
    Box bdown   = Box.createVerticalBox();
    Box bacross = Box.createHorizontalBox();
    JButton browse = new JButton("Browse...");
    browse.addActionListener(new ActionListener ()
    {
      public void actionPerformed (ActionEvent e)
      {
        final StickyFileChooser file_dialog = new StickyFileChooser();
        file_dialog.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        final int status = file_dialog.showOpenDialog(null);
        if(status == JFileChooser.APPROVE_OPTION)
          wdir.setText(file_dialog.getSelectedFile().getAbsolutePath());
      }
    });
    bacross.add(wdir);
    bacross.add(browse);
    bdown.add(bacross);

    bacross = Box.createHorizontalBox();
    JCheckBox saveDir = new JCheckBox("Save between sessions", false);
    bacross.add(saveDir);
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);
    
    
    bacross = Box.createHorizontalBox();
    JCheckBox hide = new JCheckBox("Do not show this prompt again", false);
    bacross.add(hide);
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);
    
    
    Object[] possibleValues = { "OK" };
    JOptionPane.showOptionDialog(null,
                               bdown,
                               "Set Working Directory",
                               JOptionPane.DEFAULT_OPTION,
                               JOptionPane.QUESTION_MESSAGE,null,
                               possibleValues, possibleValues[0]);

    if(hide.isSelected() || 
       Options.getOptions().getProperty("artemis.user.dir.prompt") != null )
    {
      Options.getOptions().setProperty("artemis.user.dir.prompt",
          Boolean.toString(!hide.isSelected())); 
    }
    
    if( (new File(wdir.getText().trim())).exists() )
    {
      System.setProperty("user.dir", wdir.getText().trim());
      if(saveDir.isSelected())
        save_wd_properties = true;
    }
  }

  protected static void exitApp()
  {
    if(save_wd_properties    || save_display_name || save_systematic_names || 
       Options.getOptions().getProperty("artemis.user.dir.prompt") != null)
      saveProperties();
    
    // write the user project properties
    ProjectProperty.writeProperties();
    System.exit(0);
  }

  /**
   *
   * Save properties (working directory) between sessions.
   *
   */
  private static void saveProperties()
  {
    String uhome = System.getProperty("user.home");
    String fs = System.getProperty("file.separator");
    String prop = uhome+fs+".artemis_options";

    writeProperties(prop);
  }

  /**
  *
  * Write or re-write properties and insert/update the user.dir property
  * @param jemProp      properties file
  * @param uHome        user working directory
  *
  */
  private static void writeProperties(final String prop)
  {
     File file_txt = new File(prop);
     File file_tmp = new File(prop + ".tmp");
     try
     {
       if(file_txt.exists())
       {
         BufferedReader bufferedreader = new BufferedReader(new FileReader(file_txt));
         BufferedWriter bufferedwriter = new BufferedWriter(new FileWriter(file_tmp));
         String line;
         
         boolean prompt_saved = false;
         
         while ((line = bufferedreader.readLine()) != null)
         {
           if(line.startsWith("artemis.user.dir") && save_wd_properties)
           {
             line = addEscapeChars("artemis.user.dir="+System.getProperty("user.dir"));
             save_wd_properties = false;
           }
           
           if(line.startsWith("artemis.user.dir.prompt"))
           {
             line = addEscapeChars("artemis.user.dir.prompt="+ 
                 Options.getOptions().getProperty("artemis.user.dir.prompt"));
             prompt_saved = true;
           }

           if(line.startsWith("display_name_qualifiers") && save_display_name)
           {
             String str = "display_name_qualifiers=";
             StringVector strs = Options.getOptions().getDisplayQualifierNames();
             for(int i=0; i<strs.size(); i++)
               str = str.concat(" "+strs.get(i));
             line = addEscapeChars(str);
             save_display_name = false;
           }
           
           if(line.startsWith("systematic_name_qualifiers") && save_systematic_names)
           {
             String str = "systematic_name_qualifiers=";
             StringVector strs = Options.getOptions().getSystematicQualifierNames();
             for(int i=0; i<strs.size(); i++)
               str = str.concat(" "+strs.get(i));
             line = addEscapeChars(str);
             save_systematic_names = false;
           }
           bufferedwriter.write(line);
           bufferedwriter.newLine();
         }
         
         if(save_wd_properties)
         {
           bufferedwriter.write(
               addEscapeChars("artemis.user.dir="+System.getProperty("user.dir")));
           bufferedwriter.newLine();
         }

         if(!prompt_saved && Options.getOptions().getProperty("artemis.user.dir.prompt") != null)
         {
           bufferedwriter.write("artemis.user.dir.prompt="+ 
               Options.getOptions().getProperty("artemis.user.dir.prompt"));
           bufferedwriter.newLine();
         }
         
         if(save_display_name)
         {
           String str = "display_name_qualifiers=";
           StringVector strs = Options.getOptions().getDisplayQualifierNames();
           for(int i=0; i<strs.size(); i++)
             str = str.concat(" "+strs.get(i));
           bufferedwriter.write(addEscapeChars(str));
           bufferedwriter.newLine();
         }
         
         if(save_systematic_names)
         {
           String str = "systematic_name_qualifiers=";
           StringVector strs = Options.getOptions().getSystematicQualifierNames();
           for(int i=0; i<strs.size(); i++)
             str = str.concat(" "+strs.get(i));
           bufferedwriter.write(addEscapeChars(str));
           bufferedwriter.newLine();
         }
         
         bufferedreader.close();
         bufferedwriter.close();
         file_txt.delete();
         file_tmp.renameTo(file_txt);
       }
       else // no existing options file
       {
         BufferedWriter bufferedwriter = new BufferedWriter(new FileWriter(file_txt));
         
         if(save_wd_properties)
         {
           bufferedwriter.write(
               addEscapeChars("artemis.user.dir="+System.getProperty("user.dir")));
           bufferedwriter.newLine();
         }
         
         if(save_display_name)
         {
           String str = "display_name_qualifiers=";
           StringVector strs = Options.getOptions().getDisplayQualifierNames();
           for(int i=0; i<strs.size(); i++)
             str = str.concat(" "+strs.get(i));
           bufferedwriter.write(addEscapeChars(str));
           bufferedwriter.newLine();
         }
         
         if(save_systematic_names)
         {
           String str = "systematic_name_qualifiers=";
           StringVector strs = Options.getOptions().getSystematicQualifierNames();
           for(int i=0; i<strs.size(); i++)
             str = str.concat(" "+strs.get(i));
           bufferedwriter.write(addEscapeChars(str));
           bufferedwriter.newLine();
         }
         bufferedwriter.close();
       }
     }
     catch (FileNotFoundException filenotfoundexception)
     {
       System.err.println(prop+" read error");
     }
     catch (IOException e)
     {
       System.err.println(prop+" i/o error");
     }

  }


  /**
  *
  * Add in escape chars (for windows) to the backslash chars
  * @param l    string to insert escape characters to
  *
  */
  private static String addEscapeChars(String l)
  {
    int n = l.indexOf("\\");

    while( n > -1)
    {
      l = l.substring(0,n)+"\\"+l.substring(n,l.length());
      n = l.indexOf("\\",n+2);
    }
    return l;
  }

  /**
  *
  * Construct menu for genetic code tables.
  *
  */
  private void makeGeneticCodeMenu(final JMenu options_menu)
  {
    // available genetic codes
    StringVector v_genetic_codes = Options.getOptions().getOptionValues("genetic_codes");
    String gcodes[] = (String[])v_genetic_codes.toArray(new String[v_genetic_codes.size()]);

    // get the default
    StringVector gcode_default = Options.getOptions().getOptionValues("genetic_code_default");

    // determine default genetic code table
    int default_code = 0;
    if(gcode_default != null)
    {
      String defS = (String)gcode_default.elementAt(0);
      if(defS.length() < 3)
      {
        try
        {
          int num = Integer.parseInt(defS);
          if(num > 0 && num <= gcodes.length)
            default_code = num-1;
          else
            System.err.println(defS+" is not a valid number");
        }
        catch(NumberFormatException nfe)
        {
          System.err.println(defS+" is not a valid number");
        }
      }
    }
    
    ButtonGroup gcodeGroup = new ButtonGroup();
    geneCode = new JCheckBoxMenuItem[gcodes.length];
    
    for(int i = 0; i< gcodes.length; i++)
    {
      if(gcodes[i].equals("-"))
        continue;

      int ind1; 
      while((ind1 = gcodes[i].indexOf("_")) > -1)
        gcodes[i] = gcodes[i].substring(0,ind1) + " " +
                    gcodes[i].substring(ind1+1,gcodes[i].length());

      String num = Integer.toString(i+1);
      final String gc_name = num+". "+gcodes[i];
      geneCode[i] = new JCheckBoxMenuItem(gc_name);
      gcodeGroup.add(geneCode[i]);
      final int menuNum = i;
      geneCode[i].setActionCommand(num);

      geneCode[i].addItemListener(new ItemListener()
      {
        public void itemStateChanged(ItemEvent event)
        {
          if(geneCode[menuNum].getState())
          {
            geneticCode = gc_name;
            String tab = "translation_table_"+geneCode[menuNum].getActionCommand();
            String startCodons = "start_codons_"+geneCode[menuNum].getActionCommand();

            StringVector options_file_table =
                         Options.getOptions().getOptionValues(tab);

            if(options_file_table != null)
            {
              if(options_file_table.size() == 64) 
              {
                StringBuffer sbuff = new StringBuffer();
                for(int i = 0; i < 64; ++i) 
                  sbuff.append(options_file_table.elementAt(i)+" ");

                Options.getOptions().setGeneticCode(sbuff.toString());
              }
              else
              {
                StringVector table = Options.getOptions().getOptionValues("translation_table_1");

                for(int i = 0; i < options_file_table.size(); ++i)
                {
                  String cod_plus_aa = (String)options_file_table.elementAt(i);
//                System.out.println(cod_plus_aa);
                  final int codon_index = Bases.getIndexOfBase(cod_plus_aa.charAt(0)) * 16 +
                                          Bases.getIndexOfBase(cod_plus_aa.charAt(1)) * 4 +
                                          Bases.getIndexOfBase(cod_plus_aa.charAt(2));

//                System.out.println(cod_plus_aa.substring(3)+"  "+codon_index+"  "+
//                                   table.elementAt(codon_index));
                  table.setElementAt(cod_plus_aa.substring(3), codon_index);
                }
                
                StringBuffer sbuff = new StringBuffer();
                for(int i = 0; i < 64; ++i)
                  sbuff.append(table.elementAt(i)+" ");

                Options.getOptions().setGeneticCode(sbuff.toString());
              }

              options_file_table =
                         Options.getOptions().getOptionValues(startCodons);

              if(options_file_table != null)
              {
                StringBuffer sbuff = new StringBuffer();
                for(int i = 0; i < options_file_table.size(); ++i)
                  sbuff.append(options_file_table.elementAt(i)+" ");
                
                Options.getOptions().setProperty("start_codons",sbuff.toString());
              }
            }
           
            AminoAcidSequence.setGeneCode();
 
            if(helix_canvas != null)
              helix_canvas.repaint();
          }
        }
      });
      options_menu.add(geneCode[i]);

      if(i == default_code)
        geneCode[i].setState(true);
    }
  }
  
  /**
   * Set the Genetic Code number
   * @param geneticCodeValue
   */
  public void setTranslationTable(final String geneticCodeValue)
  {
	for(int i=0;i<geneCode.length;i++)
	{
	  if(geneCode[i] == null)
		continue;

	  if( geneCode[i].getActionCommand().equals(geneticCodeValue) && 
	     !geneCode[i].getState() )
	  {  
		int val = JOptionPane.showConfirmDialog(null, 
		    		"Change translation table to:\n"+
		    		geneCode[i].getText(), 
			    	"Confirm Translation Table Change",
		    		JOptionPane.OK_CANCEL_OPTION);
		if(val == JOptionPane.CANCEL_OPTION)
		  return;
		logger4j.debug("SET GENETIC CODE translationTable="
					+ geneticCodeValue);
	    geneCode[i].setState(true);
		break; 
	  }
	}
  }
  
  /**
   *  Make a new menu item in the given menu, with its label given the
   *  String and add the given ActionListener to it.
   */
  protected static void makeMenuItem(JMenu menu, String name,
                                     ActionListener listener) 
  {
    JMenuItem new_item = new JMenuItem(name);
    menu.add(new_item);
    new_item.addActionListener(listener);
    if(name.equals("Open ..."))
      new_item.setAccelerator(KeyStroke.getKeyStroke
             (KeyEvent.VK_O, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask())); // InputEvent.CTRL_MASK));
  }

  /**
   *  Return a Logger for warnings/errors/messages.
   **/
  protected static Logger getLogger() 
  {
    return logger;
  }

  /**
   *  Return the JComponent that the Splash screen is drawing on.
   **/
  public JComponent getCanvas() 
  {
    return helix_canvas;
  }

  /**
   *  Make the LogViewer visible.
   **/
  public static void showLog() 
  {
    logger.setVisible(true);
  }

  /**
   *  An InputStreamProgressListener used to update the error label with the
   *  current number of chars read.
   **/
  private final InputStreamProgressListener stream_progress_listener =
    new InputStreamProgressListener() {
      public void progressMade(final InputStreamProgressEvent event) 
      {
        final int char_count = event.getCharCount();
        if(char_count == -1) 
          getStatusLabel().setText("");
        else 
          getStatusLabel().setText("chars read so far: " + char_count);
      }

      public void progressMade(String progress)
      {
        getStatusLabel().setText(progress);
      }

    };
}
