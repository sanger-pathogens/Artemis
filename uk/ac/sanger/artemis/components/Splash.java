/* Splash.java
 *
 * created: Wed May 10 2000
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
 *
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/Splash.java,v 1.8 2005-04-01 16:08:23 tjc Exp $
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

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.border.Border;

/**
 *  Base class that creates a generic "Splash Screen"
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: Splash.java,v 1.8 2005-04-01 16:08:23 tjc Exp $
 **/

abstract public class Splash extends JFrame 
{

  /**
   *  Do any necessary cleanup then exit.
   **/
  abstract protected void exit();

  /**
   *  A label for status and error messages.
   **/
  final private JLabel status_line = new JLabel("");

  /**
   *  The program name that was passed to the constructor.
   **/
  private String program_name;

  /**
   *  The program version that was passed to the constructor.
   **/
  private String program_version;

  /**
   *  The JComponent to draw the main splash screen into
   **/
  private JComponent helix_canvas;

  /**
   *  JMenu bar for the main window.
   **/
  private JMenuBar menu_bar;

  protected JMenu file_menu;

  protected JMenu options_menu;

  private String geneticCode;

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
    super(program_title + " " + program_version);

    this.program_name    = program_name;
    this.program_version = program_version;

    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        exit();
      }
    });

    final javax.swing.LookAndFeel look_and_feel =
      javax.swing.UIManager.getLookAndFeel();

    final javax.swing.plaf.FontUIResource font_ui_resource =
      Options.getOptions().getFontUIResource();

    java.util.Enumeration keys = UIManager.getDefaults().keys();
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

    final FontMetrics fm =
      this.getFontMetrics(status_line.getFont());

    final int font_height = fm.getHeight()+10;

    status_line.setMinimumSize(new Dimension(100, font_height));
    status_line.setPreferredSize(new Dimension(100, font_height));

    Border loweredbevel = BorderFactory.createLoweredBevelBorder();
    Border raisedbevel = BorderFactory.createRaisedBevelBorder();
    Border compound = BorderFactory.createCompoundBorder(raisedbevel,loweredbevel);
    status_line.setBorder(compound);

    getContentPane().add(helix_canvas, "Center");
    getContentPane().add(status_line, "South");

    ClassLoader cl = this.getClass().getClassLoader();
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

    final int x = 460;
    final int y = 250;

    setSize(x, y);

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation(new Point((screen.width - getSize().width) / 2,
                          (screen.height - getSize().height) / 2));
  }

  /**
   *  Return a JComponent object that will display a helix and a short
   *  copyright notice.
   **/
  private JComponent makeHelixCanvas() 
  {
    return new JPanel() {
      public void update(Graphics g) 
      {
        paint(g);
      }

      // return the program name and the program mode in one String
//    private String getNameString() 
//    {
//      if(Options.getOptions().isEukaryoticMode())
//        return program_name + "  [Eukaryotic mode]";
//      else
//        return program_name + "  [Prokaryotic mode]";
//      return geneticCode;
//    }

      /**
       *  Draws the splash screen text.
       **/
      public int textPaint(final Graphics g)
      {
        FontMetrics fm = this.getFontMetrics(g.getFont());
        final int font_height = fm.getHeight() + 3;
        g.setColor(Color.black);
        final int left_margin = 5;

        g.drawString(program_name,
                     helix_width + left_margin, font_height);
        g.drawString(program_version,
                      helix_width + left_margin, font_height * 2);
//      if(Options.getOptions().isEukaryoticMode()) 
//        g.drawString("[Eukaryotic mode]",
//                      helix_width + left_margin, font_height * 3);
//      else
//        g.drawString("[Prokaryotic mode]",
//                      helix_width + left_margin, font_height * 3);
        g.drawString(geneticCode,
                     helix_width + left_margin, font_height * 3);

        g.drawString("Copyright 1998 - 2005",
                      helix_width + left_margin, font_height * 9 / 2);
        g.drawString("Genome Research Limited",
                      helix_width + left_margin, font_height * 11 / 2);

        return font_height;
      }

      public void paint(final Graphics g) 
      {
        final boolean simple_splash_screen =
          Options.getOptions().getPropertyTruthValue("simple_splash_screen");

        g.setColor(Color.white);

        g.fillRect(0, 0, this.getSize().width, this.getSize().height);

        if(simple_splash_screen) {
          // java SIGILL bug work-around
          textPaint(g);
          return;
        }

        if(helix == null) 
        {
//        Toolkit toolkit = Toolkit.getDefaultToolkit();
//        final URL helix_url = Splash.class.getResource("/uk.ac.sanger.artemis/helix.gif");
//        helix = toolkit.getImage(helix_url);

          ClassLoader cl = this.getClass().getClassLoader();
          ImageIcon helix_icon = new ImageIcon(cl.getResource("images/helix.gif"));
          helix = helix_icon.getImage();

//        final URL sanger_url =
//          Splash.class.getResource("/uk.ac.sanger.artemis/sanger-centre.gif");
//        sanger = toolkit.getImage(sanger_url);
          ImageIcon sanger_icon = new ImageIcon(cl.getResource("images/sanger-centre.gif"));
          sanger = sanger_icon.getImage();

          tracker = new MediaTracker(this);
          tracker.addImage(helix, 0);
          tracker.addImage(sanger, 1);

          try 
          {
            tracker.waitForAll();
            helix_height = helix.getHeight(this);
            helix_width = helix.getWidth(this);
            sanger_height = sanger.getHeight(this);
          }
          catch(InterruptedException e) 
          {
            return;
          }
        }

        for(int i=0; i*helix_height<=this.getSize().height; ++i) 
          g.drawImage(helix,
                      0,
                      i * helix_height, this);

        final int font_height = textPaint(g);

        int sanger_position = this.getSize().height - sanger_height;

        if(sanger_position > font_height * 5.5) 
          g.drawImage(sanger,
                       helix_width + 5,
                       sanger_position, this);
      }

      MediaTracker tracker = null;

      /**
       *  The image of the Sanger DNA logo.  This is set in paint().
       **/
      private Image helix = null;

      /**
       *  The image of the Sanger logo.  This is set in paint().
       **/
      private Image sanger = null;

      /**
       *  The height of the Sanger logo.  This is set in paint().
       **/
      private int sanger_height;

      /**
       *  The height of the Sanger DNA logo.  This is set in paint().
       **/
      private int helix_height;

      /**
       *  The width of the Sanger DNA logo.  This is set in paint().
       **/
      private int helix_width;
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

    ActionListener menu_listener = null;

    menu_listener = new ActionListener() 
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


//  final JCheckBoxMenuItem enable_euk_mode_item = new JCheckBoxMenuItem(
//                                                    "Eukaryotic Mode");
//  enable_euk_mode_item.setState(Options.getOptions().isEukaryoticMode());
//  enable_euk_mode_item.addItemListener(new ItemListener() 
//  {
//    public void itemStateChanged(ItemEvent event) 
//    {
//      final boolean item_state = enable_euk_mode_item.getState();
//      Options.getOptions().setEukaryoticMode(item_state);
//      helix_canvas.repaint();
//    }
//  });
//  options_menu.add(enable_euk_mode_item);

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

    if(Options.getOptions().getPropertyTruthValue("sanger_options") &&
        Options.getOptions().getProperty("black_belt_mode") != null) 
    {
      final JCheckBoxMenuItem black_belt_mode_item =
        new JCheckBoxMenuItem("Black Belt Mode");
      final boolean state =
        Options.getOptions().isBlackBeltMode();
      black_belt_mode_item.setState(state);
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

    if(Options.isUnixHost()) 
    {
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
          logger.setVisible(false);
        }
      };
      makeMenuItem(options_menu, "Hide Log Window", menu_listener);
    }
  }

  /**
  *
  * Construct menu for genetic code tables.
  *
  */
  protected void makeGeneticCodeMenu(final JMenu options_menu)
  {
    // available genetic codes
    String gcodes[] = Options.getOptions().getOptionValues("genetic_codes").getArray();

    // get the default
    StringVector gcode_default = Options.getOptions().getOptionValues("genetic_code_default");

    // determine default genetic code table
    int default_code = 0;
    if(gcode_default != null)
    {
      String defS = gcode_default.elementAt(0);
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
      final JCheckBoxMenuItem geneCode = new JCheckBoxMenuItem(gc_name);
      gcodeGroup.add(geneCode);
      geneCode.setActionCommand(num);

      geneCode.addItemListener(new ItemListener()
      {
        public void itemStateChanged(ItemEvent event)
        {
          if(geneCode.getState())
          {
            geneticCode = gc_name;
            String tab = "translation_table_"+geneCode.getActionCommand();
            String startCodons = "start_codons_"+geneCode.getActionCommand();

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
                  String cod_plus_aa = options_file_table.elementAt(i);
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
      options_menu.add(geneCode);

      if(i == default_code)
        geneCode.setState(true);
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
             (KeyEvent.VK_O, InputEvent.CTRL_MASK));
  }

  /**
   *  Return a Logger for warnings/errors/messages.
   **/
  public static Logger getLogger() 
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
   *  The Logger that is returned by getLogger().
   **/
  private final static LogViewer logger = new LogViewer();

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
