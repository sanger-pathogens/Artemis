/*
 * Copyright (C) 2009  Genome Research Limited
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
 */
package uk.ac.sanger.artemis.circular.digest;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.circular.Block;
import uk.ac.sanger.artemis.circular.DNADraw;
import uk.ac.sanger.artemis.components.ArtemisMain;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.StickyFileChooser;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.OutOfRangeException;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FlowLayout;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionAdapter;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;

/**
 * 
 */
public class CircularGenomeController
{
  private Block lastBlock = null;
  private JPanel gelPanel;
  private int hgt;
  private JFrame frame = new JFrame();
  private boolean methylation = false;

  public CircularGenomeController()
  {
  }

  /**
   * Create in-silico Pulse Field Gel Electrophoresis from a restriction enzyme
   * digest and draw alongside DNAPlotter
   * @param enzymes
   * @param sequenceFiles
   * @param restrictOutputs
   * @param methylation if true then RE recognition sites will not match methylated bases
   * @throws Exception
   */
  protected void setup(String enzymes, 
                       List<File> sequenceFiles,
                       List<File> restrictOutputs,
                       boolean methylation)
      throws Exception
  {
    this.methylation = methylation;
    // add each sequence file to a different entry group
    List<EntryGroup> entries = new Vector<EntryGroup>();
    if (sequenceFiles != null && sequenceFiles.size() > 0)
    {
      for (int i = 0; i < sequenceFiles.size(); i++)
      {
        EntryGroup entryGroup = getEntryGroupFromFile(sequenceFiles.get(i));
        entries.add(entryGroup);
      }
    }
    else
    {
      EntryGroup entryGroup = getEntryGroupFromFile(null);
      File sequenceFile = new File(((File) entryGroup.getSequenceEntry()
          .getRootDocument().getLocation()).getAbsolutePath()
          + File.separator + entryGroup.getSequenceEntry().getName());
      
      if(sequenceFiles == null)
        sequenceFiles = new Vector<File>();
      sequenceFiles.add(sequenceFile);
      entries.add(entryGroup);
    }

    if (enzymes == null)
      enzymes = promptForEnzymes();

    // run restrict
    if(restrictOutputs == null)
    {
      restrictOutputs = new Vector<File>(sequenceFiles.size());
      for (int i = 0; i < sequenceFiles.size(); i++)
      {
        File sequenceFile = sequenceFiles.get(i);
        File restrictOutput = File.createTempFile("restrict_"
            + sequenceFile.getName(), ".txt");
        restrictOutputs.add(restrictOutput);
        runEmbossRestrict(sequenceFile.getAbsolutePath(), enzymes,
            restrictOutput);
      }
    }
    drawResults(restrictOutputs, entries, sequenceFiles, enzymes);
  }

  /**
   * Run the EMBOSS application restrict. This uses the EMBOSS_ROOT property to
   * define the location of EMBOSS.
   * @param fileName
   * @param cgcb
   * @param restrictOutput
   * @throws IOException
   * @throws InterruptedException
   */
  private void runEmbossRestrict(final String fileName, 
                                 final String enzymes,
                                 final File restrictOutput) throws IOException, InterruptedException
  {
    String[] args = { 
        System.getProperty("EMBOSS_ROOT") + "/bin/restrict", fileName, "-auto",
        "-limit", "y", "-enzymes", enzymes, 
        methylation ? "-methylation" : "", 
        "-out",
        restrictOutput.getCanonicalPath() };

    ProcessBuilder pb = new ProcessBuilder(args);
    pb.redirectErrorStream(true);
    Process p = pb.start();
    System.err.print("** Running restrict");
    try
    {
      InputStream is = p.getInputStream();
      int inchar;
      while ((inchar = is.read()) != -1)
      {
        char c = (char) inchar;
        System.err.print(c);
      }
      System.err.println("**");
      p.waitFor();
      System.err.println("Process exited with '" + p.exitValue() + "'");
    }
    catch (InterruptedException exp)
    {
      exp.printStackTrace();
    }
    p.waitFor();
  }

  /**
   * Display the result in DNAPlotter and the virtual digest.
   * 
   * @param restrictOutput
   * @param entryGroup
   * @param fileName
   * @param cgcb
   * @throws FileNotFoundException
   * @throws IOException
   */
  private void drawResults(final List<File> restrictOutputs,
      final List<EntryGroup> entries, final List<File> sequenceFiles,
      final String enzymes) throws FileNotFoundException, IOException
  {
    JTabbedPane tabbedPane = new JTabbedPane();
    gelPanel= new JPanel(new FlowLayout(FlowLayout.CENTER, 0, 0));
    Dimension preferredSize = null;
    for (int i = 0; i < entries.size(); i++)
    {
      final ReportDetails rd = Utils.findCutSitesFromEmbossReport(
          new FileReader(restrictOutputs.get(i).getCanonicalPath()));

      if(rd.cutSites.size() == 0)
      {
        JOptionPane.showMessageDialog(null, 
            "No cut site found for "+sequenceFiles.get(i).getName(),
            "RE Digest Results", JOptionPane.INFORMATION_MESSAGE);
      }

      final DNADraw dna = Utils.createDNADrawFromReportDetails(
          rd, entries.get(i));
      tabbedPane.add(sequenceFiles.get(i).getName(), dna);
      
      hgt = dna.getHeight();
      final InSilicoGelPanel inSilicoGelPanel = new InSilicoGelPanel(rd.length,
          rd.cutSites, hgt, restrictOutputs.get(i),
          sequenceFiles.get(i).getName());
      gelPanel.add(inSilicoGelPanel);
      preferredSize = inSilicoGelPanel.getPreferredSize();
      addMouseListener(rd, dna, inSilicoGelPanel);
    }

    frame.setTitle ("Sandpiper :: "+enzymes);
    addMenuBar(frame);
    Dimension d = frame.getToolkit().getScreenSize();

    JScrollPane jspGel = new JScrollPane(gelPanel);
    
    Dimension dgel = new Dimension(
        preferredSize.width* (entries.size()>1 ? 2 : 1), preferredSize.height);
    jspGel.setPreferredSize(dgel);
    gelPanel.setMinimumSize(dgel);
    gelPanel.setBackground(Color.white);
    
    JScrollPane jspTabbedPane = new JScrollPane(tabbedPane);
    JSplitPane mainPanel = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, false, 
        jspGel, jspTabbedPane);
    mainPanel.setBackground(Color.white);

    JScrollPane jsp = new JScrollPane(mainPanel);
    jsp.getViewport().setBackground(Color.white);
    frame.getContentPane().add(jsp);

    frame.pack();
    frame.setLocation(((int) d.getWidth() - frame.getWidth()) / 4,
        ((int) d.getHeight() - frame.getHeight()) / 2);
    frame.setVisible(true);
  }

  /**
   * Add the menu bar
   * 
   * @param f
   */
  private void addMenuBar(final JFrame f)
  {
    f.addWindowListener(new WindowAdapter()
    {
      public void windowClosing(WindowEvent event)
      {
        exitApp(f);
      }
    });

    JMenuBar menuBar = new JMenuBar();

    JMenu fileMenu = new JMenu("File");
    menuBar.add(fileMenu);
    
    JMenuItem loadExpData = new JMenuItem("Load experimental data...");
    loadExpData.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        StickyFileChooser fileChooser = new StickyFileChooser();
        int ret = fileChooser.showOpenDialog(null);
        if(ret == StickyFileChooser.CANCEL_OPTION)
          return;
        
        File expFile = fileChooser.getSelectedFile();
        FileReader reader;
        try
        {
          reader = new FileReader(expFile);
          final List<FragmentBand> bands = Utils.findCutSitesFromExperiment(reader);
          final InSilicoGelPanel inSilicoGelPanel = new InSilicoGelPanel(
              bands, hgt, expFile, expFile.getName()); 
          gelPanel.add(inSilicoGelPanel);
          frame.validate();
        }
        catch (FileNotFoundException e1)
        {
          // TODO Auto-generated catch block
          e1.printStackTrace();
        }
        catch (IOException e2)
        {
          // TODO Auto-generated catch block
          e2.printStackTrace();
        }
      }
    });
    fileMenu.add(loadExpData);
    fileMenu.addSeparator();
    
    JMenuItem exitMenu = new JMenuItem("Exit");
    exitMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        exitApp(f);
      }
    });
    fileMenu.add(exitMenu);
    
    f.setJMenuBar(menuBar);
  }

  /**
   * 
   * @param f
   */
  private void exitApp(JFrame f)
  {
    f.dispose();
    System.exit(0);
  }

  /**
   * Add mouse lister to highlight bands on the virtual digest when the mouse is
   * over the corresponding feature in the circular plot.
   * 
   * @param rd
   * @param dna
   * @param inSilicoGelPanel
   */
  private void addMouseListener(final ReportDetails rd, final DNADraw dna,
                                final InSilicoGelPanel inSilicoGelPanel)
  {
    final JPopupMenu popup = new JPopupMenu();
    final JMenuBar menuBar = dna.createMenuBar();

    JMenu[] menus = new JMenu[menuBar.getMenuCount()];
    for (int i = 0; i < menuBar.getMenuCount(); i++)
      menus[i] = menuBar.getMenu(i);

    for (int i = 0; i < menus.length; i++)
      popup.add(menus[i]);

    final JMenuItem openArtemis = new JMenuItem("Open in Artemis...");
    openArtemis.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Range range = null;
        try
        {
          if (lastBlock == null)
            range = new Range(1, dna.getArtemisEntryGroup().getBases()
                .getLength());
          else
            range = new Range(lastBlock.getBstart(), lastBlock.getBend());
        }
        catch (OutOfRangeException e1)
        {
          e1.printStackTrace();
          return;
        } 
        final ArtemisMain main_window = new ArtemisMain(null);
        main_window.setVisible(false);
        final EntryGroup entryGroup = dna.getArtemisEntryGroup().truncate(
            range);
        final EntryEdit entryEdit = new EntryEdit(entryGroup);
        entryEdit.setVisible(true);
      }
    });

    popup.add(openArtemis);

    MouseMotionListener mouseMotionListener = new MouseMotionAdapter()
    {
      public void mouseMoved(MouseEvent me)
      {
        List<CutSite> cutSites = rd.cutSites;
        final Block b = dna.getBlockAtLocation(me.getPoint());
        if (b != null && b.isOverMe(me.getX(), me.getY()))
        {
          int bend = b.getBend();
          // int bstart = b.getBstart();

          if (bend == rd.length)
          {
            ((CutSite) cutSites.get(0)).setHighlighted(true);
            for (int i = 1; i < cutSites.size(); i++)
              ((CutSite) cutSites.get(i)).setHighlighted(false);
          }
          else
          {
            for (int i = 0; i < cutSites.size(); i++)
            {
              CutSite cutSite = (CutSite) cutSites.get(i);
              if (bend == cutSite.getFivePrime())
                cutSite.setHighlighted(true);
              else
                cutSite.setHighlighted(false);
            }
          }
        }
        else
        {
          for (int i = 0; i < cutSites.size(); i++)
          {
            CutSite cutSite = (CutSite) cutSites.get(i);
            cutSite.setHighlighted(false);
          }
        }
        inSilicoGelPanel.repaint();
      }
    };
    dna.addMouseMotionListener(mouseMotionListener);

    MouseListener popupListener = new MouseAdapter()
    {
      public void mousePressed(MouseEvent e)
      {
        maybeShowPopup(e);
      }

      public void mouseReleased(MouseEvent e)
      {
        maybeShowPopup(e);
      }

      private void maybeShowPopup(MouseEvent e)
      {
        if (e.isPopupTrigger())
        {
          lastBlock = dna.getBlockAtLocation(e.getPoint());

          if (lastBlock == null)
            openArtemis.setText("Open in Artemis...");
          else
            openArtemis.setText("Open in Artemis [" + lastBlock.getBstart()
                + ".." + lastBlock.getBend() + "]...");

          popup.show(e.getComponent(), e.getX(), e.getY());
        }
      }
    };
    dna.addMouseListener(popupListener);

    MouseMotionListener gelMouseMotionListener = new MouseMotionAdapter()
    {
      Block lastBlock = null;
      Color lastBlockColour = null;
      public void mouseMoved(MouseEvent me)
      {
        final FragmentBand band = inSilicoGelPanel.getBandAtLocation(me.getPoint());
        if(lastBlock != null)
        {
          lastBlock.setColour(lastBlockColour);
          dna.repaint();
        }
        if(band != null)
        {
          CutSite cs = band.bandCutSite;
          lastBlock = dna.getBlockAtBasePosition(cs.getFivePrime());
          lastBlockColour = lastBlock.getColour();
          lastBlock.setColour(Color.GREEN);
          dna.repaint();
        }
        else
          lastBlock = null;
      }
    };
    inSilicoGelPanel.addMouseMotionListener(gelMouseMotionListener);
  }

  /**
   * Create a DNADraw panel from a file
   * 
   * @param dna_current
   * @return
   */
  private static EntryGroup getEntryGroupFromFile(File fileName)
  {
    Options.getOptions();
    final EntryGroup entryGroup;

    try
    {
      Entry entry = getEntry(fileName);

      if (entry.getBases() != null)
        entryGroup = new SimpleEntryGroup(entry.getBases());
      else
        entryGroup = new SimpleEntryGroup();
      entryGroup.add(entry);
      return entryGroup;
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    }
    catch (NoSequenceException e)
    {
      JOptionPane.showMessageDialog(null, "No sequence found!",
          "Sequence Missing", JOptionPane.WARNING_MESSAGE);
    }
    return null;
  }

  /**
   * Return an Artemis entry from a file
   * 
   * @param entryFileName
   * @param entryGroup
   * @return
   * @throws NoSequenceException
   * @throws OutOfRangeException
   */
  private static Entry getEntry(final File entryFileName)
      throws NoSequenceException, OutOfRangeException
  {

    if (entryFileName == null)
    {
      // no file - prompt for a file
      uk.ac.sanger.artemis.components.FileDialogEntrySource entrySource = new uk.ac.sanger.artemis.components.FileDialogEntrySource(
          null, null);
      Entry entry = entrySource.getEntry(true);
      return entry;
    }

    final Document entry_document = DocumentFactory.makeDocument(entryFileName
        .getAbsolutePath());
    final EntryInformation artemis_entry_information = Options
        .getArtemisEntryInformation();

    final uk.ac.sanger.artemis.io.Entry new_embl_entry = EntryFileDialog
        .getEntryFromFile(null, entry_document, artemis_entry_information,
            false);

    if (new_embl_entry == null) // the read failed
      return null;

    Entry entry = null;
    try
    {
      entry = new Entry(new_embl_entry);
    }
    catch (OutOfRangeException e)
    {
      new MessageDialog(null, "read failed: one of the features in "
          + entryFileName + " has an out of range " + "location: "
          + e.getMessage());
    }
    return entry;
  }

  /**
   * Prompt for the enzyme list.
   * @return
   */
  private String promptForEnzymes()
  {
    Box yBox = Box.createVerticalBox();
    JTextField enzymeList = new JTextField("HincII,hinfI,ppiI,hindiii");
    yBox.add(enzymeList);
    JCheckBox methylationCheckBox = new JCheckBox(
        "RE sites will not match methylated bases", methylation);
    yBox.add(methylationCheckBox);
    
    JOptionPane.showMessageDialog(null, yBox, "Enzyme", JOptionPane.QUESTION_MESSAGE);
    methylation = methylationCheckBox.isSelected();
    return enzymeList.getText().trim();
  }

  public static void main(String args[])
  {
    if (System.getProperty("EMBOSS_ROOT") == null)
    {
      String embossRoot = JOptionPane.showInputDialog(null,
          "Input the EMBOSS installation directory", "/usr/local/emboss");
      System.setProperty("EMBOSS_ROOT", embossRoot.trim());
    }

    String enzymes = null;
    final CircularGenomeController controller = new CircularGenomeController();
    boolean methylation = false;

    List<File> fileNames = null;
    List<File> restrictOutputs = null;
    if (args != null && args.length > 0)
    {
      if (args.length == 1)
      {
        if (args[0].startsWith("-h"))
        {
          System.out.println("-h\t\tshow help");
          System.out
              .println("-enz\t\tcomma separated list of digest enzymes (optional)");
          System.out
              .println("-seq\t\tspace separated list of sequences (optional)");
          System.out
              .println("-methylation\tif this is set then RE recognition sites "
                  + "will not match methylated bases.");
          System.out
              .println("-restrict\tspace separated lists of EMBOSS restrict output "
                  + "in the same order as the sequences (optional).");
          System.exit(0);
        }
        fileNames = new Vector<File>();
        fileNames.add(new File(args[0]));
      }

      for (int i = 0; i < args.length; i++)
      {
        if (args[i].startsWith("-enz"))
          enzymes = args[i + 1];
        else if (args[i].startsWith("-meth"))
          methylation = true;
        else if (args[i].startsWith("-seq"))
        {
          if (fileNames == null)
            fileNames = new Vector<File>();

          for (int j = i + 1; j < args.length; j++)
          {
            if (args[j].startsWith("-"))
              break;
            fileNames.add(new File(args[j]));
          }
        }
        else if (args[i].startsWith("-restrict"))
        {
          if (restrictOutputs == null)
            restrictOutputs = new Vector<File>();

          for (int j = i + 1; j < args.length; j++)
          {
            if (args[j].startsWith("-"))
              break;
            restrictOutputs.add(new File(args[j]));
          }
        }
      }
    }

    final FileSelectionPanel selectionPanel = new FileSelectionPanel(enzymes,
        fileNames, restrictOutputs, methylation);
    final JFrame f = new JFrame("Options and File Selection");
    ActionListener displayButtonListener = new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        try
        {
          f.setCursor(new Cursor(Cursor.WAIT_CURSOR));
          if(selectionPanel.getEmbossRootField() != null)
            System.getProperties().put("EMBOSS_ROOT",
                selectionPanel.getEmbossRootField().getText().trim());
          controller.setup(selectionPanel.getEnzymes(), 
                           selectionPanel.getSequenceFiles(), 
                           selectionPanel.getRestrictOutputs(),
                           selectionPanel.isMethylation());
          f.dispose();
        }
        catch (Exception ex)
        {
          ex.printStackTrace();
        }
        finally
        {
          f.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        }
      }
    };
    selectionPanel.showJFrame(f, displayButtonListener);

  }
};
