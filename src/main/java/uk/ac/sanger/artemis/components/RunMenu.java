/* RunMenu.java
 *
 * created: Fri Jan 22 1999
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998,1999,2000  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/RunMenu.java,v 1.13 2009-05-13 15:53:14 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.InvalidKeyException;

import java.io.IOException;
import java.io.StringWriter;
import java.util.Hashtable;
import java.awt.event.*;

import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;


/**
 *  A JMenu of external commands/functions.
 *
 *  @author Kim Rutherford
 *  @version $Id: RunMenu.java,v 1.13 2009-05-13 15:53:14 tjc Exp $
 **/

public class RunMenu extends SelectionMenu 
{
 
  /** */
  private static final long serialVersionUID = 1L;
  private JMenu fastaMenu = null;
  private JMenu fastaMenuOptions = null;
  private Hashtable<String, JMenu> blastMenu = null;
  private Hashtable<String, JMenu> blastMenuOptions = null;

  /**
   *  Create a new RunMenu object.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param menu_name The name of the new menu.
   **/
  public RunMenu(final JFrame frame, final Selection selection,
                 final String menu_name) 
  {
    super(frame, menu_name, selection);

    addNCBISearches(selection);
    addPfamSearches(selection);
    if(Options.isUnixHost())
    {
      addSeparator();
      final ExternalProgramVector external_programs = Options.getOptions()
          .getExternalPrograms();

      boolean sanger_options = Options.getOptions().getPropertyTruthValue(
          "sanger_options");

      final int external_programs_size = external_programs.size();
      for(int i = 0; i < external_programs_size; ++i)
        makeMenuItem(external_programs.elementAt(i), sanger_options);

      addSeparator();

      for(int i = 0; i < external_programs_size; ++i)
        makeOptionsMenuItem(external_programs.elementAt(i));
    }
  }

  /**
   *  Create a new RunMenu object.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   **/
  public RunMenu(final JFrame frame,
                 final Selection selection) 
  {
    this(frame, selection, "Run");
  }

  private void addPfamSearches(final Selection selection)
  {
    final JMenuItem ncbiSearchLinks = new JMenuItem("Pfam Search");
    add(ncbiSearchLinks);
    ncbiSearchLinks.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        final FeatureVector features = selection.getAllFeatures();
        
        if(features.size() != 1)
        {
          JOptionPane.showMessageDialog(RunMenu.this,
              "Please select a single feature to send to Pfam for searching.", 
              "Pfam Search", JOptionPane.INFORMATION_MESSAGE);
          return; 
        }
        final String residues = features.elementAt(0).getTranslation().toString().toUpperCase();

        RunPfamSearchThread pfamSearch = new RunPfamSearchThread(residues);
        pfamSearch.start();
      }
    });
    
    
    final JMenuItem rfam = new JMenuItem("Rfam Search");
    add(rfam);
    rfam.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        final FeatureVector features = selection.getAllFeatures();
        
        if(features.size() != 1)
        {
          JOptionPane.showMessageDialog(RunMenu.this,
              "Please select a single feature to send to Rfam for searching.", 
              "Rfam Search", JOptionPane.INFORMATION_MESSAGE);
          return; 
        }
        final String residues = features.elementAt(0).getTranslationBases();

        RunPfamSearchThread pfamSearch = new RunPfamSearchThread(
            residues, RunPfamSearchThread.rfamUrl);
        pfamSearch.start();
      }
    });
  }
  
  /**
   * Add menu for NCBI web searches
   * @param selection
   */
  private void addNCBISearches(final Selection selection)
  {
    final JMenu ncbiSearchLinks = new JMenu("NCBI Searches");
    add(ncbiSearchLinks);
    
    final ExternalProgramVector ncbi_protein = 
      Options.getOptions().getNCBIPrograms();
    
    for(int i = 0; i < ncbi_protein.size(); ++i) 
    {
      final ExternalProgram program = (ExternalProgram)ncbi_protein.elementAt(i);
      final String programName = program.getName();

      final JMenuItem programMenu = new JMenuItem(programName);
      ncbiSearchLinks.add(programMenu);
      programMenu.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent arg0)
        {
          final FeatureVector features = selection.getAllFeatures();
          
          if(features.size() != 1)
          {
            JOptionPane.showMessageDialog(RunMenu.this,
                "Selected a single feature to send to NCBI for searching.", 
                "NCBI Search", JOptionPane.INFORMATION_MESSAGE);
            return; 
          }
          
          final StringWriter writer = new StringWriter();
          try
          {
            if(program.getType() == ExternalProgram.AA_PROGRAM)
              features.elementAt(0).writeAminoAcidsOfFeature(writer);
            else
              features.elementAt(0).writeBasesOfFeature(writer);
            
            writer.close();
            final String data = RunBlastAtNCBI.setData(programName, writer.toString());
            if(data != null)
            {
              RunBlastAtNCBI blastSearch = new RunBlastAtNCBI(data);
              blastSearch.start();
            }
          }
          catch(IOException ioe){}
          //BrowserControl.displayURL(program.getProgramOptions()+residues);
        }
      });
    }
  }
  
  /**
   *  Make a new menu item for running the given ExternalProgram object.
   *  @param program Create two menu items for this program.
   **/
  private void makeMenuItem(final ExternalProgram program,
                            final boolean sanger_options) 
  {
    final JMenuItem new_menu;
    final String program_name = program.getName();  
    
    if(program.getType() == ExternalProgram.AA_PROGRAM ||
       program.getType() == ExternalProgram.DNA_PROGRAM &&
       sanger_options) 
    {
      final String options_string = program.getProgramOptions();

      if(options_string.length() > 0) 
      {
        if(program_name.startsWith("fasta") ||
           program_name.indexOf("blast")>-1)
          new_menu = new JMenuItem(options_string);
        else 
          new_menu = new JMenuItem("Run " + program_name + " (" +
                        options_string + ") on selected features");
      }
      else 
        new_menu =
          new JMenuItem("Run " + program_name + " on selected features");
    } 
    else 
      new_menu =
        new JMenuItem("Run " + program_name + " on selected features");

    new_menu.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(!checkForSelectionFeatures()) 
          return;

        final FeatureVector selection_features =
          getSelection().getAllFeatures();

        try
        {
          final ExternalProgramMonitor monitor =
                        program.run(selection_features, Splash.getLogger());

          if(monitor == null)
            return;

          monitor.addExternalProgramListener(new ExternalProgramListener() 
          {
            public void statusChanged(final ExternalProgramEvent e) 
            {
              if(e.getType() == ExternalProgramEvent.FINISHED) 
                new MessageFrame(e.getMessage()).setVisible(true);
            }
          });
          new Thread(monitor).start();
        } 
        catch(InvalidKeyException e) 
        {
          new MessageDialog(getParentFrame(),
                             "execution failed: " + e.getMessage());
        }
        catch(EntryInformationException e)
        {
          new MessageDialog(getParentFrame(),
                             "execution of " + program_name +
                             " failed because: " + e.getMessage());
        }
        catch(ReadOnlyException e)
        {
          new MessageDialog(getParentFrame(),
                             "execution of " + program_name +
                             " failed because one of the features is " +
                             "read only");
        }
        catch(IOException e)
        {
          new MessageDialog(getParentFrame(),
                             "execution of " + program_name +
                             " failed because of an I/O error: " +
                             e);
        }
        catch(ExternalProgramException e) 
        {
          new MessageDialog(getParentFrame(),
                            "execution of " + program_name +
                            " failed: " + e.getMessage());
        }
      }
    });

    if(program_name.startsWith("fasta"))
    {
      if(fastaMenu == null)
      {
        fastaMenu = new JMenu("Run fasta on selected features against");
        add(fastaMenu);
      }
      fastaMenu.add(new_menu);
    }
    else if(program.getName().indexOf("blast")>-1)
    {
      if(blastMenu == null)
        blastMenu = new Hashtable<String, JMenu>();

      if(!blastMenu.containsKey(program.getName()))
      {
        JMenu topMenu = new JMenu("Run "+program.getName()+
                          " on selected features against");
        blastMenu.put(program.getName(), topMenu);
        add(topMenu);
      }
      
      JMenu topMenu = blastMenu.get(program.getName());
      topMenu.add(new_menu);
    }
    else
      add(new_menu);
  }

  /**
   *  Make a new options menu item for the given ExternalProgram object.
   *  @param program Create two menu items for this program.
   **/
  private void makeOptionsMenuItem(final ExternalProgram program)
  {
    if(!(program.getType() == ExternalProgram.AA_PROGRAM ||
         program.getType() == ExternalProgram.DNA_PROGRAM)) 
      return;

    final JMenuItem new_options_menu;
    final String program_name = program.getName();

    if(program_name.startsWith("fasta"))
    {
      if(fastaMenuOptions == null)
      {
        fastaMenuOptions = new JMenu("Set " + program_name + " options");
        add(fastaMenuOptions);
      }
      new_options_menu = new JMenuItem(program.getProgramOptions());
      fastaMenuOptions.add(new_options_menu);
    }  
    else if(program_name.indexOf("blast")>-1)
    {
      if(blastMenuOptions == null)
        blastMenuOptions = new Hashtable<String, JMenu>();
      
      String menuStr = "Set " + program_name + " options";
      if(!blastMenuOptions.containsKey(menuStr))
      {
        JMenu topMenu = new JMenu(menuStr);
        blastMenuOptions.put(menuStr, topMenu);
        add(topMenu);
      }
      
      JMenu topMenu = blastMenuOptions.get(menuStr);
      new_options_menu = new JMenuItem(program.getProgramOptions());
      topMenu.add(new_options_menu);
    }
    else
    {
      new_options_menu = new JMenuItem("Set " + program_name + " options");
      add(new_options_menu);
    }

    new_options_menu.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event)
      {
        new ExternalProgramOptions(program);
      }
    });
  }

}
