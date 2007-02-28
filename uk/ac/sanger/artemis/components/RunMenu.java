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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/RunMenu.java,v 1.9 2007-02-28 15:47:56 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.InvalidKeyException;

import java.io.IOException;
import java.awt.event.*;

import javax.swing.*;

/**
 *  A JMenu of external commands/functions.
 *
 *  @author Kim Rutherford
 *  @version $Id: RunMenu.java,v 1.9 2007-02-28 15:47:56 tjc Exp $
 **/

public class RunMenu extends SelectionMenu 
{
 
  /** */
  private static final long serialVersionUID = 1L;
  private JMenu fastaMenu = null;
  private JMenu fastaMenuOptions = null;
  private JMenu blastpMenu = null;
  private JMenu blastpMenuOptions = null;

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

    final ExternalProgramVector external_programs =
            Options.getOptions().getExternalPrograms();

    boolean sanger_options = 
            Options.getOptions().getPropertyTruthValue("sanger_options");

    final int external_programs_size = external_programs.size();
    for(int i = 0; i < external_programs_size; ++i) 
      makeMenuItem(external_programs.elementAt(i), sanger_options);

    addSeparator();

    for(int i = 0; i < external_programs_size; ++i) 
      makeOptionsMenuItem(external_programs.elementAt(i));

//  if(Options.getOptions().getProperty("jcon_min_jobs") != null) 
//  {
//    addSeparator();
//    final JMenuItem jcon_status = new JMenuItem("Show Job Status ...");

//    jcon_status.addActionListener(new ActionListener() 
//    {
//      public void actionPerformed(ActionEvent event) 
//      {
//        try 
//        {
//          final int ids[] = getIds();
//          TaskViewerFrame tvf = new TaskViewerFrame(ids);

//          tvf.setSize(400, 600);
//          tvf.setVisible(true);
//        }
//        catch(Exception e) 
//        {
//          e.printStackTrace();
//          new MessageDialog(frame, "unable to view job status: " + e);
//        }
//      }
//    });

//    add(jcon_status);
//  }
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
           program_name.startsWith("blastp"))
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
    else if(program.getName().startsWith("blastp"))
    {
      if(blastpMenu == null)
      {
        blastpMenu = new JMenu("Run blastp on selected features against");
        add(blastpMenu);
      }
      blastpMenu.add(new_menu);
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
    else if(program_name.startsWith("blastp"))
    {
      if(blastpMenuOptions == null)
      {
        blastpMenuOptions = new JMenu("Set " + program_name + " options");
        add(blastpMenuOptions);
      }
      new_options_menu = new JMenuItem(program.getProgramOptions());
      blastpMenuOptions.add(new_options_menu);
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

  /**
   *  
   **/
//private int[] getIds() 
//{
//  final FeatureVector selected_features = getSelection().getAllFeatures();
//  final Vector ids_vector = new Vector();

//  for(int feature_index = 0; feature_index < selected_features.size();
//      ++feature_index) 
//  {
//    final Feature this_feature = selected_features.elementAt(feature_index);

//    try
//    {
//      final Qualifier job_qualifier =
//        this_feature.getQualifierByName("job");
//      
//      final StringVector values = job_qualifier.getValues();

//      if(values != null && values.size() > 0) 
//      {
//        for(int value_index=0; value_index<values.size();
//            ++value_index) 
//        {
//          final String job_value = values.elementAt(value_index);
//          final StringVector bits = StringVector.getStrings(job_value);

//          if(bits.size() > 4 && bits.elementAt(2).equals("task:")) 
//          {
//            try 
//            {
//              final Integer task_id = Integer.valueOf(bits.elementAt(3));

//              if(!ids_vector.contains(task_id)) 
//                ids_vector.add(task_id);
//            }
//            catch(NumberFormatException e) {}
//          }
//        }
//      }
//    } 
//    catch(InvalidRelationException e) {}
//  }

//  final int[] ids = new int[ids_vector.size()];

//  for(int i=0 ; i<ids.length ; ++i) 
//    ids[i] =((Integer)ids_vector.elementAt(i)).intValue();

//  return ids;
//}
}
