/* WriteMenu.java
 *
 * created: Mon Jan 11 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/WriteMenu.java,v 1.12 2009-05-06 10:53:51 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.ReadOnlyException;

import uk.ac.sanger.artemis.circular.TextFieldInt;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.StreamSequence;
import uk.ac.sanger.artemis.io.FastaStreamSequence;
import uk.ac.sanger.artemis.io.RawStreamSequence;
import uk.ac.sanger.artemis.io.EmblStreamSequence;
import uk.ac.sanger.artemis.io.GenbankStreamSequence;
import uk.ac.sanger.artemis.io.StreamSequenceFactory;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

/**
 *  A menu of commands for writing out protein and bases.
 *
 *  @author Kim Rutherford
 *  @version $Id: WriteMenu.java,v 1.12 2009-05-06 10:53:51 tjc Exp $
 **/
public class WriteMenu extends SelectionMenu 
{
  private static final long serialVersionUID = 1L;

  /**
   *  Create a new WriteMenu component.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param entry_group The EntryGroup object to use when writing a sequence.
   *  @param menu_name The name of the new menu.
   **/
  public WriteMenu(final JFrame frame, final Selection selection,
                   final EntryGroup entry_group, final String menu_name)
  {
    super(frame, menu_name, selection);

    this.entry_group = entry_group;

    final JMenuItem aa_item = new JMenuItem("Amino Acids Of Selected Features");
    aa_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeAminoAcids();
      }
    });

    add(aa_item);

    
    final JMenuItem aa_to_qualifier_item = new JMenuItem("Amino Acids Of Selected Features to Qualifier");
    aa_to_qualifier_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        final FeatureVector features_to_write =
          getSelection().getAllFeatures();

        try
        {
          for(int i = 0; i < features_to_write.size(); ++i) 
          {
            final Feature selection_feature = features_to_write.elementAt(i);
            final String translation_string =
              selection_feature.getTranslation().toString().toUpperCase();
            selection_feature.setQualifier(new Qualifier("translation", translation_string));
          }
        }
        catch(ReadOnlyException e)
        {
          e.printStackTrace();
        }
        catch(EntryInformationException e)
        {
          e.printStackTrace();
        }
      }
    });

    add(aa_to_qualifier_item);
    
    final JMenuItem pir_item =
      new JMenuItem("PIR Database Of Selected Features");
    pir_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        writePIRDataBase();
      }
    });

    add(pir_item);
    addSeparator();

    
    final JMenu bases_menu = new JMenu("Bases Of Selection");
    final JMenuItem raw_bases_item = new JMenuItem("Raw Format");
    raw_bases_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeBasesOfSelection(StreamSequenceFactory.RAW_FORMAT);
      }
    });

    bases_menu.add(raw_bases_item);

    final JMenuItem fasta_bases_item = new JMenuItem("FASTA Format");
    fasta_bases_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeBasesOfSelection(StreamSequenceFactory.FASTA_FORMAT);
      }
    });

    bases_menu.add(fasta_bases_item);

    final JMenuItem embl_bases_item = new JMenuItem("EMBL Format");
    embl_bases_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeBasesOfSelection(StreamSequenceFactory.EMBL_FORMAT);
      }
    });

    bases_menu.add(embl_bases_item);

    final JMenuItem genbank_bases_item = new JMenuItem("Genbank Format");
    genbank_bases_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeBasesOfSelection(StreamSequenceFactory.GENBANK_FORMAT);
      }
    });

    bases_menu.add(genbank_bases_item);
    add(bases_menu);

//  final JMenu exons_menu = new JMenu("Exons Of Selection");
//  final JMenuItem fasta_exons_item = new JMenuItem("Multiple FASTA Format");
//  fasta_exons_item.addActionListener(new ActionListener()
//  {
//    public void actionPerformed(ActionEvent event)
//    {
//      writeExonsOfSelection(StreamSequenceFactory.FASTA_FORMAT);
//    }
//  });
//  exons_menu.add(fasta_exons_item);
//  add(exons_menu);

    final JMenu upstream_bases_menu =
      new JMenu("Upstream Bases Of Selected Features");

    final JMenuItem raw_upstream_bases_item = new JMenuItem("Raw Format");
    raw_upstream_bases_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeUpstreamBases(StreamSequenceFactory.RAW_FORMAT);
      }
    });

    upstream_bases_menu.add(raw_upstream_bases_item);

    final JMenuItem fasta_upstream_bases_item = new JMenuItem("FASTA Format");
    fasta_upstream_bases_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeUpstreamBases(StreamSequenceFactory.FASTA_FORMAT);
      }
    });

    upstream_bases_menu.add(fasta_upstream_bases_item);
    final JMenuItem embl_upstream_bases_item = new JMenuItem("EMBL Format");
    embl_upstream_bases_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        writeUpstreamBases(StreamSequenceFactory.EMBL_FORMAT);
      }
    });

    upstream_bases_menu.add(embl_upstream_bases_item);

    final JMenuItem genbank_upstream_bases_item = new JMenuItem("Genbank Format");
    genbank_upstream_bases_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeUpstreamBases(StreamSequenceFactory.GENBANK_FORMAT);
      }
    });

    upstream_bases_menu.add(genbank_upstream_bases_item);

    add(upstream_bases_menu);

    final JMenu downstream_bases_menu =
      new JMenu("Downstream Bases Of Selected Features");

    final JMenuItem raw_downstream_bases_item = new JMenuItem("Raw Format");
    raw_downstream_bases_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeDownstreamBases(StreamSequenceFactory.RAW_FORMAT);
      }
    });

    downstream_bases_menu.add(raw_downstream_bases_item);

    final JMenuItem fasta_downstream_bases_item = new JMenuItem("FASTA Format");
    fasta_downstream_bases_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeDownstreamBases(StreamSequenceFactory.FASTA_FORMAT);
      }
    });

    downstream_bases_menu.add(fasta_downstream_bases_item);

    final JMenuItem embl_downstream_bases_item = new JMenuItem("EMBL Format");
    embl_downstream_bases_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeDownstreamBases(StreamSequenceFactory.EMBL_FORMAT);
      }
    });

    downstream_bases_menu.add(embl_downstream_bases_item);

    final JMenuItem genbank_downstream_bases_item =
      new JMenuItem("Genbank Format");
    genbank_downstream_bases_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeDownstreamBases(StreamSequenceFactory.GENBANK_FORMAT);
      }
    });

    downstream_bases_menu.add(genbank_downstream_bases_item);

    add(downstream_bases_menu);

    addSeparator();
    
    final JMenuItem write_combo_bases = new JMenuItem("Upstream+Feature+Downstream Bases ...");
    write_combo_bases.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeComboBases();
      }
    });
    add(write_combo_bases);
    
    addSeparator();

    final JMenu write_all_bases_menu = new JMenu("All Bases");

    final JMenuItem write_raw_item = new JMenuItem("Raw Format");
    write_raw_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeAllSequence(StreamSequenceFactory.RAW_FORMAT);
      }
    });

    write_all_bases_menu.add(write_raw_item);

    final JMenuItem write_fasta_item = new JMenuItem("FASTA Format");
    write_fasta_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeAllSequence(StreamSequenceFactory.FASTA_FORMAT);
      }
    });

    write_all_bases_menu.add(write_fasta_item);

    final JMenuItem write_embl_item = new JMenuItem("EMBL Format");
    write_embl_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeAllSequence(StreamSequenceFactory.EMBL_FORMAT);
      }
    });

    write_all_bases_menu.add(write_embl_item);

    final JMenuItem write_genbank_item = new JMenuItem("Genbank Format");
    write_genbank_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        writeAllSequence(StreamSequenceFactory.GENBANK_FORMAT);
      }
    });

    write_all_bases_menu.add(write_genbank_item);

    add(write_all_bases_menu);

    addSeparator();

    final JMenuItem codon_usage_item =
      new JMenuItem("Codon Usage of Selected Features");
    codon_usage_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        writeCodonUsage();
      }
    });

    add(codon_usage_item);
  }

  /**
   *  Create a new WriteMenu component.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param entry_group The EntryGroup object to use when writing a sequence.
   **/
  public WriteMenu(final JFrame frame, final Selection selection,
                   final EntryGroup entry_group) 
  {
    this(frame, selection, entry_group, "Write");
  }
  
  /**
   *  Write a PIR database of the selected features to a file choosen by the
   *  user.
   **/
  private void writePIRDataBase()
  {
    if(!checkForSelectionFeatures()) 
      return;

    final File write_file =
      getWriteFile("Select a PIR output file name ...", "cosmid.pir", null);

    if(write_file == null)
      return;

    try 
    {
      final FileWriter writer = new FileWriter(write_file);

      final FeatureVector features_to_write =
        getSelection().getAllFeatures();

      for(int i = 0; i < features_to_write.size() ; ++i)
      {
        final Feature selection_feature = features_to_write.elementAt(i);
        selection_feature.writePIROfFeature(writer);
      }

      writer.close();
    }
    catch(IOException e) 
    {
      new MessageDialog(getParentFrame(),
                        "error while writing: " + e.getMessage());
    }
  }

  /**
   *  Write the amino acid symbols of the selected features to a file choosen
   *  by the user.
   **/
  private void writeAminoAcids()
  {
    if(!checkForSelectionFeatures())
      return;

    final File write_file =
      getWriteFile("Select an output file name..", "amino_acids", null);

    if(write_file == null)
      return;

    try
    {
      final FileWriter writer = new FileWriter(write_file);

      final FeatureVector features_to_write =
        getSelection().getAllFeatures();

      for(int i = 0; i < features_to_write.size(); ++i) 
      {
        final Feature selection_feature = features_to_write.elementAt(i);
        selection_feature.writeAminoAcidsOfFeature(writer);
      }

      writer.close();
    } 
    catch(IOException e) 
    {
      new MessageDialog(getParentFrame(),
                        "error while writing: " + e.getMessage());
    }
  }

  /**
   *  Write the bases of the selection to a file choosen by the user.
   *  @param output_type One of EMBL_FORMAT, RAW_FORMAT etc.
   **/
  private void writeBasesOfSelection(final int output_type) 
  {
    final MarkerRange marker_range = getSelection().getMarkerRange();

    if(marker_range == null) 
      writeFeatureBases (output_type);
    else  
    {
      final String selection_bases = Strand.markerRangeBases(marker_range);

      writeBases(selection_bases, "selected bases",
                 getSequenceFileName(output_type),
                 output_type);
    }
  }


  /**
   *  Write the bases of the selected features to a file choosen by the user
   *  or show a message and return immediately if there are no selected
   *  features.
   *  @param sequence The sequence to be written.
   *  @param header The header line that will be used on those output formats
   *    that need it.
   *  @param default_output_filename The filename that is passed to
   *    getWriteFile()
   *  @param output_type One of EMBL_FORMAT, RAW_FORMAT etc.
   **/
  private void writeBases(final String sequence, String header,
                          final String default_output_filename,
                          final int output_type) 
  {
    JTextField headerField = null;

    if(output_type == StreamSequenceFactory.FASTA_FORMAT ||
       output_type == StreamSequenceFactory.EMBL_FORMAT)
      headerField = new JTextField(header);

    final File write_file = getWriteFile("Select an output file name ...",
                                    default_output_filename, headerField);

    if(write_file == null)
      return;

    try 
    {
      final FileWriter writer = new FileWriter(write_file);
      final StreamSequence stream_sequence;

      switch(output_type) 
      {
        case StreamSequenceFactory.FASTA_FORMAT:
          if(headerField != null && !headerField.getText().equals(""))
            header = headerField.getText().trim();

          stream_sequence = new FastaStreamSequence(sequence, header);
          break;

        case StreamSequenceFactory.EMBL_FORMAT:
          stream_sequence = new EmblStreamSequence(sequence);

          if(headerField != null && !headerField.getText().equals(""))
          {
            header = "ID   "+headerField.getText().trim();
            header = header.concat("\nFH   Key             "+
                                   "Location/Qualifiers\nFH\n");
            ((EmblStreamSequence)stream_sequence).setHeader(header);
          }

          break;

        case StreamSequenceFactory.GENBANK_FORMAT:
          stream_sequence = new GenbankStreamSequence(sequence);
          break;

        case StreamSequenceFactory.RAW_FORMAT:
        default:
          stream_sequence = new RawStreamSequence(sequence);
          break;
      }

      stream_sequence.writeToStream(writer);
      writer.close();
    } 
    catch(IOException e) 
    {
      new MessageDialog(getParentFrame(),
                        "error while writing: " + e.getMessage());
    }
  }

  /**
  *
  * Write selected exons to a mutiple FASTA file.
  *
  */
/*
  private void writeExonsOfSelection(final int output_type)
  {
    if(!checkForSelectionFeatures())
      return;

    final File write_file =
      getWriteFile("Select an output file name ...", "exons", null);

    if(write_file == null)
      return;

    try
    {
      final FileWriter writer = new FileWriter(write_file);
      final FeatureVector features_to_write =
                           getSelection().getAllFeatures();
      for(int i = 0; i < features_to_write.size(); ++i)
      {
        final Feature selection_feature = features_to_write.elementAt(i);
        final StringBuffer header_buffer = new StringBuffer();

        header_buffer.append(selection_feature.getSystematicName());
        header_buffer.append(" ");
        header_buffer.append(selection_feature.getIDString());
        header_buffer.append(" ");

        final String product = selection_feature.getProductString();

        if(product == null)
          header_buffer.append("undefined product");
        else
          header_buffer.append(product);

        int seg_size = selection_feature.getSegments().size();
        for(int j = 0; j < seg_size; ++j)
        {
          String bases = selection_feature.getSegments().elementAt(j).getBases();
          Range range = selection_feature.getSegments().elementAt(j).getRawRange();

          String s_range = " "+range.getStart()+":"+range.getEnd();
          if(selection_feature.isForwardFeature())
            s_range = s_range+" forward";
          else
            s_range = s_range+" reverse";

          final StreamSequence stream_sequence =
            getStreamSequence(bases,
                              header_buffer.toString()+s_range,
                              output_type);

          stream_sequence.writeToStream(writer);
        }
      }
      writer.close();
    }
    catch(IOException e)
    {
      new MessageDialog(getParentFrame(),
                        "error while writing: " + e.getMessage());
    }
  }
*/

  /**
   *  Write the bases of the selected features to a file choosen by the user
   *  or show a message and return immediately if there are no selected
   *  features.
   *  @param output_type One of EMBL_FORMAT, RAW_FORMAT etc.
   **/
  private void writeFeatureBases(final int output_type) 
  {
    if(!checkForSelectionFeatures())
      return;


    final File write_file =
      getWriteFile("Select an output file name ...", "bases", null);

    if(write_file == null)
      return;

    try 
    {
      final FileWriter writer = new FileWriter(write_file);
      final FeatureVector features_to_write =
        getSelection().getAllFeatures();

      for(int i = 0; i < features_to_write.size(); ++i) 
      {
        final Feature selection_feature = features_to_write.elementAt(i);
        final StringBuffer header_buffer = new StringBuffer();

        header_buffer.append(selection_feature.getSystematicName());
        header_buffer.append(" ");
        header_buffer.append(selection_feature.getIDString());
        header_buffer.append(" ");

        final String product = selection_feature.getProductString();

        if(product == null) 
          header_buffer.append("undefined product");
        else 
          header_buffer.append(product);
        
        header_buffer.append(" ").append(selection_feature.getWriteRange());

        final StreamSequence stream_sequence =
          getStreamSequence(selection_feature.getBases(),
                            header_buffer.toString(),
                            output_type);

        stream_sequence.writeToStream(writer);
      }

      writer.close();
    } 
    catch(IOException e) 
    {
      new MessageDialog(getParentFrame(),
                        "error while writing: " + e.getMessage());
    }
  }

  /**
   *  Return a StreamSequence object for the given sequence string and the
   *  given type.
   *  @param sequence A String containing the sequence.
   *  @param header The header line that will be used on those output formats
   *    that need it.
   *  @param output_type One of EMBL_FORMAT, RAW_FORMAT etc.
   **/
  private StreamSequence getStreamSequence(final String sequence,
                                           final String header,
                                           final int output_type) 
  {
    switch(output_type)
    {
      case StreamSequenceFactory.FASTA_FORMAT:
        return new FastaStreamSequence (sequence, header);

      case StreamSequenceFactory.EMBL_FORMAT:
        return new EmblStreamSequence (sequence);

      case StreamSequenceFactory.GENBANK_FORMAT:
        return new GenbankStreamSequence (sequence);

      case StreamSequenceFactory.RAW_FORMAT:
      default:
        return new RawStreamSequence (sequence);
    }
  }

  /**
   *  Write the upstream bases of the selected features to a file choosen by
   *  the user.  The user can also choose the number of upstream base to
   *  write.
   *  @param output_type One of EMBL_FORMAT, RAW_FORMAT etc.
   **/
  private void writeUpstreamBases(final int output_type) 
  {
    if(!checkForSelectionFeatures()) 
      return;

    final TextRequester text_requester =
      new TextRequester("write how many bases upstream of each feature?",
                        18, "");

    text_requester.addTextRequesterListener(new TextRequesterListener() 
    {
      public void actionPerformed(final TextRequesterEvent event) 
      {
        if(event.getType() == TextRequesterEvent.CANCEL)
          return;

        final String base_count_string = event.getRequesterText().trim();

        if(base_count_string.length() == 0)
        {
          new MessageDialog(getParentFrame(), "no bases written");
          return;
        }

        try 
        {
          final int base_count =
            Integer.valueOf(base_count_string).intValue();

          final File write_file =
            getWriteFile("Select an output file name ...",
                         "upstream_" + getSequenceFileName(output_type), null);

          if(write_file == null)
            return;

          try 
          {
            final FileWriter writer = new FileWriter(write_file);

            final FeatureVector features_to_write =
              getSelection().getAllFeatures();

            for(int i = 0; i < features_to_write.size(); ++i)
            {
              final Feature selection_feature =
                features_to_write.elementAt(i);

              final String sequence_string =
                selection_feature.getUpstreamBases(base_count);


              final String header_line =
                selection_feature.getIDString() + " - " +
                sequence_string.length() + " bases upstream";

              final StreamSequence stream_sequence =
                getStreamSequence(sequence_string,
                                   header_line,
                                   output_type);
              
              stream_sequence.writeToStream(writer);
            }

            writer.close();
          }
          catch(IOException e) 
          {
            new MessageDialog(getParentFrame(),
                              "error while writing: " + e.getMessage());
          }
        } 
        catch(NumberFormatException e)
        {
          new MessageDialog(getParentFrame(),
                            "this is not a number: " + base_count_string);
        }
      }
    });

    text_requester.setVisible(true);
  }
  
  /**
   * Write a combination of upstream + feature + downstream bases for 
   * selected features.
   */
  private void writeComboBases()
  {
    TextFieldInt upstreamBases   = new TextFieldInt();
    upstreamBases.setValue(0);
    Dimension preferredSize = new Dimension(200, 
         upstreamBases.getPreferredSize().height);
    upstreamBases.setPreferredSize(preferredSize);
    TextFieldInt downstreamBases = new TextFieldInt();
    downstreamBases.setValue(0);
    downstreamBases.setPreferredSize(preferredSize);
    TextFieldInt featureBases = new TextFieldInt();
    featureBases.setValue(0);
    featureBases.setPreferredSize(preferredSize);
    JCheckBox allFeatureBases = new JCheckBox("All bases of the feature(s)", false);
    
    int row = 0;
    JPanel gridPanel = new JPanel(new GridBagLayout());
    GridBagConstraints c = new GridBagConstraints();
    c.gridx = 0;
    c.gridy = row;
    c.anchor = GridBagConstraints.WEST;
    
    gridPanel.add(new JLabel("Number of bases upstream of each feature to write:"), c);
    c.gridy = ++row;
    gridPanel.add(upstreamBases, c);
    c.gridy = ++row;
    gridPanel.add(Box.createVerticalStrut(15), c);
    
    c.gridy = ++row;
    gridPanel.add(new JLabel("Number of feature bases to write:"), c);
    c.gridy = ++row;
    gridPanel.add(featureBases, c);
    c.gridy = ++row;
    gridPanel.add(allFeatureBases, c);
    c.gridy = ++row;
    gridPanel.add(Box.createVerticalStrut(15), c);
    
    c.gridy = ++row;
    gridPanel.add(new JLabel("Number of bases downstream of each feature to write:"),c);
    c.gridy = ++row;
    gridPanel.add(downstreamBases, c);
    
    JOptionPane.showMessageDialog(null, gridPanel, 
                   "Bases to write to FASTA file for selected features", 
                   JOptionPane.QUESTION_MESSAGE);
    
    int output_type = StreamSequenceFactory.FASTA_FORMAT;
    final File write_file =
      getWriteFile("Select an output file name ...",
                   "bases_" + getSequenceFileName(output_type), null);
    if(write_file == null)
      return;
    
    try 
    {
      final FileWriter writer = new FileWriter(write_file);
      final FeatureVector features_to_write =
        getSelection().getAllFeatures();

      int upstream_base_count   = upstreamBases.getValue();
      int downstream_base_count = downstreamBases.getValue();
      int feature_base_count    = featureBases.getValue();
      
      for(int i = 0; i < features_to_write.size(); ++i) 
      {
        final Feature selection_feature =
          features_to_write.elementAt(i);
 
        String header_line =
          selection_feature.getIDString();
        StringBuffer sequenceBuff = new StringBuffer();
        
        if(upstream_base_count > 0)
        {
          sequenceBuff.append(
              selection_feature.getUpstreamBases(upstream_base_count));
          header_line = header_line.concat(" : upstream - "+upstream_base_count);
        }
        
        if(allFeatureBases.isSelected())       
        {
          String bases = selection_feature.getBases();
          sequenceBuff.append(bases);
          header_line = header_line.concat(" : complete feature - "+bases.length());
        }
        else if (feature_base_count > 0)
        {
          sequenceBuff.append(
              selection_feature.getBases().substring(0, feature_base_count));
          header_line = header_line.concat(" : feature - "+feature_base_count);
        }
        
        if(downstream_base_count > 0)
        {
          sequenceBuff.append(
            selection_feature.getDownstreamBases(downstream_base_count));
          header_line = header_line.concat(" : downstream - "+downstream_base_count);
        }

        final StreamSequence stream_sequence =
          getStreamSequence(sequenceBuff.toString(),
                            header_line,
                            output_type);
        
        stream_sequence.writeToStream(writer);
      }

      writer.close();
    } 
    catch(IOException e) 
    {
      new MessageDialog(getParentFrame(),
                        "error while writing: " + e.getMessage());
    }
    
    if(write_file == null)
      return;
    
  }

  /**
   *  Write the downstream bases of the selected features to a file choosen by
   *  the user.  The user can also choose the number of downstream base to
   *  write.
   *  @param output_type One of EMBL_FORMAT, RAW_FORMAT etc.
   **/
  private void writeDownstreamBases(final int output_type) 
  {
    if(!checkForSelectionFeatures())
      return;

    final TextRequester text_requester =
      new TextRequester("write how many bases downstream of each feature?",
                         18, "");

    text_requester.addTextRequesterListener(new TextRequesterListener()
    {
      public void actionPerformed(final TextRequesterEvent event) 
      {
        if(event.getType() == TextRequesterEvent.CANCEL)
          return;

        final String base_count_string = event.getRequesterText().trim();

        if(base_count_string.length() == 0) 
        {
          new MessageDialog(getParentFrame(), "no bases written");
          return;
        }

        try 
        {
          final int base_count =
            Integer.valueOf(base_count_string).intValue();

          final File write_file =
            getWriteFile("Select an output file name ...",
                         "downstream_" + getSequenceFileName(output_type), null);

          if(write_file == null)
            return;

          try 
          {
            final FileWriter writer = new FileWriter(write_file);
            final FeatureVector features_to_write =
              getSelection().getAllFeatures();

            for(int i = 0; i < features_to_write.size(); ++i) 
            {
              final Feature selection_feature =
                features_to_write.elementAt(i);

              final String sequence_string =
                selection_feature.getDownstreamBases(base_count);


              final String header_line =
                selection_feature.getIDString() + " - " +
                sequence_string.length() + " bases downstream ";

              final StreamSequence stream_sequence =
                getStreamSequence(sequence_string,
                                  header_line,
                                  output_type);
              
              stream_sequence.writeToStream(writer);
            }

            writer.close();
          } 
          catch(IOException e) 
          {
            new MessageDialog(getParentFrame(),
                              "error while writing: " + e.getMessage());
          }
        }
        catch(NumberFormatException e) 
        {
          new MessageDialog(getParentFrame(),
                            "this is not a number: " + base_count_string);
        }
      }
    });

    text_requester.setVisible(true);
  }

  /**
   *  Helper method for writeFeatureBases(), writeUpstreamBases() and
   *  writeDownstreamBases(). 
   *  @return A string of the form "100:200 reverse" or "1:2222 forward".
   **/
  /*private static String getWriteRange(final Feature feature) 
  {
    return (feature.isForwardFeature() ?
            feature.getFirstCodingBaseMarker().getRawPosition() + ":" +
            feature.getLastBaseMarker().getRawPosition() +
            " forward" :
            feature.getLastBaseMarker().getRawPosition() + ":" +
            feature.getFirstCodingBaseMarker().getRawPosition() +
            " reverse");
  }*/

  
  /**
   *  Return a sensible file name to write sequence of the given type to.
   *  @param output_type One of EMBL_FORMAT, RAW_FORMAT etc.
   **/
  private String getSequenceFileName(final int output_type) 
  {
    switch(output_type) 
    {
      case StreamSequenceFactory.FASTA_FORMAT:
        return "sequence.dna";

      case StreamSequenceFactory.EMBL_FORMAT:
        return "sequence.embl";

      case StreamSequenceFactory.GENBANK_FORMAT:
        return "sequence.genbank";

      case StreamSequenceFactory.RAW_FORMAT:
      default:
        return "sequence.seq";
    }
  }
  
  /**
   *  Write the bases from entry_group as raw sequence.
   *  @param output_type One of EMBL_FORMAT, RAW_FORMAT etc.
   **/
  private void writeAllSequence(final int output_type) 
  {
    final String file_name = getSequenceFileName(output_type);

    String name = "all_bases";
    try
    {
      name = entry_group.getSequenceEntry().getName()+" all_bases";
      if( entry_group.getSequenceEntry().getHeaderText() != null &&
         (entry_group.getSequenceEntry().getHeaderText().startsWith("ID   ") ||
          entry_group.getSequenceEntry().getHeaderText().startsWith("LOCUS       ")) )
      {
        StringTokenizer tok = new StringTokenizer(
            entry_group.getSequenceEntry().getHeaderText(),
            " ;");
        tok.nextElement();

        if(tok.hasMoreTokens())
        {
          String ID = tok.nextToken();
          if(ID.length() > 1 && !(ID.indexOf("\n") > -1) && !(ID.indexOf("\r") > -1))
            name = ID + " all_bases";
        }
      }
    }
    catch(Exception e){}
    writeBases(entry_group.getBases().toString(),
               name, file_name, output_type);
  }

  /**
   *  Popup a requester and ask for a file name to write to.  If the file
   *  exists ask the user whether to overwrite.
   *  @param title The title of the new FileDialog JFrame.
   *  @param default_name The name to put in the FileDialog as a default.
   **/
  private File getWriteFile(final String title, final String default_name,
                            final JTextField idField) 
  {
    final JFileChooser dialog = new StickyFileChooser();

    JCheckBox header = new JCheckBox("Add Header", false);
    if(idField != null)
      dialog.setAccessory(header);

    dialog.setDialogTitle(title);
    dialog.setFileSelectionMode(JFileChooser.FILES_ONLY);
    dialog.setSelectedFile(new File(default_name));
    dialog.setDialogType(JFileChooser.SAVE_DIALOG);
    final int status = dialog.showSaveDialog(getParentFrame());

    if(status != JFileChooser.APPROVE_OPTION)
      return null;

    if(dialog.getSelectedFile() == null) 
      return null;
    else
    {
      final File write_file = dialog.getSelectedFile();

      if(write_file.exists())
      {
        final YesNoDialog yes_no_dialog =
          new YesNoDialog(getParentFrame(),
                          "this file exists: " + write_file +
                          " overwrite it?");
        
        if(!yes_no_dialog.getResult())
          return null;
      }

      // request user to provide header
      if(idField != null)
      {
        if(header.isSelected())
        {
          Box bdown = Box.createVerticalBox();
          bdown.add(idField);

          int n = JOptionPane.showConfirmDialog(null, bdown,
                            "Enter the entry ID",
                            JOptionPane.OK_CANCEL_OPTION,
                            JOptionPane.QUESTION_MESSAGE);

          if(n == JOptionPane.CANCEL_OPTION)
            idField.setText("");
        }
        else
          idField.setText("");
      }
      return write_file;
    }
  }

  /**
   *  Write a table of codon usage for the selected features to a file choosen
   *  by the user.
   **/
  private void writeCodonUsage()
  {
    if(!checkForSelectionFeatures())
      return;

    final File write_file =
      getWriteFile("Select a codon usage output file ...", "usage", null);

    if(write_file == null)
      return;

    try 
    {
      final PrintWriter writer = new PrintWriter(new FileWriter(write_file));

      final int[][][] codon_counts = new int[4][4][4];

      final FeatureVector features_to_write =
        getSelection().getAllFeatures();

      int codon_total = 0;

      for(int i = 0; i < features_to_write.size(); ++i) 
      {
        final Feature selection_feature = features_to_write.elementAt(i);

        for(int base1 = 0 ; base1 < 4 ; ++base1) 
        {
          for(int base2 = 0 ; base2 < 4 ; ++base2) 
          {
            for(int base3 = 0 ; base3 < 4 ; ++base3) 
            {
              codon_counts[base1][base2][base3] +=
                selection_feature.getCodonCount(base1, base2, base3);
            }
          }
        }

        codon_total += selection_feature.getTranslationBasesLength() / 3;
      }

      for(int base1 = 0 ; base1 < 4 ; ++base1) 
      {
        for(int base3 = 0 ; base3 < 4 ; ++base3) 
        {
          final StringBuffer buffer = new StringBuffer();

          for(int base2 = 0 ; base2 < 4 ; ++base2) 
          {
            buffer.append(Bases.letter_index[base1]);
            buffer.append(Bases.letter_index[base2]);
            buffer.append(Bases.letter_index[base3]);
            buffer.append(' ');

            final float per_thousand;

            if(codon_total > 0)
            {
              per_thousand =
                10000 * codon_counts[base1][base2][base3] / codon_total *
                1.0F / 10;
            } 
            else
              per_thousand = 0.0F;

            buffer.append(per_thousand);
            buffer.append("( ").append(codon_counts[base1][base2][base3]);
            buffer.append(")  ");
          }

          writer.println(buffer.toString());
        }
      }

      writer.close();
    } 
    catch(IOException e)
    {
      new MessageDialog(getParentFrame(),
                        "error while writing: " + e.getMessage());
    }
  }

  /**
   *  The EntryGroup object that was passed to the constructor.
   **/
  private EntryGroup entry_group;
}
