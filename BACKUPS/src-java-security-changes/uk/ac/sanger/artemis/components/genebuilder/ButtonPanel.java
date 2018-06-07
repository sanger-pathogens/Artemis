package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.Iterator;
import java.util.Set;

import javax.swing.JButton;
import javax.swing.JOptionPane;
import javax.swing.JToolBar;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.TransferAnnotationTool;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;


public class ButtonPanel extends JToolBar
{
  private static final long serialVersionUID = 1L;
  private BasicGeneBuilderFrame gbFrame;
  
  public ButtonPanel(final BasicGeneBuilderFrame gbFrame, 
                     final EntryGroup entry_group)
  {
    super();
    
    this.gbFrame = gbFrame;
    
    final JButton complement_button = new JButton("Complement");
    add(complement_button);
    complement_button.addActionListener(new ActionListener () 
    {
      public void actionPerformed(ActionEvent e) 
      {
        complementLocation();
      }
    });
    
    if(GeneUtils.isDatabaseEntry(getFeature().getEmblFeature()))
    {
      final JButton refresh_button = new JButton("Refresh");
      add(refresh_button);
      refresh_button.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          refresh();
        }
      });
    }

    final JButton grab_button = new JButton("Grab Range");
    add(grab_button);
    grab_button.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        grabSelectedRange();
      }
    });

    final JButton remove_button = new JButton("Remove Range");
    add(remove_button);
    remove_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e)
      {
        removeSelectedRange();
      }
    });

    final JButton select_button = new JButton("Select Feature");
    add(select_button);
    select_button.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e) 
      {
        gbFrame.getSelection().set(getFeature());
      }
    });
    
    final JButton transferAnnotationBbutton = new JButton("TAT");
    add(transferAnnotationBbutton);
    transferAnnotationBbutton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e) 
      {      
        new TransferAnnotationTool(getFeature(), entry_group, gbFrame.getMatchForm());
      }
    });

  }
  
  private Feature getFeature()
  {
    return gbFrame.getFeature();
  }
  
  /**
   *   Complement the current location_text.
   **/
  private void complementLocation() 
  {
    if(GeneUtils.isDatabaseEntry(getFeature().getEmblFeature()))
    {
      final ChadoCanonicalGene chadoGene = 
        ((GFFStreamFeature)getFeature().getEmblFeature()).getChadoGene();
      GeneUtils.complementGeneModel(chadoGene);
      gbFrame.dispose(true);
      return;  
    }
    
    if(rationalizeLocation()) 
    {
      if(gbFrame.getLocationText().getText().startsWith("complement(")) 
      {
        final String new_text = gbFrame.getLocationText().getText().substring(11);
        if (new_text.endsWith(")")) 
        {
          final String new_text2 =
            new_text.substring(0, new_text.length () - 1);
          gbFrame.getLocationText().setText(new_text2);
        } 
        else 
          gbFrame.getLocationText().setText(new_text);
      }
      else 
      {
        final String new_text = gbFrame.getLocationText().getText ();
        gbFrame.getLocationText().setText("complement(" + new_text + ")");
      }
    } 
    else 
      new MessageDialog(null, "complement failed - " +
                        "current location cannot be parsed");
  }
  
  /**
   *  Attempt to parse the current location_text as a Location.  If it can be
   *  parsed it will be canonicalized (ie. the complement, if any, will be
   *  outermost).  Returns true if and only if the location_text could be
   *  parsed.
   **/
  private boolean rationalizeLocation()
  {
    try 
    {
      final Location location = new Location(gbFrame.getLocationText().getText());
      gbFrame.getLocationText().setText(location.toStringShort());
      return true;
    }
    catch(LocationParseException e) 
    {
      return false;
    }
  }
  
  /**
   * Refresh the annotation for the active feature
   */
  private void refresh()
  {
    Feature edit_feature = gbFrame.getFeature();
    // refresh from the database
    final DatabaseDocument originalDocument =
      (DatabaseDocument)((DocumentEntry)edit_feature.getEmblFeature().getEntry()).getDocument();

    final Set uniquenames = 
      ((GFFStreamFeature)edit_feature.getEmblFeature()).getSegmentRangeStore().keySet();
    final Iterator it = uniquenames.iterator();
    final String uniquename = (String)it.next();
    final DatabaseDocument newDocument = new DatabaseDocument(originalDocument,
        uniquename, null, true, null);
    newDocument.setLazyFeatureLoad(false);
    newDocument.setReadChildren(false);
    
    try
    {
      DatabaseDocumentEntry dbentry = new DatabaseDocumentEntry(newDocument, null);
      uk.ac.sanger.artemis.io.Feature databaseFeature = dbentry.getAllFeatures().featureAt(0);
      
      // compare timelastmodified
      Qualifier qualifier = edit_feature.getQualifierByName("timelastmodified");
      String active_timelastmodified = (String)qualifier.getValues().get(0);
      qualifier = databaseFeature.getQualifierByName("timelastmodified");
      String database_timelastmodified = (String)qualifier.getValues().get(0);
      
      if(active_timelastmodified.equals(database_timelastmodified))
      {
        JOptionPane.showMessageDialog(this, 
            "No new changes found for the feature\n"+
            uniquename+"\n"+
            "in the database since:\n"+database_timelastmodified, 
            "No Updates", JOptionPane.INFORMATION_MESSAGE);
        return;
      }
      else
      {
        JOptionPane.showMessageDialog(this, 
            "Changes found for the feature\n"+
            uniquename+"\n"+
            "in the database at:\n"+database_timelastmodified, 
            "Changes Found", JOptionPane.INFORMATION_MESSAGE);
      }
      
      final QualifierVector db_qv = databaseFeature.getQualifiers();
      final QualifierVector qv = edit_feature.getQualifiers();
      final QualifierVector new_qv = new QualifierVector();
      
      for(int i=0; i<qv.size(); i++)
      {
        Qualifier q = (Qualifier)qv.get(i);
        if(q.getName().equals("Parent") ||
           q.getName().equals("Derives_from") ||
           q.getName().equals("ID")  )
          new_qv.addQualifierValues(q);
      }
      
      
      for(int i=0; i<db_qv.size(); i++)
      {
        Qualifier q = (Qualifier)db_qv.get(i);
        if(q.getName().equals("Parent") ||
           q.getName().equals("Derives_from") ||
           q.getName().equals("ID")  )
          continue;
        new_qv.add(q);  
      }
      
      edit_feature.getQualifiers().removeAllElements();
      edit_feature.getQualifiers().addAll(new_qv);
      gbFrame.setQualifiers();
    }
    catch(EntryInformationException e) {}
    catch(IOException e) {}    
  }
  
  
  /**
   *  Add the currently selected range to location_text.
   **/
  private void grabSelectedRange() 
  {
    if(!rationalizeLocation ()) 
    {
      new MessageDialog(null,
                        "grab failed - current location cannot be parsed");
      return;
    }

    final Range selected_range = gbFrame.getSelection().getSelectionRange();

    if(selected_range == null) 
    {
      new MessageDialog(null, "grab failed - nothing is selected");
      return;
    }

    // save it in case it gets mangled
    final String old_location_text = gbFrame.getLocationText().getText();

    if (old_location_text.endsWith ("))")) 
    {
      final String new_text =
        old_location_text.substring(0, old_location_text.length () - 2);

      gbFrame.getLocationText().setText(new_text + "," + selected_range.getStart () +
                            ".." + selected_range.getEnd () + "))");
    } 
    else
    {
      if(old_location_text.endsWith (")")) 
      {
        final String new_text =
          old_location_text.substring(0, old_location_text.length () - 1);

        gbFrame.getLocationText().setText(new_text + "," + selected_range.getStart () +
                              ".." + selected_range.getEnd () + ")");
      } 
      else
        gbFrame.getLocationText().setText(old_location_text + "," +
                              selected_range.getStart () +
                              ".." + selected_range.getEnd ());
    }

    if(!rationalizeLocation())
    {
      gbFrame.getLocationText().setText (old_location_text);
      new MessageDialog(null, "grab failed - location cannot be parsed after " +
                        "grabbing");
    }
  }

  /**
   *  Remove the currently selected range of bases from location_text.
   **/
  private void removeSelectedRange()
  {
    if(!rationalizeLocation())
    {
      new MessageDialog(null,
                        "grab failed - current location cannot be parsed");
      return;
    }

    final MarkerRange selected_marker_range =
                                       gbFrame.getSelection().getMarkerRange();

    if(selected_marker_range == null) 
    {
      new MessageDialog(null, "remove range failed - no bases are selected");
      return;
    }

    final Range selected_range = selected_marker_range.getRawRange();

    if (selected_marker_range.getStrand() != getFeature().getStrand()) 
    {
      new MessageDialog(null, "remove range failed - you need to select " +
                         "some bases on the other strand");
      return;
    }

    final Location location;

    try 
    {
      location = new Location(gbFrame.getLocationText().getText ());
    }
    catch (LocationParseException e) 
    {
      // this shouldn't happen because we called rationalizeLocation ()
      throw new Error("internal error - unexpected exception: " + e);
    }

    final Range location_total_range = location.getTotalRange();

    if(!selected_range.overlaps(location_total_range))
    {
      new MessageDialog(null, "remove range failed - the range you " +
                        "selected does not overlap the feature");
      return;
    }

    if(selected_range.contains(location_total_range)) 
    {
      new MessageDialog(null, "remove range failed - the range you " +
                        "selected overlaps the whole feature");
      return;
    }

    final RangeVector location_ranges = location.getRanges();
    final boolean location_is_complemented = location.isComplement();

    final RangeVector new_ranges = new RangeVector();

    // if the selected_range completely covers a range remove the
    // range. otherwise if the selected_range is completely within one of the
    // ranges two new ranges are created.  if the selected_range is not
    // completely contained then the appropriate end of the range is truncated
    for(int i = 0; i < location_ranges.size(); ++i) 
    {
      final Range this_range = (Range)location_ranges.elementAt(i);

      if(selected_range.overlaps(this_range)) 
      {
        try 
        {
          if(this_range.contains(selected_range) &&
             this_range.getStart() != selected_range.getStart() &&
             this_range.getEnd() != selected_range.getEnd()) 
          {
            // chop a piece out of the middle and make two new ranges
            final Range new_start_range =
                    this_range.change(selected_range.getEnd() + 1,
                                      this_range.getEnd());
            new_ranges.add(new_start_range);
            final Range new_end_range =
              this_range.change(this_range.getStart(),
                                selected_range.getStart() - 1);
            new_ranges.add(new_end_range);
          } 
          else
          {
            if(selected_range.contains(this_range)) {
              // delete (ie. don't copy) the range
            } 
            else
            {
              if(this_range.getStart() < selected_range.getStart()) 
              {
                // truncate the end of the range
                final Range new_start_range =
                  this_range.change(this_range.getStart(),
                                    selected_range.getStart() - 1);
                new_ranges.add(new_start_range);
              } 
              else
              {
                if(this_range.getEnd() > selected_range.getEnd())
                {
                  // truncate the start of the range
                  final Range new_end_range =
                    this_range.change(selected_range.getEnd() + 1,
                                      this_range.getEnd());
                  new_ranges.add(new_end_range);
                } 
                else
                  throw new Error ("internal error - can't remove range");
              }
            }
          }
        }
        catch(OutOfRangeException e)
        {
          throw new Error ("internal error - unexpected exception: " + e);
        }
      }
      else  
        new_ranges.add(this_range); // copy it unchanged
    }

    final Location new_location =
      new Location(new_ranges, location_is_complemented);

    gbFrame.getLocationText().setText(new_location.toStringShort());
  }
  
}