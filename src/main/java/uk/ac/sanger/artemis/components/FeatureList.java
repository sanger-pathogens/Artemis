/* FeatureList.java
 *
 * created: Fri Oct  9 1998
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998,1999,2000,2001,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeatureList.java,v 1.28 2008-10-23 14:44:28 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.plot.CodonUsageAlgorithm;

import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.StreamQualifier;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

import java.awt.event.MouseEvent;
import java.awt.event.InputEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.Container;
import java.awt.Color;
import java.awt.Point;
import java.awt.Graphics;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.text.NumberFormat;

import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.JViewport;
import javax.swing.JComponent;

/**
 *  This component gives the user a list containing the details the current
 *  Features.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureList.java,v 1.28 2008-10-23 14:44:28 tjc Exp $
 *
 **/

public class FeatureList extends EntryGroupPanel
  implements EntryGroupChangeListener,
             EntryChangeListener, FeatureChangeListener,
             SelectionChangeListener, DisplayComponent
{
  private static final long serialVersionUID = 1L;

  /** true if correlation scores should be shown */
  private boolean show_correlation_scores = false;

  /**
   *  This is set to true by selectionChanged() and used by paintComponent().
   **/
  private boolean selection_changed_flag = false;

  /** colour used to draw the background. */
  private Color background_colour = Color.white;

  /**
   *  If true this component will show Feature.getIDString() (ie /gene or
   *  /label) instead of the key.
   **/
  private boolean show_gene_names = false;

  private String user_defined_qualifier = null;

  /** show Feature.getSystematicName() */
  private boolean show_systematic_names = false;
   
  /** show the /product qualifier instead of /note field. */
  private boolean show_products = false;

  /** show all the qualifiers after the note. */
  private boolean show_qualifiers = false;

  /**
   *  The is the maximum width of the strings containing the feature start and
   *  stop positions.  Set in the constructor.
   **/
  private int max_base_pos_width;

  /** JScrollPane viewport that this panel is the view of */
  private JViewport viewport = null;
  
  /**
   *  Create a new FeatureList with the default number of rows.
   *  @param entry_group The EntryGroup that this component will display.
   *  @param selection The Selection object for this component.  Selected
   *    objects will be highlighted.
   *  @param goto_event_source The object to use when we need to call
   *    gotoBase().
   **/
  public FeatureList(final EntryGroup entry_group,
                     final Selection selection,
                     final GotoEventSource goto_event_source,
                     final BasePlotGroup base_plot_group)
  {
    super(entry_group, selection, goto_event_source, base_plot_group);

    addMouseListener(new MouseAdapter() 
    {
      private FeaturePopup popup = null;

      /**
       *  Listen for mouse press events so that we can do popup menus and
       *  selection.
       **/
      public void mousePressed(MouseEvent event)
      {
        if(isMenuTrigger(event)) 
        {
          if(popup == null)
            popup = new FeaturePopup(FeatureList.this,
                                     getEntryGroup(),
                                     getSelection(),
                                     getGotoEventSource(),
                                     getBasePlotGroup());
          final JComponent parent = (JComponent)event.getSource();

          popup.show(parent, event.getX(), event.getY());
        } 
        else 
          handleCanvasMousePress(event);
      }
    });

    getSelection().addSelectionChangeListener(this);

    // changes to the EntryGroup will be noticed by listening for EntryChange
    // and FeatureChange events.

    getEntryGroup().addEntryGroupChangeListener(this);
    getEntryGroup().addEntryChangeListener(this);
    getEntryGroup().addFeatureChangeListener(this);

    // find the maximum posible width for the high and low positions
    final int sequence_length = getEntryGroup().getSequenceLength();
    max_base_pos_width = (int)(Math.log(sequence_length)/Math.log(10)) + 1;

    if(max_base_pos_width < 4) 
      max_base_pos_width = 4;

    setBackground(background_colour);
  }

  /**
   *  Remove this component from all the listener lists it is on.
   **/
  void stopListening()
  {
    getSelection().removeSelectionChangeListener(this);
    getEntryGroup().removeEntryGroupChangeListener(this);
    getEntryGroup().removeEntryChangeListener(this);
    getEntryGroup().removeFeatureChangeListener(this);
  }

  /**
   *  Set value of the show correlation scores flag.
   *  @param show_correlation_scores Show correlation scores in the list if
   *    and only if this argument is true.
   **/
  protected void setCorrelationScores(final boolean show_correlation_scores) 
  {
    if(this.show_correlation_scores != show_correlation_scores) 
    {
      this.show_correlation_scores = show_correlation_scores;
      repaint();
    } 
  }

  /**
   *  Get the value of the "show correlation scores" flag.
   **/
  protected boolean getCorrelationScores() 
  {
    return show_correlation_scores;
  }

  /**
   *  Set value of the show /gene flag.
   *  @param show_gene_names If true this component will show the /gene (really
   *    Feature.getIDString()) instead of the key.
   **/
  protected void setShowGenes(final boolean show_gene_names) 
  {
    if(this.show_gene_names != show_gene_names) 
    {
      this.show_gene_names = show_gene_names;
      repaint();
    }
  }


  /**
   *  Set value of the show /systematic_id flag.
   *  @param show_systematic_names If true this component will show the /gene (really
   *    Feature.getSystematicName()) instead of the key.
   **/
  protected void setShowSystematicID(final boolean show_systematic_names)
  {
    if(this.show_systematic_names != show_systematic_names)
    {
      this.show_systematic_names = show_systematic_names;
      repaint();
    }
  }


  protected void setShowUserDefinedQualifier(final String user_defined_qualifier)
  {
    this.user_defined_qualifier = user_defined_qualifier;
    repaint();
  }
  
  protected StringVector getShowUserDefinedQualifier()
  {
    if(user_defined_qualifier == null)
      return null;
    return StringVector.getStrings(user_defined_qualifier);
  }

  /**
   *  Get the value of the "show genes" flag.
   **/
  protected boolean getShowGenes() 
  {
    return show_gene_names;
  }


  /**
   *  Get the value of the "show systematic id" flag.
   **/
  protected boolean getShowSysID()
  {
    return show_systematic_names;
  }


  /**
   *  Set value of the show qualifiers flag.
   *  @param show_quailfiers If true this component will show all the
   *    qualifiers after the note.
   **/
  protected void setShowQualifiers(final boolean show_qualifiers) 
  {
    if(this.show_qualifiers != show_qualifiers) 
    {
      if(show_qualifiers)
        user_defined_qualifier = null;

      this.show_qualifiers = show_qualifiers;
      repaint();
    }
  }

  /**
   *  Get the value of the "show qualifiers" flag.
   **/
  protected boolean getShowQualifiers() 
  {
    return show_qualifiers;
  }

  /**
   *  Set value of the show /product flag.
   *  @param show_products If true this component will show the /product
   *    qualifier instead of the /note.
   **/
  protected void setShowProducts(final boolean show_products) 
  {
    if(this.show_products != show_products) 
    {
      if(show_products)
        user_defined_qualifier = null;

      this.show_products = show_products;
      repaint();
    }
  }

  /**
   *  Get the value of the "show products" flag.
   **/
  protected boolean getShowProducts() 
  {
    return show_products;
  }
  

  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can update the display if entries
   *  are added or deleted.
   **/
  public void entryGroupChanged(EntryGroupChangeEvent event) 
  {
    final int hgt = getEntryGroup().getAllFeaturesCount() *
                               getLineHeight();

    setPreferredSize(new Dimension(getSize().width*4,hgt));
    revalidate();
    repaint();
  }

  /**
   *  Implementation of the FeatureChangeListener interface.
   **/
  public void featureChanged(FeatureChangeEvent event) 
  {
    if(!isVisible()) 
      return;

    repaint();
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so that we can update the list if features are added
   *  or deleted.
   **/
  public void entryChanged(EntryChangeEvent event) 
  {
    if(!isVisible()) 
      return;

    repaint();
  }

  /**
   *  Implementation of the SelectionChangeListener interface.  We listen to
   *  SelectionChange events so that we can update the list to reflect the
   *  current selection.
   **/
  public void selectionChanged(SelectionChangeEvent event) 
  {
    if(!isVisible())
      return;

    // don't bother with events we sent ourself
    if(event.getSource() == this) 
      return;

    // if the selected range changes we don't care
    if(getSelection().getMarkerRange() != null &&
       event.getType() == SelectionChangeEvent.OBJECT_CHANGED) 
      return;

    selection_changed_flag = true;

    onSelectionChange();
    repaint();
  }

  /**
   *  Return a vector containing the text that is shown in the list - one
   *  String per line.
   **/
  protected StringVector getListStrings() 
  {
    final StringVector return_vector = new StringVector();
    final FeatureEnumeration test_enumerator = getEntryGroup().features();

    while(test_enumerator.hasMoreFeatures()) 
    {
      final Feature this_feature = test_enumerator.nextFeature();
      return_vector.add(makeFeatureString(this_feature, true));
    }

    return return_vector;
  }


  /**
  * Return the JViewport that this component is contained in.
  */
  protected JViewport getViewport()
  {
    if(viewport != null)
      return viewport;

    Container container = getParent();
    while(!(container instanceof JScrollPane))
      container = container.getParent();

    viewport = ((JScrollPane)container).getViewport();
    return viewport;
  }

  /**
  *
  * Find the point at the top right hand corner of the
  * scroll pane.
  *
  */
  private Point getScrollPoint()
  {
    return getViewport().getViewPosition();
  }

  /**
   *  Handle a mouse press event on the drawing canvas - select on click,
   *  select and broadcast it on double click.
   **/
  private void handleCanvasMousePress(final MouseEvent event)
  {
    if(event.getID() != MouseEvent.MOUSE_PRESSED) 
      return;

    requestFocus();

    if(!event.isShiftDown()) 
      getSelection().clear();

    final int clicked_feature_index = event.getY()/getLineHeight();

    if(clicked_feature_index < getEntryGroup().getAllFeaturesCount())
    {
      final FeatureVector selected_features =
        getSelection().getAllFeatures();

      final Feature clicked_feature =
        getEntryGroup().featureAt(clicked_feature_index);

      if(selected_features.contains(clicked_feature)) 
      {
        getSelection().remove(clicked_feature);
        getSelection().removeSegmentsOf(clicked_feature);
      } 
      else 
        getSelection().add(clicked_feature);

      if(event.getClickCount() == 2) 
      {
        makeSelectionVisible();

        if((event.getModifiers() & InputEvent.BUTTON2_MASK) != 0 ||
            event.isAltDown()) 
        {
          if(Options.readWritePossible()) 
          {
            final JFrame frame = new JFrame("Artemis Feature Edit: " + 
                clicked_feature.getIDString() +
                (clicked_feature.isReadOnly() ?
                    "  -  (read only)" :
                    ""));
            
            final FeatureEdit fe = new FeatureEdit(clicked_feature, getEntryGroup(),
                                       getSelection(), getGotoEventSource(), frame);
            frame.addWindowListener(new WindowAdapter() 
            {
              public void windowClosing(WindowEvent event) 
              {
                fe.stopListening();
                frame.dispose();
              }
            });
            
            frame.getContentPane().add(fe);
            frame.pack();

            final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
            frame.setLocation(new Point((screen.width - getSize().width)/2,
                                        (screen.height - getSize().height)/2));
            frame.setVisible(true);
          }
        }
      }

    }
  }

  private void onSelectionChange()
  {
    if(!selection_changed_flag)
      return;

    selection_changed_flag = false;
    final FeatureVector selected_features =
                         getSelection().getAllFeatures();

    if(selected_features.size() > 0)
    {
      Point viewPoint = getScrollPoint();
      final EntryGroup entry_group = getEntryGroup();
      final int feature_count = entry_group.getAllFeaturesCount();

      // set to true if any of the selected features is visible
      boolean a_selected_feature_is_visible = false;

      int first_line_in_view = viewPoint.y/getLineHeight();

      if(first_line_in_view == -1)
        first_line_in_view = 0;

      int numberLines = linesInView();
      for(int i = first_line_in_view;
          i < feature_count && i <= first_line_in_view + numberLines;
          ++i)
      {
        final Feature this_feature = entry_group.featureAt(i);
        if(selected_features.contains(this_feature))
        {
          a_selected_feature_is_visible = true;
          break;
        }
      }

      if(!a_selected_feature_is_visible)
      {
        // make the first selected feature visible
        final Feature first_selected_feature =
          selected_features.elementAt(0);

        final int index_of_first_selected_feature =
                         entry_group.indexOf(first_selected_feature);

        if( index_of_first_selected_feature > -1 &&
           (index_of_first_selected_feature < first_line_in_view ||
            index_of_first_selected_feature >= first_line_in_view + numberLines))
        {
          getViewport().setViewPosition(new Point(0,
                                index_of_first_selected_feature * getLineHeight()));
        }
      }
    }
  }

  /**
   *  The main paint function for the canvas.  An off screen image used for
   *  double buffering when drawing the canvas.
   *  @param g The Graphics object of the canvas.
   **/
  protected void paintComponent(Graphics g) 
  {
    super.paintComponent(g);

    if(!isVisible()) 
      return;

    Point viewPoint = getScrollPoint();
    final int feature_count = getEntryGroup().getAllFeaturesCount();

    // must check size in case new features/entries added
    if(feature_count*getLineHeight() > getPreferredSize().height)
    {
      final int hgt = feature_count * getLineHeight();
      setPreferredSize(new Dimension(getSize().width*4,hgt));
      revalidate();
    }

    if(feature_count != 0) 
    {
      final int lines_in_view = linesInView()+1;
      int first_index_in_view = (viewPoint.y/getLineHeight());

      if(first_index_in_view == -1) 
        first_index_in_view = 0;

      int last_index_in_view;

      if(lines_in_view < feature_count - first_index_in_view) 
        last_index_in_view = first_index_in_view + lines_in_view;
      else 
        last_index_in_view = feature_count - 1;

      final FeatureVector features_in_view =
        getEntryGroup().getFeaturesInIndexRange(first_index_in_view,
                                                last_index_in_view);

      g.setFont(getFont());

      final int features_in_view_size = features_in_view.size();
      for(int i = 0; i < features_in_view_size; i++)
      {
        final Feature this_feature  = features_in_view.elementAt(i);
        final String feature_string = makeFeatureString(this_feature, false);
        drawFeatureLine(g, this_feature, feature_string);
      }
    }
  }


  /**
   *  Return the number of visible text lines on canvas.
   **/
  private int linesInView() 
  {
    return getViewport().getExtentSize().height/getLineHeight();
  }

  /**
   *  Draw the given Feature at the given line of the list, taking the
   *  selection into account.
   **/
  private void drawFeatureLine(final Graphics g,
                               final Feature feature,
                               final String feature_string)
  {
    // width of coloured blob at the left of the text
    final int BOX_WIDTH = getLineHeight();
    final int y_pos = getEntryGroup().indexOf(feature)*BOX_WIDTH;

    final Color feature_colour = feature.getColour();

    // default colour is white
    if(feature_colour == null) 
      g.setColor(Color.white);
    else 
      g.setColor(feature_colour);
    
    g.fillRect(1, y_pos+1,
               BOX_WIDTH, BOX_WIDTH - 1);

    g.setColor(Color.black);
    
    final FeatureVector selected_features = getSelection().getAllFeatures();
    if(selected_features.contains(feature)) 
    {
      // draw in reverse
      g.fillRect(BOX_WIDTH + 4, y_pos,
                 getSize().width + getScrollPoint().x,
                 getLineHeight());
      g.setColor(background_colour);
    } 

    if( feature.getEmblFeature() instanceof GFFStreamFeature &&
        !getSelection().contains(feature) &&
        !((GFFStreamFeature)feature.getEmblFeature()).isVisible() )
    {
      //
      // use gray for the key if the feature is NOT visible
      g.setColor(Color.gray);
      int ind = feature_string.indexOf(' ');
      final String keyString = feature_string.substring(0, ind);
      g.drawString(keyString,
          BOX_WIDTH + 5,
          y_pos + getFontAscent());
      
      g.setColor(Color.black);
      g.drawString(feature_string.substring(ind),
          BOX_WIDTH + 5 + getFontMetrics(getFont()).stringWidth(keyString),
          y_pos + getFontAscent());
    }
    else
      g.drawString(feature_string,
                 BOX_WIDTH + 5,
                 y_pos + getFontAscent());
  }

  /**
   *  Return a String object suitable for displaying in the list of features.
   *  @param dont_truncate if true the gene name / key field won't be
   *    truncated if it is longer than the field width
   **/
  private String makeFeatureString(final Feature feature,
                                   final boolean dont_truncate) 
  {
    String key_string;
    final int KEY_FIELD_WIDTH = 15;

    if(show_gene_names)
    {
      key_string = feature.getGeneName();

      if(key_string == null)
        key_string = feature.getSystematicName();

      if(key_string.length() > KEY_FIELD_WIDTH && !dont_truncate) 
        key_string = key_string.substring(0, KEY_FIELD_WIDTH);
    } 
    else if(show_systematic_names)
    {
      key_string = feature.getSystematicName();

      if(key_string.length() > KEY_FIELD_WIDTH && !dont_truncate)
        key_string = key_string.substring(0, KEY_FIELD_WIDTH);
    }
    else 
      key_string = feature.getKey().toString();

    final Marker low_marker  = feature.getFirstBaseMarker();
    final Marker high_marker = feature.getLastBaseMarker();

    final StringBuffer description_string_buffer = new StringBuffer();

    if(user_defined_qualifier != null && !user_defined_qualifier.equals(""))
    {
      try
      { 
        final StringVector sv = StringVector.getStrings(user_defined_qualifier);
        for(int i=0; i<sv.size(); i++)
        {
          final Qualifier q = feature.getQualifierByName((String)sv.get(i));
          if(q != null)
          {
            final StringVector values = q.getValues();
            if(values != null)
            {
              if(values.size() == 0)
                description_string_buffer.append("/"+sv.get(i)+" ");
              else
                for(int j=0; j<values.size(); j++)  // show multiple values
                  description_string_buffer.append("/"+sv.get(i)+
                    (values.get(j) == null ? "" : "="+values.get(j))+" ");
            }
            else
              description_string_buffer.append("/"+sv.get(i)+" ");
          }
        }
      }
      catch(InvalidRelationException ire){}
    }
    else if(show_products) 
    {
      final String product_string = feature.getProductString();

      if(product_string == null) 
      {
        // description is not blank
        if(feature.isCDS())
          description_string_buffer.append("[no /product]");
      } 
      else 
        description_string_buffer.append(product_string);
    }
    else 
    {
      String note = null;
      if(feature.getEmblFeature() instanceof GFFStreamFeature)
      {
        try
        {
          if(feature.getValueOfQualifier("isObsolete") != null &&
             feature.getValueOfQualifier("isObsolete").equals("true"))
            note = "obsolete";
          else
            note = feature.getValueOfQualifier("comment");
        }
        catch(InvalidRelationException e){}
      }
      
      if(note == null)
        note = feature.getNote();

      if(note != null && note.length() != 0) 
      {
        final int QUALIFIER_COLUMN = 10;

        final String note_string =
          padRightWithSpaces(note, QUALIFIER_COLUMN);

        description_string_buffer.append(note_string);
        description_string_buffer.append("   ");
      }

      if(show_qualifiers) 
        description_string_buffer.append(getQualifierString(feature));
    }


    String low_pos;
    String high_pos;
    if(low_marker == null || high_marker == null) 
    {
      low_pos  = "unknown";
      high_pos = "unknown";
    }
    else 
    {
      if(low_marker.getRawPosition() < high_marker.getRawPosition()) 
      {
        low_pos = String.valueOf(low_marker.getRawPosition());
        high_pos = String.valueOf(high_marker.getRawPosition());
      } 
      else
      {
        low_pos = String.valueOf(high_marker.getRawPosition());
        high_pos = String.valueOf(low_marker.getRawPosition());
      }
    }

    if(feature.getEmblFeature() instanceof GFFStreamFeature)
    {
      try
      {
        if(feature.getQualifierByName("Start_range") != null)
          low_pos = "<"+low_pos;
        if(feature.getQualifierByName("End_range") != null)
          high_pos = ">"+high_pos;
      }
      catch (InvalidRelationException e){}
    }
    else
    {
      if(feature.getLocation().isPartial(true)) // 5prime
      {
        if(feature.isForwardFeature()) 
          low_pos = "<"+low_pos;
        else
          high_pos = ">"+high_pos;
      }
      if(feature.getLocation().isPartial(false)) // 3prime
      {
        if(feature.isForwardFeature())
          high_pos = ">"+high_pos;
        else
          low_pos = "<"+low_pos;
      }
    }

    StringBuffer new_list_line = new StringBuffer();

    new_list_line.append(padRightWithSpaces(key_string, KEY_FIELD_WIDTH));
    new_list_line.append(" ");

    new_list_line.append(padLeftWithSpaces(low_pos, max_base_pos_width));
    new_list_line.append(" ");
    new_list_line.append(padLeftWithSpaces(high_pos, max_base_pos_width));
    new_list_line.append(" ");

    if(feature.isForwardFeature()) 
      new_list_line.append("   ");
    else
      new_list_line.append("c  ");

    if(show_correlation_scores)
    {
      if(feature.isCDS() || 
         feature.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL)) 
      {
        new_list_line.append(getScoresString(feature));
        new_list_line.append("  ");
      } 
      else 
      {
        new_list_line.append("                         ");
        if(getBasePlotGroup().getCodonUsageAlgorithm() != null) 
          new_list_line.append("      ");
      }
    }

    new_list_line.append(description_string_buffer.toString());

    return new_list_line.toString();
  }


  /**
   *  Return a String containing the given Qualifier and it's values (in EMBL
   *  format).
   *  @param start_index ignore the values before this index
   **/
  private String formatQualifier(final String qualifier_name,
                                 final Feature feature,
                                 final int start_index) 
  {
    final StringBuffer buffer = new StringBuffer();

    try 
    {
      final Qualifier qualifier = feature.getQualifierByName(qualifier_name);

      if(qualifier != null) 
      {
        final EntryInformation entry_information =
          feature.getEntry().getEntryInformation();

        final QualifierInfo qualifier_info =
          entry_information.getQualifierInfo(qualifier_name);

        final StringVector qualifier_strings =
          StreamQualifier.toStringVector(qualifier_info,
                                         qualifier);

        final int qualifier_strings_size = qualifier_strings.size();
        for(int i = start_index; i < qualifier_strings_size; ++i)
        {
          final String qualifier_string = (String)qualifier_strings.elementAt(i);
          buffer.append(qualifier_string + " ");
        }
      }
    } 
    catch(InvalidRelationException e) {}

    return buffer.toString();
  }

  /**
   *  Return a String containing all the qualifiers of the given Feature
   *  (except /note) in EMBL format.  Any /similarity qualifier will come
   *  first.
   **/
  private String getQualifierString(final Feature feature) 
  {
    final StringBuffer buffer = new StringBuffer();

    final QualifierVector qualifiers = feature.getQualifiers();

    // if there is a /note and it has more than one value put it next (without
    // the first value)
    final Qualifier note_qualifier =
      qualifiers.getQualifierByName("note");

    if(note_qualifier != null && note_qualifier.getValues().size() > 1) 
    {
      buffer.append(formatQualifier("note", feature, 1));
      buffer.append(" ");
    }

    // put /similarity before all but the /note qualifier
    final Qualifier similarity_qualifier =
      qualifiers.getQualifierByName("similarity");

    if(similarity_qualifier != null) 
    {
      buffer.append(formatQualifier("similarity", feature, 0));
      buffer.append(" ");
    }

    final int qualifiers_size = qualifiers.size();
    for(int i = 0 ; i < qualifiers_size; ++i) 
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(i);
      final String this_qualifier_name = this_qualifier.getName();

      if(!this_qualifier_name.equals("note") &&
         !this_qualifier_name.equals("similarity")) 
      {
        buffer.append(formatQualifier(this_qualifier_name, feature, 0));
        buffer.append(" ");
      }
    }

    return buffer.toString();
  }

  /**
   *  Return a String containing the correlation scores.
   **/
  private String getScoresString(final Feature feature)
  {
    //final int base_total = feature.getTranslationBases().length();
    final int c_total = feature.getBaseCount(Bases.getIndexOfBase('c'));
    final int g_total = feature.getBaseCount(Bases.getIndexOfBase('g'));

    final int g1_count =
      feature.getPositionalBaseCount(0, Bases.getIndexOfBase('g'));

    final int c3_count =
      feature.getPositionalBaseCount(2, Bases.getIndexOfBase('c'));
    final int g3_count =
      feature.getPositionalBaseCount(2, Bases.getIndexOfBase('g'));

    final double c3_score = 100.0 * (3 * c3_count - c_total) / c_total;
    final double g1_score = 100.0 * (3 * g1_count - g_total) / g_total;
    final double g3_score = 100.0 * (3 * g3_count - g_total) / g_total;

    final double cor1_2_score = feature.get12CorrelationScore();

    final NumberFormat number_format = NumberFormat.getNumberInstance();

    number_format.setMaximumFractionDigits(1);
    number_format.setMinimumFractionDigits(1);

    final String cor1_2_score_string = number_format.format(cor1_2_score);
    final String c3_score_string;
    final String g1_score_string;
    final String g3_score_string;

    if(c_total == 0) 
      c3_score_string = "ALL";
    else 
      c3_score_string = number_format.format(c3_score);

    if(g_total == 0) 
      g1_score_string = "ALL";
    else 
      g1_score_string = number_format.format(g1_score);

    if(g_total == 0)
      g3_score_string = "ALL";
    else
      g3_score_string = number_format.format(g3_score);

    String codon_usage_score_string = "";

    final CodonUsageAlgorithm codon_usage_alg =
      getBasePlotGroup().getCodonUsageAlgorithm();

    if(codon_usage_alg != null) 
    {
      number_format.setMaximumFractionDigits(3);
      number_format.setMinimumFractionDigits(3);

      codon_usage_score_string =
        number_format.format(codon_usage_alg.getFeatureScore(feature)) + " ";
    }

    return codon_usage_score_string +
           padRightWithSpaces(cor1_2_score_string, 5) + " " +
           padRightWithSpaces(c3_score_string, 5) + " " +
           padRightWithSpaces(g1_score_string, 5) + " " +
           padRightWithSpaces(g3_score_string, 5);
  }

  /**
   *  Return the given string padded with spaces to the given width.  The
   *  spaces are added on the right of the string.
   **/
  private String padRightWithSpaces(final String string, final int width) 
  {
    final int len = string.length();
    if(len == width)
      return string;

    final StringBuffer buffer = new StringBuffer(string);
    for(int i = 0 ; i < width - len; ++i) 
      buffer.append(' ');

    return buffer.toString();
  }

  /**
   *  Return the given string padded with spaces to the given width.  The
   *  spaces are added on the left of the string.
   **/
  private String padLeftWithSpaces(final String string, final int width) 
  {
    final int len = string.length();
    if(len == width) 
      return string;

    final StringBuffer buffer = new StringBuffer();

    for(int i = 0; i < width - len; ++i) 
      buffer.append(' ');

    buffer.append(string);
    return buffer.toString();
  }

  /**
   *  Return the height each line of the display should be.  Each feature will
   *  be drawn into one line.
   **/
  protected int getLineHeight() 
  {
    return getFontAscent() + 2;
  }

}
