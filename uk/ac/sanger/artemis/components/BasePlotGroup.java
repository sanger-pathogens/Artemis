/* BasePlotGroup.java
 *
 * created: Tue Dec 15 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/BasePlotGroup.java,v 1.3 2004-11-09 16:21:16 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.plot.*;

import java.awt.*;
import java.awt.event.*;
import java.io.IOException;
import java.io.File;
import java.util.Vector;

import javax.swing.*;

/**
 *  This is a super-component containing several BasePlot components, each of
 *  which can toggled off and on.
 *
 *  @author Kim Rutherford
 *  @version $Id: BasePlotGroup.java,v 1.3 2004-11-09 16:21:16 tjc Exp $
 **/

public class BasePlotGroup extends JPanel
                           implements DisplayAdjustmentListener 
{
  /**
   *  Create a new BasePlotGroup component for the given EntryGroup.
   *  @param selection Used to set and display the current selection in the
   *    BasePlot.
   *  @param goto_event_source The object the we will call gotoBase () on.
   *    This allows the user to double click on a base in a BasePlot and have
   *    the FeatureDisplay follow.
   **/
  public BasePlotGroup(final EntryGroup entry_group,
                       final Component owning_component,
                       final Selection selection,
                       final GotoEventSource goto_event_source) 
  {
    this.owning_component = owning_component;
    this.entry_group = entry_group;
    this.selection = selection;
    this.goto_event_source = goto_event_source;

    final Strand forward_strand =
      entry_group.getBases().getForwardStrand();

    final Strand reverse_strand =
      entry_group.getBases().getReverseStrand();

    // the following arrays are from failed tests
    final float [] [] test_weights = {
      {1, 0,39,99,11}, // a
      {76,8,15,1,45}, // c
      {2,0,42,0,6}, // g
      {21,91,4,0,38}  // t
    };

    final float [] [] test_weights2 = {
      {11,11,10, 8,11,10,11,11, 7, 8,25, 3,100,  0,27},
      {29,33,30,30,32,34,37,38,39,36,26,75,  0,  0,14},
      {14,12,10,10, 9,11,10, 9, 7, 6,26, 1,  0,100,49},
      {46,44,50,52,48,45,42,43,47,51,23,21,  0,  0,10}
    };

    final float [] [] test_weights3 = {
      {0,0,1,0,0,0},
      {0,0,0,1,0,0},
      {1,0,0,0,1,0},
      {0,1,0,0,0,1},
    };

    final float [] cai_test = {
      0.113F, 1.0F,   0.117F, 1.0F,
      1.0F,   0.693F, 0.036F, 0.005F,
      0.071F, 1.0F,   0.0F,   0.0F,
      1.0F,   0.077F, 0.0F,   1.0F,

      0.006F, 0.003F, 0.039F, 0.003F,
      0.047F, 0.009F, 1.0F,   0.002F,
      0.245F, 1.0F,   1.0F,   0.007F,
      0.137F, 0.002F, 0.002F, 0.002F,

      0.823F, 1.0F,   0.003F, 1.0F,
      0.921F, 1.0F,   0.012F, 0.006F,
      0.053F, 1.0F,   0.135F, 1.0F,
      0.021F, 0.031F, 1.0F,   0.003F,

      1.0F,   0.831F, 0.002F, 0.018F,
      1.0F,   0.316F, 0.015F, 0.001F,
      0.554F, 1.0F,   1.0F,   0.016F,
      1.0F,   0.02F,  0.002F, 0.004F
    };

//    plot_algorithms_temp.add(new CAIWindowAlgorithm(forward_strand, cai_test));
//    plot_algorithms_temp.add(new UserBaseAlgorithm(forward_strand, "test", test_weights));
//    plot_algorithms_temp.add(new UserBaseAlgorithm(forward_strand, "test2", test_weights2));
//    plot_algorithms_temp.add(new UserBaseAlgorithm(forward_strand, "test3", test_weights3));

    setLayout(gridbag);

    c.fill = GridBagConstraints.HORIZONTAL;
    c.anchor = GridBagConstraints.NORTH;
    c.gridwidth = GridBagConstraints.REMAINDER;
    c.gridheight = 1;
    c.weightx = 1;
    c.weighty = 0;

    c.insets = new Insets(0,0,5,0);

    addAlgorithm(new GCWindowAlgorithm(forward_strand));
    addAlgorithm(new GCSDWindowAlgorithm(forward_strand));
    addAlgorithm(new AGWindowAlgorithm(forward_strand));
    addAlgorithm(new GCFrameAlgorithm(forward_strand));
    addAlgorithm(new GCFrameAlgorithm(reverse_strand));
    addAlgorithm(new Codon12CorrelationAlgorithm(forward_strand));
    addAlgorithm(new Codon12CorrelationAlgorithm(reverse_strand));
    addAlgorithm(new GCDeviationAlgorithm(forward_strand));
    addAlgorithm(new ATDeviationAlgorithm(forward_strand));
    addAlgorithm(new KarlinSigAlgorithm(forward_strand));
  }

 
  /**
  * 
  * Routine for printing graphs in artemis
  *
  */
  protected void printComponent(Graphics g)
  {
    final Component[] children = getComponents();
    for(int i = 0 ; i<children.length ; ++i)
      if(children[i] instanceof BasePlot)
      {
        BasePlot bp = (BasePlot)children[i];
        if(!bp.isVisible())
          continue;
        bp.paintCanvas(g);
        g.translate(0,bp.getHeight());
      }
  }


  /**
  *
  * Return the number of visible BasePlot objects
  *
  */
  protected int getNumberBasePlots()
  {
    final Component[] children = getComponents();
    int num = 0;

    for(int i = 0 ; i<children.length ; ++i)
      if(children[i] instanceof BasePlot)
      {
        BasePlot bp = (BasePlot)children[i];
        if(bp.isVisible())
          num++;
      }
    return num;
  }

  /**
   *  Implementation of the DisplayAdjustmentListener interface.  Invoked when
   *  a component (FeatureDisplay) scrolls or changes the scale.  Sends the
   *  event to all the BasePlot components in this BasePlotGroup.
   **/
  public void displayAdjustmentValueChanged(DisplayAdjustmentEvent event) 
  {
    final Component[] children = getComponents();

    for(int i = 0 ; i<children.length ; ++i)
    {
      if(children[i] instanceof BasePlot) 
        ((BasePlot)children[i]).displayAdjustmentValueChanged(event);
    }
  }

  /**
   *  Add a new BasePlot component for the given algorithm.
   *  @return The new BasePlot
   **/
  public BasePlot addAlgorithm(final BaseAlgorithm algorithm) 
  {
    plot_value_producers.addElement(algorithm);

    return makePlot(algorithm, gridbag, c);
  }

  /**
   *  Return the first CodonUsageAlgorithm or null if there isn't one in the
   *  group.
   **/
  public CodonUsageAlgorithm getCodonUsageAlgorithm() 
  {
    for(int i = 0 ; i < plot_value_producers.size() ; ++i) 
    {
      final BaseAlgorithm this_algorithm =
        (BaseAlgorithm) plot_value_producers.elementAt(i);

      if(this_algorithm instanceof CodonUsageAlgorithm) 
        return (CodonUsageAlgorithm) this_algorithm;
    }

    return null;
  }

  /**
   *  Return an array containing the Algorithm objects of the BasePlot
   *  components in this BasePlotGroup.
   **/
  public BaseAlgorithm[] getPlotAlgorithms() 
  {
    final BaseAlgorithm[] return_array =
      new BaseAlgorithm[plot_value_producers.size ()];

    for(int i = 0 ; i < plot_value_producers.size () ; ++i) 
    {
      final BaseAlgorithm this_algorithm =
        (BaseAlgorithm) plot_value_producers.elementAt (i);
      return_array[i] = this_algorithm;
    }

    return return_array;
  }

  /**
   *  Return true if and only if the BasePlot for the given Algorithm is
   *  visible.
   **/
  public boolean basePlotIsVisible(final Algorithm algorithm) 
  {
    final Component base_plot = findPlotByAlgorithm(algorithm);
    return base_plot.isVisible ();
  }

  /**
   *  Given an Algorithm, find and set the visibility of the corresponding
   *  BasePlot component.
   **/
  public void setVisibleByAlgorithm(final Algorithm algorithm,
                                    final boolean visible) 
  {
    final JComponent base_plot = findPlotByAlgorithm(algorithm);

    base_plot.setVisible(visible);
    if (getParent () != null) {
      // XXX change to revalidate().
      getParent ().validate ();
    }
  }

  /**
   *  Find the BasePlot component in this BasePlotGroup object that is
   *  plotting the given Algorithm or null if no such BasePlot object exists.
   **/
  private JComponent findPlotByAlgorithm(final Algorithm algorithm) 
  {
    final Component[] children = getComponents();

    for(int i = 0 ; i < children.length ; ++i) 
    {
      if(children[i] instanceof BasePlot) 
      {
        final Algorithm component_algorithm =
                           ((BasePlot)children[i]).getAlgorithm ();
        if(component_algorithm == algorithm)
          return (JComponent)children[i];
      }
    }

    return null;
  }


  /**
   *  Create a Plot component for the given Algorithm and then make a button
   *  for it so that it can be shown and hidden.  Note that button making is
   *  currently disabled
   *  @param algorithm The Algorithm to create a Plot of.
   *  @param gridbag The GridBagLayout to use to lay out the Plot components.
   *  @param constraints The GridBagConstraints object to use to lay out the
   *    Plot components.
   *  @return The new BasePlot
   **/
  private BasePlot makePlot(BaseAlgorithm algorithm,
                            GridBagLayout gridbag,
                            GridBagConstraints constraints) 
  {
    final BasePlot new_base_plot =
      new BasePlot(algorithm, getSelection(), getGotoEventSource());

    gridbag.setConstraints(new_base_plot, constraints);
    add(new_base_plot);
    new_base_plot.setVisible(false);

    getSelection().addSelectionChangeListener(new_base_plot);

    if (getParent () != null) {
      // XXX change to revalidate().
      getParent ().validate ();
    }

    return new_base_plot;
  }

  /**
   *  Return the Selection object that was passed to the constructor.
   **/
  private Selection getSelection() 
  {
    return selection;
  }

  /**
   *  Return the EntryGroup object that was passed to the constructor.
   **/
  private EntryGroup getEntryGroup() 
  {
    return entry_group;
  }

  /**
   *  Return the GotoEventSource object that was passed to the constructor.
   **/
  private GotoEventSource getGotoEventSource()
  {
    return goto_event_source;
  }

  /**
   *  The EntryGroup that contains the sequence that this JComponent is
   *  displaying.
   **/
  private final EntryGroup entry_group;

  /**
   *  This array contains the Algorithm objects of the BasePlot components in
   *  this BasePlotGroup, as set by the constructor.
   **/
  private final Vector plot_value_producers = new Vector ();

  /**
   *  The layout object used by this component.
   **/
  private final GridBagLayout gridbag = new GridBagLayout();

  /**
   *  The constraints object used by this component.
   **/
  private final GridBagConstraints c = new GridBagConstraints();

  /**
   *  The Selection that was passed to the constructor.
   **/
  private Selection selection;

  /**
   *  The GotoEventSource that was passed to the constructor.
   **/
  private GotoEventSource goto_event_source;

  /**
   *  This is a reference to the parent Container of this BasePlot object.
   **/
  private Component owning_component;
}
