/* BasePlotGroup.java
 *
 * created: Tue Dec 15 1998
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998-2008  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/BasePlotGroup.java,v 1.9 2008-06-16 12:11:01 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.plot.*;

import java.awt.*;
import java.util.Vector;

import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSplitPane;

/**
 *  This is a super-component containing several BasePlot components, each of
 *  which can toggled off and on.
 *
 *  @author Kim Rutherford
 *  @version $Id: BasePlotGroup.java,v 1.9 2008-06-16 12:11:01 tjc Exp $
 **/

public class BasePlotGroup extends JPanel
                           implements DisplayAdjustmentListener 
{
  private static final long serialVersionUID = 1L;

  /**
   *  The EntryGroup that contains the sequence that this JComponent is
   *  displaying.
   **/
  private final EntryGroup entry_group;

  /**
   *  This array contains the Algorithm objects of the BasePlot components in
   *  this BasePlotGroup, as set by the constructor.
   **/
  private final Vector<BaseAlgorithm> plot_value_producers = new Vector<BaseAlgorithm> ();

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
    this.entry_group = entry_group;
    this.selection = selection;
    this.goto_event_source = goto_event_source;

    final Strand forward_strand =
      entry_group.getBases().getForwardStrand();

    // obtain the reverse complement
    final Strand reverse_strand =
      entry_group.getBases().getReverseStrand();

    setLayout(gridbag);
    setBackground(Color.WHITE);

    c.fill = GridBagConstraints.BOTH;
    c.anchor = GridBagConstraints.NORTH;
    c.gridwidth = GridBagConstraints.REMAINDER;
    c.weightx = 1;
    c.weighty = 1;
    c.insets = new Insets(0,0,2,0);

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
    

    //CumulativeATSkewAlgorithm
    addAlgorithm(new CumulativeATSkewAlgorithm(forward_strand));
    addAlgorithm(new CumulativeGCSkewAlgorithm(forward_strand));

    //Positional Asymmetry
    addAlgorithm(new PositionalAsymmetryAlgorithm(forward_strand));

    //Informational Entropy
    addAlgorithm(new EntropyAlgorithm(forward_strand));

    //Scaled Chi
    addAlgorithm(new ScaledChiAlgorithm(forward_strand));
    addAlgorithm(new ScaledChiAlgorithm(reverse_strand));

    //Corrected Scaled Chi Square
    addAlgorithm(new CSCSAlgorithm(forward_strand));
    addAlgorithm(new CSCSAlgorithm(reverse_strand));

    //Mutational Response Index
    addAlgorithm(new MRIAlgorithm(forward_strand));
    addAlgorithm(new MRIAlgorithm(reverse_strand));

    //Effective Codon Number
    addAlgorithm(new NcAlgorithm(forward_strand));
    addAlgorithm(new NcAlgorithm(reverse_strand));

    //Intrinsic Codon Deviation Index
    addAlgorithm(new ICDIAlgorithm(forward_strand));
    addAlgorithm(new ICDIAlgorithm(reverse_strand));
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
        bp.paintComponent(g);
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
    final StringBuffer closingPlots = new StringBuffer();
    final Component[] children = getComponents();
    int nvis = 0;
    for(int i = 0 ; i<children.length ; ++i)
    {
      if(children[i] instanceof BasePlot) 
      {
        // if this is an indexed sequence change hide any
        // userplots that are not indexed
        if(event.getType() == DisplayAdjustmentEvent.IDX_SEQUENCE_CHANGE)
        {
          Algorithm alg = ((BasePlot)children[i]).getAlgorithm();
          if(alg instanceof UserDataAlgorithm &&
            ((UserDataAlgorithm)alg).FORMAT !=  UserDataAlgorithm.TABIX_INDEXED_FORMAT &&
            findPlotByAlgorithm(alg).isVisible())
          {
            closingPlots.append(alg.getAlgorithmName()+"\n");
            children[i].setVisible(false);
            continue;
          }
        }

        if(children[i].isVisible())
          nvis++;
        ((BasePlot)children[i]).displayAdjustmentValueChanged(event);
      }
    }
    
    if(event.getType() == DisplayAdjustmentEvent.IDX_SEQUENCE_CHANGE &&
       closingPlots.length() > 0)
    {
      if(nvis == 0 && getParent() instanceof JSplitPane)
      {
        JSplitPane splitPane = (JSplitPane) getParent();
        splitPane.setDividerSize(0);
        splitPane.setDividerLocation(0);
      }

      JOptionPane.showMessageDialog(this, 
          closingPlots.toString()+
          "\nAs the sequence is changing the above user plot(s) are closing as they are\n"+
          "not indexed with multiple sequences. You can load in the corresponding plot\n"+
          "for the new sequence.", 
          "Closing Userplot", JOptionPane.INFORMATION_MESSAGE);
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
      final BaseAlgorithm this_algorithm = plot_value_producers.elementAt(i);
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
      final BaseAlgorithm this_algorithm = plot_value_producers.elementAt (i);
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
   * Return the number of visible plots
   * @return
   */
  public int getVisibleCount()
  {
    int cnt = 0;
    Component comp[] = getComponents();
    for(int i = 0 ; i<comp.length ; ++i)
      if(comp[i] instanceof BasePlot)
      {
        if(comp[i].isVisible())
          cnt++;
      }
    return cnt;
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
      new BasePlot(algorithm, getSelection(), getGotoEventSource(), entry_group);

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
   *  Return the GotoEventSource object that was passed to the constructor.
   **/
  private GotoEventSource getGotoEventSource()
  {
    return goto_event_source;
  }
}
