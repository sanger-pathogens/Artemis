/* SelectionMenu.java
 *
 * created: Sun Jan 10 1999
 *
 * This file is part of Artemis
 * 
 * Copyright(C) 1998,1999,2000,2001,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/SelectionMenu.java,v 1.1 2004-06-09 09:47:44 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.*;

import java.awt.event.*;
import javax.swing.*;

/**
 *  This a super class for EditMenu, ViewMenu and GotoMenu.  It is a JMenu that
 *  knows how to get hold of the Selection.  It also has the method
 *  getParentFrame() to find the owning JFrame of the menu.
 *
 *  @author Kim Rutherford
 *  @version $Id: SelectionMenu.java,v 1.1 2004-06-09 09:47:44 tjc Exp $
 **/

public class SelectionMenu extends JMenu 
{

  /** The Selection that was passed to the constructor. */
  /* final */ private Selection selection;

  /** The JFrame reference that was passed to the constructor. */
  private JFrame frame = null;

  /**
   *  Create a new SelectionMenu object.
   *  @param frame The JFrame that owns this JMenu.
   *  @param name The string to use as the menu title.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   **/
  public SelectionMenu(final JFrame frame,
                        final String menu_name,
                        final Selection selection) 
  {
    super(menu_name);
    this.frame = frame;
    this.selection = selection;
  }

  /**
   *  Return the JFrame that was passed to the constructor.
   **/
  public JFrame getParentFrame() 
  {
    return frame;
  }

  /**
   *  Check that the are some Features in the current selection.  If there are
   *  some features then return true.  If there are no features selected then
   *  popup a message saying so and return false.
   **/
  protected boolean checkForSelectionFeatures() 
  {
    return checkForSelectionFeatures(getParentFrame(), getSelection());
  }

  /**
   *  Check that the are some Features in the given selection.  If there are
   *  some features then return true.  If there are no features selected then
   *  popup a message saying so and return false.
   *  @param frame The JFrame to use for MessageDialog components.
   *  @param selection The selected features to check.
   **/
  static boolean checkForSelectionFeatures(final JFrame frame,
                                           final Selection selection) 
  {
    final FeatureVector features_to_check = selection.getAllFeatures();

    if(features_to_check.size() == 0) 
    {
      new MessageDialog(frame, "No features selected");
      return false;
    } 
    else 
      return true;
  }

  /**
   *  Check that the are some Features in the current selection and that the
   *  selection isn't too big.  If there is more than one feature in the
   *  selection and less than the given maximum then return true.  If there
   *  are no features selected then popup a message saying so and return
   *  false.  If there are more selected features than the given maximum
   *  than display the message in a YesNoDialog component and return the
   *  result of the dialog.
   **/
  protected boolean checkForSelectionFeatures(final int maximum_size,
                                              final String message) 
  {
    return checkForSelectionFeatures(getParentFrame(), getSelection(),
                                      maximum_size, message);
  }

  /**
   *  Check that the are some Features in the given selection and that the
   *  selection isn't too big.  If there is more than one feature in the
   *  selection and less than the given maximum then return true.  If there
   *  are no features selected then popup a message saying so and return
   *  false.  If there are more selected features than the given maximum
   *  than display the message in a YesNoDialog component and return the
   *  result of the dialog.
   *  @param frame The JFrame to use for MessageDialog components.
   *  @param selection The selected features to check.
   **/
  static boolean checkForSelectionFeatures(final JFrame frame,
                                           final Selection selection,
                                           final int maximum_size,
                                           final String message) 
  {
    final FeatureVector features_to_check = selection.getAllFeatures();

    if(features_to_check.size() == 0) 
    {
      new MessageDialog(frame, "No features selected");
      return false;
    }
    else 
    {
      if(features_to_check.size() > maximum_size) 
      {
        final YesNoDialog dialog =
          new YesNoDialog(frame, message);

        return dialog.getResult();
      }
      else 
        return true;
    }
  }

  /**

   *  Check that there are only CDS features in the given selection, that there
   *  are some features selected and that the selection isn't too big.  If
   *  there is more than one feature in the selection and less than the given
   *  maximum then return true.  If there are no features selected then popup
   *  a message saying so and return false.  If there are more selected
   *  features than the given maximum than display the message in a
   *  YesNoDialog component and return the result of the dialog.
   *  @param frame The Frame to use for MessageDialog components.
   *  @param selection The selected features to check.
   **/
  static boolean checkForSelectionCDSFeatures(final JFrame frame,
                                               final Selection selection,
                                               final int maximum_size,
                                               final String max_message) 
  {
    final FeatureVector features_to_check = selection.getAllFeatures();

    if(features_to_check.size() == 0) 
    {
      new MessageDialog(frame, "No CDS features selected");
      return false;
    }
    else
    {
      for(int i = 0 ; i < features_to_check.size() ; ++i) 
      {
        final Feature this_feature = features_to_check.elementAt(i);
        if(!this_feature.isCDS()) 
        {
          final String message =
            "a non-CDS feature(" + this_feature.getIDString() +
            ") is selected - can't continue";
          final MessageDialog dialog =
            new MessageDialog(frame, message);
          return false;
        }
      }

      if(features_to_check.size() > maximum_size) 
      {
        final YesNoDialog dialog =
          new YesNoDialog(frame, max_message);

        return dialog.getResult();
      }
      else 
        return true;
    }
  }

  /**
   *  Check that the are some FeatureSegments in the current selection and
   *  that the selection isn't too big.  If there is more than one
   *  FeatureSegments in the selection and less than the given maximum then
   *  return true.  If there are no FeatureSegments selected then popup a
   *  message saying so and return false.  If there are more selected
   *  FeatureSegments than the given maximum than display the message in a
   *  YesNoDialog component and return the result of the dialog.
   **/
  protected boolean checkForSelectionFeatureSegments(final int maximum_size,
                                                     final String message) 
  {
    final FeatureSegmentVector segments_to_check =
      getSelection().getSelectedSegments();

    if(segments_to_check.size() == 0) 
    {
      new MessageDialog(getParentFrame(), "No exons selected");
      return false;
    } 
    else 
    {
      if(segments_to_check.size() > maximum_size) 
      {
        final YesNoDialog dialog =
          new YesNoDialog(getParentFrame(), message);

        return dialog.getResult();
      } 
      else 
        return true;
    }
  }

  /**
   *  Check that the current selection contains a MarkerRange.  If it does
   *  then return true.  If not then popup a message saying so and return
   *  false.
   **/
  protected boolean checkForSelectionRange() 
  {
    return checkForSelectionRange(getParentFrame(), getSelection());
  }

  /**
   *  Check that the current selection contains a MarkerRange.  If it does
   *  then return true.  If not then popup a message saying so and return
   *  false.
   *  @param frame The JFrame to use for MessageDialog components.
   *  @param selection The selected features to check.
   **/
  static boolean checkForSelectionRange(final JFrame frame,
                                        final Selection selection) 
  {
    final MarkerRange marker_range = selection.getMarkerRange();

    if(marker_range == null) 
    {
      new MessageDialog(frame, "No bases selected");
      return false;
    } 
    else 
      return true;
  }

  /**
   *  Return the Selection object that was passed to the constructor.
   **/
  protected Selection getSelection() 
  {
    return selection;
  }
  
  /**
   *
   **/
  protected static KeyStroke makeMenuKeyStroke(final int key_code)
  {
    return KeyStroke.getKeyStroke(key_code, InputEvent.CTRL_MASK);
  }

}

