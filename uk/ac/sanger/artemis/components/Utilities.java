/* Utilities.java
 *
 * created: Tue Sep 18 2001
 *
 * This file is part of Artemis
 * 
 * Copyright(C) 2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/Utilities.java,v 1.4 2007-02-22 19:45:18 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.EntrySourceVector;
import uk.ac.sanger.artemis.Options;

import java.awt.*;
import java.net.MalformedURLException;
import javax.swing.*;

/**
 *  Utilities methods used by the uk.ac.sanger.artemis.components package.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: Utilities.java,v 1.4 2007-02-22 19:45:18 tjc Exp $
 **/

public class Utilities 
{
  /**
   *  Return the JFrame that contains the given component.
   **/
  public static JFrame getComponentFrame(final JComponent component) 
  {
    if(component.getTopLevelAncestor() instanceof JFrame) 
      return(JFrame) component.getTopLevelAncestor();
    else 
      return null;
  }

  /**
   *  Move the given JFrame to the centre of the screen.
   **/
  public static void centreFrame(final JFrame frame) 
  {
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    final int x_position =(screen.width - frame.getSize().width) / 2;
    int y_position =(screen.height - frame.getSize().height) / 2;

    if(y_position < 10) 
      y_position = 10;

    frame.setLocation(new Point(x_position, y_position));
  }

   /**
   *  Move the given JFrame to the right of the screen.
   **/
  public static void rightJustifyFrame(final JFrame frame)
  {
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    int y_position =(screen.height - frame.getSize().height) / 2;

    if(y_position < 10)
      y_position = 10;

    frame.setLocation(new Point(0, y_position));
  }

  /**
   *  Move the given JFrame to the centre of the screen.
   **/
  public static void centreJustifyFrame(final JFrame frame, final int y_position)
  {
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    final int x_position =(screen.width - frame.getSize().width) / 2;

    frame.setLocation(new Point(x_position, y_position));
  }



  /**
   *  Find the parent Frame of the given component and re-centre it on the
   *  screen.
   **/
  public static void centreOwningFrame(final JComponent component) 
  {
    final JFrame frame = getComponentFrame(component);
    centreFrame(frame);
  }
  
  /**
   * Test to find if this is a multi-display/device environment
   * @return true if multiple displays found
   */
  public static boolean isMultiDisplay()
  {
    GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
    
    if(ge.getScreenDevices().length>1)
      return true;
    
    return false;
  }

  /**
   *  Returns a Vector containing the possible Entry sources for this
   *  application.
   *  @param frame The component that is creating the EntrySource objects.
   *   (Used for requesters.)
   *  @param listener InputStreamProgressEvent objects will be sent to this
   *    listener as progress on reading is made.
   **/
  public static EntrySourceVector getEntrySources(final JFrame frame,
                     final InputStreamProgressListener listener) 
  {
    final EntrySourceVector return_vector = new EntrySourceVector();
    
    return_vector.add(new FileDialogEntrySource(frame, listener));
    return_vector.add(new DbfetchEntrySource(frame));

    // return_vector.add(new BioJavaEntrySource());

    return return_vector;
  }
}
