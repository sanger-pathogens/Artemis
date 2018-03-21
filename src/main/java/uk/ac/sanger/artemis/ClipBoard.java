/* ClipBoard.java
 *
 * created: Mon Oct 26 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/ClipBoard.java,v 1.1 2004-06-09 09:44:14 tjc Exp $
 **/

package uk.ac.sanger.artemis;

/**
 *  A simple clipboard for Diana which can interact with the system clipboard.
 *  Any Object can be stored in this clipboard but only some will be useful to
 *  other programs.
 *
 *  @author Kim Rutherford
 *  @version $Id: ClipBoard.java,v 1.1 2004-06-09 09:44:14 tjc Exp $
 **/

public class ClipBoard {
  /**
   *  Create a new ClipBoard object.
   **/
  public ClipBoard () {

  }

  /**
   *  Set the contents of the clipboard.
   *  @param clip The new clipboard contents.
   **/
  public void setClip (Selection clip) {
    this.clip = clip;
  }

  
  /**
   *  Return the current contents of the clipboard.
   **/
  public Selection getClip () {
    return clip;
  }


  /**
   *  Holds whatever was last clipped.
   **/
  Selection clip;
}


