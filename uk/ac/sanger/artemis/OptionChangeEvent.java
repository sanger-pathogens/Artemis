/* OptionChangeEvent.java
 *
 * created: Wed Jun 13 2001
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/OptionChangeEvent.java,v 1.1 2004-06-09 09:44:58 tjc Exp $
 */

package uk.ac.sanger.artemis;

/**
 *  This event is sent when an option changes.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: OptionChangeEvent.java,v 1.1 2004-06-09 09:44:58 tjc Exp $
 **/

public class OptionChangeEvent extends ChangeEvent {
  /**
   *  Create a new OptionChangeEvent object.
   *  @param source The object that generated the event.
   *  @param option_name The name of option that changed
   **/
  public OptionChangeEvent (final Object source,
                            final String option_name)
  {
    super (source);
    this.option_name = option_name;
  }

  /**
   *  Return the name of the option that changed.
   **/
  public String getOptionName () {
    return option_name;
  }

  /**
   *  The option name that was passed to the constructor.
   **/
  private String option_name;
}
