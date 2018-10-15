/* CodonWeight.java
 *
 * created: Tue Apr 13 1999
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/CodonWeight.java,v 1.1 2004-06-09 09:51:23 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

/**
 *  Objects of this class are the source of the codon usage data that is used
 *  by the methods in the CodonWindowAlgorithm class.
 *
 *  @author Kim Rutherford
 *  @version $Id: CodonWeight.java,v 1.1 2004-06-09 09:51:23 tjc Exp $
 **/

abstract public class CodonWeight {
  /**
   *  Return the codon weight for the given codon.
   *  @param codon XXX A lowercase string containing the bases of the codon to
   *    look up.
   **/
  abstract public float getCodonValue (final char base1, final char base2,
                                       final char base3);
  

  /**
   *  Returns the name of the file that the weights were read from.
   **/
  abstract public String getName ();
}
