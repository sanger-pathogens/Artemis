/* GFF3AttributeAggregator.java
 * *
 * This file is part of Artemis
 *
 * Copyright (C) 2014 Genome Research Limited
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
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.StringVector;

/**
 * Abstraction for processing various Artemis qualifiers into usable GFF3 strings.
 */
public interface GFF3AttributeAggregator {
  /**
   * Prepare a set of string values as a GFF attribute value.
   * @param values the value set to convert to a <code>String</code>
   * @return  the <code>String</code> representation
   */
  public abstract String process(StringVector values);
  
}
