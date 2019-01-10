/* GFF3Encoder.java
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

import org.apache.commons.lang.StringUtils;

public class GFF3Encoder {
  public static final String mapChars[] = { "%", "&", ",", ";", "=", "\t", "\n", "\r" };
  public static final String mappedEscapes[] = { "%25", "%26", "%2C", "%3B", "%3D", "%09", "%0A", "%0D" };

  public static String decode(String s) {
    return StringUtils.replaceEach(s, mappedEscapes, mapChars);
  }

  public static String encode(String s) {
    return StringUtils.replaceEach(s, mapChars, mappedEscapes);
  }
}
