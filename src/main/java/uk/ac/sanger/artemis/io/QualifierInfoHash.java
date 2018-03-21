/* QualifierInfoHash.java
 *
 * created: Fri Nov 21 2003
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2003  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/QualifierInfoHash.java,v 1.1 2004-06-09 09:50:10 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.StringVector;

import java.util.Hashtable;

/**
 *  A Hashtable of QualifierInfo objects.  The qualifier names are the keys.
 *  The QualifierInfo objects are the values.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: QualifierInfoHash.java,v 1.1 2004-06-09 09:50:10 tjc Exp $
 **/

public class QualifierInfoHash {
  /**
   *  Create a new (empty) QualifierInfoHash object.
   **/
  public QualifierInfoHash () {

  }

  /**
   *  Put the QualifierInfo object to the Hashtable with the qualifier name
   *  as the key and the QualifierInfo as the value.
   */
  public void put (final QualifierInfo info) {
    hash.put (info.getName (), info);

    name_cache = null;
  }

  /**
   *  Return the QualifierInfo that describes the Qualifier with the given
   *  name.
   **/
  public QualifierInfo get (final String name) {
    return (QualifierInfo) hash.get (name);
  }

  /**
   *  Return a StringVector containing the names of all the QualifierInfo
   *  objects.
   **/
  public StringVector names () {
    if (name_cache == null) {
      name_cache = new StringVector ();

      for (java.util.Enumeration e = hash.keys () ; e.hasMoreElements() ;) {
        name_cache.add ((String) e.nextElement());
      }
    }

    return name_cache;
  }

  /**
   *  Return a new copy of this object.
   **/
  public QualifierInfoHash copy () {
    final QualifierInfoHash new_hash = new QualifierInfoHash ();

    new_hash.hash = (Hashtable) hash.clone ();

    return new_hash;
  }

  /**
   *  Performs the same function as Vector.size ()
   */
  public int size () {
    return hash.size ();
  }

  /**
   *  Storage for QualifierInfo objects.
   */
  private Hashtable hash = new Hashtable ();

  /**
   *  A cache used by names() to
   **/
  private StringVector name_cache = null;
}
