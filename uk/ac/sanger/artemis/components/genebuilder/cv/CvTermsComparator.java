/* CvTermsComparator.java
 *
 * created: 2007
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2007  Genome Research Limited
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
 **/

package uk.ac.sanger.artemis.components.genebuilder.cv;

import java.util.Comparator;

import org.gmod.schema.cv.CvTerm;

public class CvTermsComparator implements Comparator
{
  public int compare(Object o1, Object o2)
  {
    CvTerm cvTerm1 = (CvTerm)o1;
    CvTerm cvTerm2 = (CvTerm)o2;
    
    return cvTerm1.getName().compareTo( cvTerm2.getName() );
  }   
}