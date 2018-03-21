/*
 * Copyright (c) 2006 Genome Research Limited.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published
 * by  the Free Software Foundation; either version 2 of the License or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with this program; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation Inc., 59 Temple Place - Suite 330,
 * Boston, MA  02111-1307 USA
 */

package org.gmod.schema.dao;

import java.util.List;

import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.pub.Pub;
import org.gmod.schema.pub.PubProp;
import org.gmod.schema.pub.PubDbXRef;

public interface PubDaoI extends BaseDaoI {

    /**
     * Retrieve the publication with this primary key
     * 
     * @param id the publication id
     * @return the corresponding publication, or null
     */
    public abstract Pub getPubById(int id);

    /**
     * Retrieve the publication with this primary key
     * 
     * @param uniqueName
     * @return the publication with this unique name, or null
     */
    public abstract Pub getPubByUniqueName(String uniqueName);

    
    public Pub getPubByDbXRef(DbXRef dbXRef);
    
    /**
     * Get a list of all PubDbXRef's
     * @return list of PubDbXRef's
     */
    public List<PubDbXRef> getPubDbXRef();
    
        /**
     * Retrieve the publication property with Pub and Cvterm
     * 
     * @param pub the Publication
     * @param cvTerm the cvTerm
     * @return the publication or null
     */
    public abstract List<PubProp> getPubPropByPubAndCvTerm(Pub pub,CvTerm cvTerm);
}
