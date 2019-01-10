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

public interface BaseDaoI {

    
    /**
     * Save the object to the database (at the end of the current transaction, 
     * or depending upon flush mode). This method is available in all the DAOs. 
     * It's recommended to call it through an appropriate one eg SequenceDaoI 
     * for FeatureI 
     * @param o the new instance to persist
     */
    public abstract void persist(Object o);
    
    /**
     * Delete an object from the database
     * @param o the object to remove
     * 
     */
    public abstract void delete(Object o);
    
    /**
     * Merge (update) an already persistent object back to the database (at the end of 
     * the current transaction, or depending upon flush mode). This method is defined in 
     * all the DAOs. It's recommended to call it through an appropriate one eg SequenceDaoI
     * for FeatureI 
     * 
     * @param o the existing feature to merge
     */
    public abstract void merge(Object o);
    
}
