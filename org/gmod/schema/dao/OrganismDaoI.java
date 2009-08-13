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

import org.gmod.schema.organism.Organism;

import java.util.List;

public interface OrganismDaoI extends BaseDaoI {

    /**
     * Get the organism corresponding to this id
     * 
     * @param id the organism id (primary key) to lookup by
     * @return the corresponding organism, or null
     */
    public abstract Organism getOrganismById(int id);

    /**
     * Get the organism corresponding to this common name 
     * 
     * @param commonName the short name to look up
     * @return the corresponding organism, or null
     */
    public abstract Organism getOrganismByCommonName(String commonName);

    /**
     * Get a list of the common name of all the organisms.
     * 
     * @return a (possibly empty) List<String> of all the organisms' common names
     */
    public abstract List<String> findAllOrganismCommonNames();
    
    /**
     * Get a list of all the  organisms 
     * 
     * @return a (possibly empty) List<Organism> of all the organisms'
     */
    public abstract List<Organism> getOrganisms();

}
