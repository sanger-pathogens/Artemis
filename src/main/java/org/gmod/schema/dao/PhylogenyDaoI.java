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

import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.phylogeny.Phylonode;
import org.gmod.schema.phylogeny.Phylotree;

import java.util.List;



public interface PhylogenyDaoI extends BaseDaoI {

    public Phylotree getPhyloTreeByName(String name);
    
    public List<Phylonode> getPhyloNodesByCvTermInTree(CvTerm type, Phylotree tree);
    
    public List<Phylonode> getAllPhylonodes();
    
    public List<Phylonode> getPhylonodeByDepthAndParent(double depth,Phylonode parent);
    
    public List<Phylonode> getPhylonodeByName(String name);
    
    public List<Phylonode> getPhylonodesByParent(Phylonode parent);
    
    //public List<Phylonode> getPhylonodeByDepth(double depth);
}
