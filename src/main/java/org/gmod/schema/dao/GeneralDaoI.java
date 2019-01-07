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

import org.gmod.schema.analysis.AnalysisFeature;
import org.gmod.schema.general.Db;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.sequence.Feature;



public interface GeneralDaoI extends BaseDaoI{

    /**
     * Retrieve a database by name
     * 
     * @param name the name to lookup
     * @return the corresponding db, or null
     */
    public abstract Db getDbByName(String name);

    /**
     * Retrieve the db xref corresponding to a given DB and accession number
     * 
     * @param db the db the dbxref refers to
     * @param accession the accession "number" the dbxref refers to
     * @return the dbxref, or null
     */
    public abstract DbXRef getDbXRefByDbAndAcc(Db db, String accession);

    /**
     * Retrieve the analysisfeature corresponding to the given feature
     * 
     * @param feature the feature whose analysisfeature has to be found
     * @return the analysisfeature
     */
    public abstract AnalysisFeature getAnalysisFeatureFromFeature(Feature feature);
}
