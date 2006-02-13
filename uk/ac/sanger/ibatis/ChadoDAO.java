/* 
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2006  Genome Research Limited
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

package uk.ac.sanger.ibatis;

import java.sql.*;
import java.io.*;
import java.util.List;

public interface ChadoDAO
{

  public Feature getSequence(final int feature_id,
                             final String schema)
                        throws SQLException;

  public String getFeatureName(final int feature_id,
                               final String schema)
                       throws SQLException;

  public List getGff(final int parentFeatureID,
                     final String schema)
                     throws SQLException;

  public List getResidueFeatures(List cvterm_ids,
                                 final String schema)
                     throws SQLException;

  public List getResidueType(final String schema)
                     throws SQLException;

  public List getSchema()
              throws SQLException;

  public List getCvterm()
              throws SQLException;
}
