/* ChadoTestConstants.java
 *
 * created: 2018
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2018  Genome Research Limited
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
package uk.ac.sanger.artemis.chado;

/**
 * Purely a holder for shared unit test constants.
 * @author kp11
 *
 */
public interface ChadoTestConstants
{
	// Connection parameters for tests...
	final String CHADO_CONNECTION_USERNAME = "<Enter username here>";
	final String CHADO_CONNECTION_URL = "vm-pg-path-dev-02.internal.sanger.ac.uk:5432/pathogens_pg10?" + CHADO_CONNECTION_USERNAME;
	final String CHADO_CONNECTION_PASSWORD = "<Enter password here>";
	
	// Test parameters that may need to be tuned per database instance...
	final Integer testOrganismId = 27; // Pfalciparum
	final Integer testFeatureIdForOrg = 20760246;
	final String  testFeatureName = "PF3D7_0528500.1:pep";
	final String  testFeatureType = "polypeptide";
	final Integer numCvTerms = 6;
	final Integer numDbXRefs = 4;
	final Integer numSynonyms = 0;
	final Integer numFeatureCvTermDbXRef = 5;
	final Integer numCvTermPub = 0;
	final Integer numSchemas = 0;
	final Integer numAllCvTerms = 142926;
	final Integer numAllPubDbXRefs = 9122;
}
