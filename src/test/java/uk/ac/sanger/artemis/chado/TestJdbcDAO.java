/* TestJdbcDAO.java
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

import static org.junit.Assert.*;

import java.sql.SQLException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.JPasswordField;

import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.pub.PubDbXRef;
import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;
import org.gmod.schema.sequence.FeatureCvTermPub;
import org.gmod.schema.sequence.Synonym;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Test suite for testing actual (non-mocked) JdbDAO database connectivity in situ.
 * It should be run on its own with suitable database connection parameters
 * specified, rather than as part of a build where any given database may not
 * necessarily be accessible.
 * NOTE 1: Testing values are tuned to a specific Sanger test database.
 * NOTE 2: I get the feeling JdbDAO methods are deprecated as ibatis is set by default
 * 		   in scripts etc, and there seem to be some bugs in it.
 * 
 * @author kp11
 *
 */
@SuppressWarnings({"unchecked"})
public class TestJdbcDAO implements ChadoTestConstants
{
	/** DAO to use */
	private static JdbcDAO dao;
	
	@BeforeClass
	public static void init() throws Exception
	{
		System.setProperty("chado", CHADO_CONNECTION_URL);
	
		final JPasswordField passwordField = new JPasswordField(CHADO_CONNECTION_PASSWORD);
		
		
		// Connection is created in the constructor bythe looks of it
		// and not closed until an application exit.
		dao = new JdbcDAO(CHADO_CONNECTION_URL, passwordField); 
		
	}
	
	@Test
	public void testIsFeatureCvTermRank()
	{
		// Not defined for JdbcDAO
	}
	
	@Test
	public void testGetLazyFeatureNoResiduesById()
	{
		// Not implemented in JdbcDAO
	}
	
	@Test
	public void testGetResidueFeatures()
	{
		// Not implemented in JdbcDAO
	}
	
	@Test
	public void testGetFeatureDbXRefsByFeatureId()
	{
		// Not implemented in JdbcDAO
	}
	
	@Test
	public void testGetFeatureById()
	{
		Feature feature = dao.getFeatureById(testFeatureIdForOrg);
		assertNotNull(feature);
		assertEquals(testFeatureIdForOrg.intValue(), feature.getFeatureId());
		assertEquals(testFeatureName, feature.getUniqueName());
		
		// This values does not appear to get set for JdbcDAO unlike IBatisDao !!
		// assertEquals(testFeatureType, feature.getCvTerm().getName());
	}
	
	@Test
	public void testGetFeatureByUniqueName()
	{
		Feature feature = dao.getFeatureByUniqueName(testFeatureName, testFeatureType);
		assertNotNull(feature);
		assertEquals(testFeatureIdForOrg.intValue(), feature.getFeatureId());
		assertEquals(testFeatureName, feature.getUniqueName());
		
		// This values does not appear to get set for JdbcDAO unlike IBatisDao !!
		// assertEquals(testFeatureType, feature.getCvTerm().getName());
		//assertEquals(testFeatureType, feature.getCvTerm().getName());
	}
	
	@Test
	public void testGetFeaturesByAnyCurrentName()
	{
		List<Feature> features = dao.getFeaturesByAnyCurrentName(testFeatureName) ;
		assertNotNull(features);
		assertEquals(1, features.size());
		assertEquals(testFeatureIdForOrg.intValue(), features.get(0).getFeatureId());
		assertEquals(testFeatureName, features.get(0).getUniqueName());
		
		// This values does not appear to get set for JdbcDAO unlike IBatisDao !!
		// assertEquals(testFeatureType, features.get(0).getCvTerm().getName());
	}
	
	@Test
	public void testGetFeatureCvTermsByFeature()
	{
		Feature feature = dao.getFeatureById(testFeatureIdForOrg);
		List<CvTerm> cvterms = dao.getFeatureCvTermsByFeature(feature);
		assertNotNull(cvterms);
		assertEquals(numCvTerms.intValue(), cvterms.size());
	}
	
	@Test
	public void testGetFeatureDbXRefsByFeatureUniquename()
	{
		List<DbXRef> refs = dao.getFeatureDbXRefsByFeatureUniquename(testFeatureName);
		assertNotNull(refs);
		assertEquals(numDbXRefs.intValue(), refs.size());
	}
	
	@Test
	public void testGetFeatureSynonymsByFeatureUniquename()
	{
		List<Synonym> synonyms = dao.getFeatureSynonymsByFeatureUniquename(testFeatureName);
		assertNotNull(synonyms);
		assertEquals(numSynonyms.intValue(), synonyms.size());
	}
	
	@Test
	public void testGetAllFeatureSynonymsAsFeature()
	{
		System.out.println("getAllFeatureSynonymsAsFeature fails for Ibatis but the method is not currently used.");
		
		// This method is not used anywhere so skip this for the moment.
		
		//List <Synonym> synonyms = dao. getAllFeatureSynonymsAsFeature();
		//assertNotNull(synonyms);
	}
	
	@Test
	public void testGetFeatureCvTermDbXRefByFeature()
	{
		Feature feature = dao.getFeatureById(testFeatureIdForOrg);
		List<FeatureCvTermDbXRef> results = dao.getFeatureCvTermDbXRefByFeature(feature);
		assertNotNull(results);
		assertEquals(numFeatureCvTermDbXRef.intValue(), results.size());
		
		Set<Integer> ids = new HashSet<Integer>();
		for (FeatureCvTermDbXRef obj : results)
		{
			ids.add(obj.getFeatureCvTerm().getFeatureCvTermId());
		}
		
		// WARNING: Database specific IDs!
		assertTrue(ids.contains(63075007));
		assertTrue(ids.contains(63075008));
		assertTrue(ids.contains(1162004));
		assertTrue(ids.contains(1162017));
		assertTrue(ids.contains(1162018));
	}
	
	@Test
	public void testGetFeatureCvTermPubByFeature()
	{
		Feature feature = dao.getFeatureById(testFeatureIdForOrg);
		List<FeatureCvTermPub> results = dao.getFeatureCvTermPubByFeature(feature);
		
		assertNotNull(results);
		assertEquals(numCvTermPub.intValue(), results.size());
	}
	
	@Test
	public void testGetSchema() throws SQLException
	{
		List<String> schemas = dao.getSchema();
		assertNotNull(schemas);
		assertEquals(numSchemas.intValue(), schemas.size());
	}
	
	@Test
	public void testGetCvTerms()
	{
		List<CvTerm> cvTerms = dao.getCvTerms();
		assertNotNull(cvTerms);
		assertEquals(numAllCvTerms.intValue(), cvTerms.size());
	}
	
	@Test
	public void testGetPubDbXRef()
	{
		List<PubDbXRef> refs =  dao.getPubDbXRef();
		assertNotNull(refs);
		assertEquals(numAllPubDbXRefs.intValue(), refs.size());
	}

}
