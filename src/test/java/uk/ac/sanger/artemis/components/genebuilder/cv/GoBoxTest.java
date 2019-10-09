/* GoBoxTest.java
 *
 * created: 2019
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2019  Genome Research Limited
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

import static org.junit.Assert.*;

import org.gmod.schema.cv.Cv;
import org.gmod.schema.cv.CvTerm;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.mockito.Mockito;
import org.mockito.MockitoAnnotations;

import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.StringVector;

import static org.mockito.Mockito.*;

import java.awt.Dimension;


/**
 * Unit test for the GoBox class.
 * The tests are currently targeted specifically 
 * at testing RT ticket 621414.
 * 
 * @author kp11
 *
 */
public class GoBoxTest
{
	// Test Go term strings
	
	private final String GO_TERM_STRING_1 = 
			"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity";
	
	private final String GO_TERM_STRING_ALL_TERMS = 
			"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;with=with-value1;db_xref=CG42234;qualifier=qual-value1;evidence=Inferred from Sequence Model;" +
					AbstractCvBox.ASSIGNEDBY_QUALIFIER +
					"EuPathDB;date=20190220";
	
	private final String GO_TERM_STRING_WITH_DBXREF_FIELD = 
		"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;db_xref=CG42234;evidence=Inferred from Sequence Model;date=20190220";
	
	private final String GO_TERM_STRING_WITH_WITH_FIELD = 
		"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;with=with-value1;evidence=Inferred from Sequence Model;date=20190220";
	
	private final String GO_TERM_STRING_WITH_QUALIFIER_FIELD = 
		"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;qualifier=qual-value1;evidence=Inferred from Sequence Model;date=20190220";
	
	private final String GO_TERM_STRING_WITH_EVIDENCE_FIELD = 
			"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;evidence=Inferred from High Throughput Direct Assay;date=20190220";
	
	private final String GO_TERM_STRING_WITH_SOURCE_FIELD = 
			"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;"+AbstractCvBox.ASSIGNEDBY_QUALIFIER+"EuPathDB;date=20190220";
	
	private final static String MOLECULAR_FUNCTION_KEY = "molecular_function";
	
	private final static String GENEDB_SOURCE = "GeneDB";
	
	
	/**
	 * Force headless mode for these tests.
	 */
	@BeforeClass
	public static void setHeadlessMode()
	{
		System.setProperty("java.awt.headless", "true");
	}
	
	/**
	 * Helper method to create a GO Box spy for testing.
	 * 
	 * @param qualifierText String
	 * @param cvName String
	 * @return GoBox spy object
	 */
	private GoBox getGoBoxForTest(final String qualifierText, final String cvName)
	{
		final GoBox goBox = new GoBox(
				new Qualifier("GO", qualifierText),
				qualifierText,
				0,
				null,
				new Dimension(100,100)
				)
		{
			@Override
			public CvTerm getGOCvTermFromChado(String term)
			{
				CvTerm cvterm = new CvTerm();
				Cv cv = new Cv();
				cv.setName(cvName);
				cvterm.setCv(cv);
				return cvterm;
			}
		};
		
		GoBox spy = Mockito.spy(goBox);
		
		return spy;
		
	}
	
	/**
	 * Initialisation for each test.
	 * @throws Exception
	 */
	@Before
	public void setUp() throws Exception
	{
		MockitoAnnotations.initMocks(this);
	}
	
	/**
	 * Test load of database sources from options file.
	 */
	@Test
	public void testSourcesLoad()
	{
		assertNotNull(GoBox.sources);
		assertTrue(GoBox.sources.size() > 1);
		assertEquals(GoBox.SOURCE_DEFAULT_OPTION, GoBox.sources.elementAt(0));
		assertTrue(GoBox.sources.contains(GENEDB_SOURCE));
	}
	
	/**
	 * Test the getSourceIndex() method.
	 */
	@Test
	public void testGetSourceIndex() 
	{
		assertEquals(-1, GoBox.getSourceIndex(null));
		assertEquals(0, GoBox.getSourceIndex(GoBox.SOURCE_DEFAULT_OPTION));
		assertEquals(-1, GoBox.getSourceIndex("FOOBAR"));
		assertTrue(GoBox.getSourceIndex(GENEDB_SOURCE) > -1);
	}
	
	// ==============================================================
	// isQualifierChanged tests
	// ==============================================================
	
	/**
	 * Test case when no GO string change has occurred.
	 */
	@Test
	public void testIsQualifierChangedForNoChange()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_1, MOLECULAR_FUNCTION_KEY);
		
		// When
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertFalse(isChanged);
	}
	
	// With field
	
	/**
	 * Test case when we have no current "with" value but user has entered one
	 */
	@Test
	public void testIsQualifierChangedForWithField1()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_1, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getWithTextFieldValue()).thenReturn("with-value1");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	/**
	 * Test case when we have a current "with" value that has not been changed by user
	 */
	@Test
	public void testIsQualifierChangedForWithField2()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_WITH_FIELD, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getWithTextFieldValue()).thenReturn("with-value1");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertFalse(isChanged);
	}
	
	/**
	 * Test case when we have a current "with" value that has been changed by user
	 */
	@Test
	public void testIsQualifierChangedForWithField3()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_WITH_FIELD, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getWithTextFieldValue()).thenReturn("with-value2");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	// Dbxref field
	
	/**
	 * Test case when we have no current "Dbxref" value but user has entered one
	 */
	@Test
	public void testIsQualifierChangedForDbxrefField1()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_1, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getDbxrefTextFieldValue()).thenReturn("dbxref-value1");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	/**
	 * Test case when we have a current "Dbxref" value that has not been changed by user
	 */
	@Test
	public void testIsQualifierChangedForDbxrefField2()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_DBXREF_FIELD, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getDbxrefTextFieldValue()).thenReturn("CG42234");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertFalse(isChanged);
	}
	
	/**
	 * Test case when we have a current "Dbxref" value that has been changed by user
	 */
	@Test
	public void testIsQualifierChangedForDbxrefField3()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_DBXREF_FIELD, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getDbxrefTextFieldValue()).thenReturn("CG42235");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	// Qualifier field
	
	/**
	 * Test case when we have no current "Qualifier" value but user has entered one
	 */
	@Test
	public void testIsQualifierChangedForQualifierField1()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_1, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getQualifierTextFieldValue()).thenReturn("qual-value1");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	/**
	 * Test case when we have a current "Qualifier" value that has not been changed by user
	 */
	@Test
	public void testIsQualifierChangedForQualifierField2()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_QUALIFIER_FIELD, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getQualifierTextFieldValue()).thenReturn("qual-value1");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertFalse(isChanged);
	}
	
	/**
	 * Test case when we have a current "Qualifier" value that has been changed by user
	 */
	@Test
	public void testIsQualifierChangedForQualifierField3()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_QUALIFIER_FIELD, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getQualifierTextFieldValue()).thenReturn("qual-value2");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	// Evidence field
	
	/**
	 * Test case when we have no current "Evidence" value but user has entered one
	 */
	@Test
	public void testIsQualifierChangedForEvidenceField1()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_1, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getEvidenceListSelectedIndex()).thenReturn(2);
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	/**
	 * Test case when we have a current "Evidence" value that has not been changed by user
	 */
	@Test
	public void testIsQualifierChangedForEvidenceField2()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_EVIDENCE_FIELD, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getEvidenceListSelectedIndex()).thenReturn(1);
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertFalse(isChanged);
	}
	
	/**
	 * Test case when we have a current "Evidence" value that has been changed by user
	 */
	@Test
	public void testIsQualifierChangedForEvidenceField3()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_EVIDENCE_FIELD, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getEvidenceListSelectedIndex()).thenReturn(5);
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	// Source field
	
	/**
	 * Test case when we have no current "Source" value but user has entered one
	 */
	@Test
	public void testIsQualifierChangedForSourceField1()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_1, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getSourceListSelectedIndex()).thenReturn(2);
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	/**
	 * Test case when we have no current "Source" value 
	 * and user has default "none" list option.
	 */
	@Test
	public void testIsQualifierChangedForSourceField1a()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_1, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getSourceListSelectedIndex()).thenReturn(0);
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertFalse(isChanged);
	}
	
	/**
	 * Test case when we have a current "Source" value that has not been changed by user
	 */
	@Test
	public void testIsQualifierChangedForSourceField2()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_SOURCE_FIELD, MOLECULAR_FUNCTION_KEY);
		
		int idx = GoBox.getSourceIndex("EuPathDB");
		
		// When
		when (goBox.getSourceListSelectedIndex()).thenReturn(idx);
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertFalse(isChanged);
	}
	
	/**
	 * Test case when we have a current "Source" value that has been changed by user
	 */
	@Test
	public void testIsQualifierChangedForSourceField3()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_SOURCE_FIELD, MOLECULAR_FUNCTION_KEY);
		
		int idx = GoBox.getSourceIndex(GENEDB_SOURCE);
		
		// When
		when (goBox.getSourceListSelectedIndex()).thenReturn(idx);
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	// Date field
	
	/**
	 * Test case when we have no current "Date" value but user has entered one
	 */
	@Test
	public void testIsQualifierChangedForDateField1()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_1, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getDateFieldValue()).thenReturn("20190220");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	/**
	 * Test case when we have a current "Date" value that has not been changed by user
	 */
	@Test
	public void testIsQualifierChangedForDateField2()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_SOURCE_FIELD, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getDateFieldValue()).thenReturn("20190220");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertFalse(isChanged);
	}
	
	/**
	 * Test case when we have a current "Date" value that has been changed by user
	 */
	@Test
	public void testIsQualifierChangedForDateField3()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_WITH_SOURCE_FIELD, MOLECULAR_FUNCTION_KEY);
		
		// When
		when (goBox.getDateFieldValue()).thenReturn("20190222");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	/**
	 * Test case when we have multiple user changes to fields
	 */
	@Test
	public void testIsQualifierChangedForMultiChanges()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_ALL_TERMS, MOLECULAR_FUNCTION_KEY);
		
		int idx = GoBox.getSourceIndex(GENEDB_SOURCE);
		
		// When
		when (goBox.getQualifierTextFieldValue()).thenReturn("newval");
		when (goBox.getSourceListSelectedIndex()).thenReturn(idx);
		when (goBox.getDateFieldValue()).thenReturn("20190225");
		boolean isChanged = goBox.isQualifierChanged();
		
		// Then
		assertTrue(isChanged);
	}
	
	// ==============================================================
	// updateQualifierString tests
	// ==============================================================
	
	/**
	 * Test case when no GO string change has occurred.
	 */
	@Test
	public void testUpdateQualifierStringForNoChange()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_ALL_TERMS, MOLECULAR_FUNCTION_KEY);
		
		int evIdx = GoBox.getEvidenceIndex("Inferred from Sequence Model");
		int srcIdx = GoBox.getSourceIndex("EuPathDB");
		
		// When
		when (goBox.getWithTextFieldValue()).thenReturn("with-value1");
		when (goBox.getQualifierTextFieldValue()).thenReturn("qual-value1");
		when (goBox.getSourceListSelectedIndex()).thenReturn(srcIdx);
		when (goBox.getEvidenceListSelectedIndex()).thenReturn(evIdx);
		when (goBox.getDateFieldValue()).thenReturn("20190220");
		
		String newVersion = goBox.updateQualifierString();
		
		// Then
		assertEquals(GO_TERM_STRING_ALL_TERMS, newVersion);
	}
	
	/**
	 * Test case when we have no values in the GO string and user has
	 * specified all new values on fields.
	 */
	@Test
	public void testUpdateQualifierStringForAllChanged1()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_1, MOLECULAR_FUNCTION_KEY);
		
		int evIdx = GoBox.getEvidenceIndex("Inferred from Sequence Model");
		int srcIdx = GoBox.getSourceIndex("EuPathDB");
		
		// When
		when (goBox.getWithTextFieldValue()).thenReturn("new-with-val1");
		when (goBox.getDbxrefTextFieldValue()).thenReturn("CG42255");
		when (goBox.getQualifierTextFieldValue()).thenReturn("new-qual-value1");
		when (goBox.getSourceListSelectedIndex()).thenReturn(srcIdx);
		when (goBox.getEvidenceListSelectedIndex()).thenReturn(evIdx);
		when (goBox.getDateFieldValue()).thenReturn("20190225");
		
		String newVersion = goBox.updateQualifierString();
		
		// Then
		assertEquals(
			"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;with=new-with-val1;db_xref=CG42255;evidence=Inferred from Sequence Model;qualifier=new-qual-value1;assigned_by=EuPathDB;date=20190225", 
			newVersion);
	}
	
	/**
	 * Test case when we have all values in the GO string and user has
	 * specified to override them all.
	 */
	@Test
	public void testUpdateQualifierStringForAllChanged2()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_ALL_TERMS, MOLECULAR_FUNCTION_KEY);
		
		int evIdx = GoBox.getEvidenceIndex("Inferred from High Throughput Direct Assay");
		int srcIdx = GoBox.getSourceIndex(GENEDB_SOURCE);
		
		// When
		when (goBox.getWithTextFieldValue()).thenReturn("new-with-val1");
		when (goBox.getDbxrefTextFieldValue()).thenReturn("CG42255");
		when (goBox.getQualifierTextFieldValue()).thenReturn("new-qual-value1");
		when (goBox.getSourceListSelectedIndex()).thenReturn(srcIdx);
		when (goBox.getEvidenceListSelectedIndex()).thenReturn(evIdx);
		when (goBox.getDateFieldValue()).thenReturn("20190225");
		
		String newVersion = goBox.updateQualifierString();
		
		// Then
		assertEquals(
			"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;with=new-with-val1;db_xref=CG42255;qualifier=new-qual-value1;evidence=Inferred from High Throughput Direct Assay;assigned_by="+GENEDB_SOURCE+";date=20190225", 
			newVersion);
	}
	
	/**
	 * 
	 */
	@Test
	public void testUpdateQualifierStringForRemovedFields()
	{
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_ALL_TERMS, MOLECULAR_FUNCTION_KEY);
		
		int evIdx = GoBox.getEvidenceIndex("Inferred from High Throughput Direct Assay");
		
		// When
		when (goBox.getWithTextFieldValue()).thenReturn("");
		when (goBox.getDbxrefTextFieldValue()).thenReturn("");
		when (goBox.getQualifierTextFieldValue()).thenReturn("");
		when (goBox.getSourceListSelectedIndex()).thenReturn(0);
		when (goBox.getEvidenceListSelectedIndex()).thenReturn(evIdx);
		when (goBox.getDateFieldValue()).thenReturn("20190225");
		
		String newVersion = goBox.updateQualifierString();
		
		// Then 
		// I don't think ;; should really be appearing but it doesn't do any harm  
		// and this behaviour has not been changed, so leave as is to avoid breaking other code.
		assertEquals(
			"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;;;;evidence=Inferred from High Throughput Direct Assay;;date=20190225", 
			newVersion);
	}
	
	// ==============================================================
	// updateQualifier tests
	// ==============================================================
	
	/**
	 * Test updating of all GO fields.
	 */
	@Test
	public void testUpdateQualifierAllFieldsChanged()
	{
		final int PROD_INDEX = 0;
		final int GO_INDEX = PROD_INDEX + 1;
		
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_ALL_TERMS, MOLECULAR_FUNCTION_KEY);
		
		StringVector productValues = new StringVector();
		productValues.add("term=putative proteasome regulatory particle base subunit RPT3;");
		
		StringVector goValues = new StringVector();
		goValues.add("GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;evidence=Inferred from Sequence Model;date=20190222");
		
		QualifierVector qv = new QualifierVector();
		qv.add(0, new Qualifier("product", productValues));
		qv.add(GO_INDEX, new Qualifier("GO", goValues));
		
		int evIdx = GoBox.getEvidenceIndex("Inferred from High Throughput Direct Assay");
		int srcIdx = GoBox.getSourceIndex(GENEDB_SOURCE);
		
		// When
		when (goBox.getWithTextFieldValue()).thenReturn("new-with-val1");
		when (goBox.getDbxrefTextFieldValue()).thenReturn("CG42255");
		when (goBox.getQualifierTextFieldValue()).thenReturn("new-qual-value1");
		when (goBox.getSourceListSelectedIndex()).thenReturn(srcIdx);
		when (goBox.getEvidenceListSelectedIndex()).thenReturn(evIdx);
		when (goBox.getDateFieldValue()).thenReturn("20190225");
		
		goBox.updateQualifier(qv);
		Qualifier newVersion = qv.get(GO_INDEX);
		
		// Then
		assertEquals(
			"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;with=new-with-val1;db_xref=CG42255;qualifier=new-qual-value1;evidence=Inferred from High Throughput Direct Assay;assigned_by="+GENEDB_SOURCE+";date=20190225", 
			newVersion.getValues().elementAt(0));
	}
	
	/**
	 * Test removal of GO fields.
	 */
	@Test
	public void testUpdateQualifierForFieldRemoval()
	{
		final int PROD_INDEX = 0;
		final int GO_INDEX = PROD_INDEX + 1;
		
		// Given
		GoBox goBox = getGoBoxForTest(GO_TERM_STRING_ALL_TERMS, MOLECULAR_FUNCTION_KEY);
		
		StringVector productValues = new StringVector();
		productValues.add("term=putative proteasome regulatory particle base subunit RPT3;");
		
		StringVector goValues = new StringVector();
		goValues.add("GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;evidence=Inferred from Sequence Model;date=20190222");
		
		QualifierVector qv = new QualifierVector();
		qv.add(0, new Qualifier("product", productValues));
		qv.add(GO_INDEX, new Qualifier("GO", goValues));
		
		int evIdx = GoBox.getEvidenceIndex("Inferred from High Throughput Direct Assay");
		
		// When
		when (goBox.getWithTextFieldValue()).thenReturn("");
		when (goBox.getDbxrefTextFieldValue()).thenReturn("");
		when (goBox.getQualifierTextFieldValue()).thenReturn("");
		when (goBox.getSourceListSelectedIndex()).thenReturn(0);
		when (goBox.getEvidenceListSelectedIndex()).thenReturn(evIdx);
		when (goBox.getDateFieldValue()).thenReturn("20190225");
		
		goBox.updateQualifier(qv);
		Qualifier newVersion = qv.get(GO_INDEX);
		
		// Then
		
		assertTrue(qv.size() == 2);
		
		assertEquals(
			"GOid=GO:0102201;aspect=molecular_function;term=(+)-2-epi-prezizaene synthase activity;;;;evidence=Inferred from High Throughput Direct Assay;;date=20190225", 
			newVersion.getValues().elementAt(0));
	}
	
	
	
}
