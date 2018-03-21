/* 
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
package uk.ac.sanger.artemis.components;


import static org.junit.Assert.*;

import org.junit.Test;

/**
 * JUnit test class for the DbfetchEntrySource class.
 * 
 * @author kp11
 *
 */
public class DbfetchEntrySourceTest 
{
  /**
   * Test the constructEbiUrl method.
   */
  @Test
  public void testConstructEbiUrl() 
  {
	  String accessionNumber = null;
	  
	  DbfetchEntrySource ebifetch = new DbfetchEntrySource(null); 
	  
	  // ==== REF SEQ ====
	  
	  accessionNumber = "NC_017633";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_REF_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_FASTA),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  accessionNumber = "NC_1";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_REF_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_FASTA),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  accessionNumber = "NC_017633.1";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_REF_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_FASTA),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  accessionNumber = "nc_017633.1";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_REF_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_FASTA),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  accessionNumber = "nc_017633567.12222222";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_REF_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_FASTA),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  // ==== ENA SEQ ======
	  
	  accessionNumber = "FN649414";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_ENA_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_EMBL),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  accessionNumber = "F649414";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_ENA_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_EMBL),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  accessionNumber = "FN1";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_ENA_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_EMBL),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  accessionNumber = "FN649414.1";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_ENA_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_EMBL),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  accessionNumber = "fn649414.1";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_ENA_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_EMBL),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  accessionNumber = "fn649433333.13333333";
	  assertEquals(
			  String.format(DbfetchEntrySource.EBI_URL, DbfetchEntrySource.EBI_ENA_SEQ_DB, accessionNumber, DbfetchEntrySource.EBI_FORMAT_EMBL),
			  ebifetch.constructEbiUrl(accessionNumber) 
	  );
	  
	  // ==== Invalids =====
	  
	  assertNull(ebifetch.constructEbiUrl("N_017633"));
	  assertNull(ebifetch.constructEbiUrl("N_017633.123"));
	  assertNull(ebifetch.constructEbiUrl("_017633"));
	  assertNull(ebifetch.constructEbiUrl("649414"));
	  assertNull(ebifetch.constructEbiUrl("DDDDDDDDD649414"));
	  assertNull(ebifetch.constructEbiUrl("NFGFG64941488"));
	  assertNull(ebifetch.constructEbiUrl("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
	  assertNull(ebifetch.constructEbiUrl(""));
	  
  }
}