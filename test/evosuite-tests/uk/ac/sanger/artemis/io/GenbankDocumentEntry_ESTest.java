/*
 * This file was automatically generated by EvoSuite
 * Fri Jan 12 14:53:36 GMT 2018
 */

package uk.ac.sanger.artemis.io;

import org.junit.Test;
import static org.junit.Assert.*;
import static org.evosuite.runtime.EvoAssertions.*;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import javax.swing.JPasswordField;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.evosuite.runtime.mock.java.io.MockFile;
import org.evosuite.runtime.mock.java.net.MockURL;
import org.junit.runner.RunWith;
import uk.ac.sanger.artemis.FilteredEntryGroup;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.LogReadListener;
import uk.ac.sanger.artemis.io.BlastEntryInformation;
import uk.ac.sanger.artemis.io.Entry;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.GenbankDocumentEntry;
import uk.ac.sanger.artemis.io.GenbankStreamSequence;
import uk.ac.sanger.artemis.io.MSPcrunchDocumentEntry;
import uk.ac.sanger.artemis.io.MSPcrunchEntryInformation;
import uk.ac.sanger.artemis.io.PublicDBDocumentEntry;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.URLDocument;
import uk.ac.sanger.artemis.util.ZipFileDocument;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, separateClassLoader = true, useJEE = true) 
public class GenbankDocumentEntry_ESTest extends GenbankDocumentEntry_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test00()  throws Throwable  {
      MSPcrunchEntryInformation mSPcrunchEntryInformation0 = new MSPcrunchEntryInformation();
      LogReadListener logReadListener0 = new LogReadListener("v");
      GenbankDocumentEntry genbankDocumentEntry0 = null;
      try {
        genbankDocumentEntry0 = new GenbankDocumentEntry(mSPcrunchEntryInformation0, (Document) null, logReadListener0);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.io.SimpleDocumentEntry", e);
      }
  }

  @Test(timeout = 4000)
  public void test01()  throws Throwable  {
      EntryInformation entryInformation0 = SimpleEntryInformation.getDefaultEntryInformation();
      URL uRL0 = MockURL.getHttpExample();
      URLDocument uRLDocument0 = new URLDocument(uRL0);
      LogReadListener logReadListener0 = new LogReadListener("");
      GenbankDocumentEntry genbankDocumentEntry0 = null;
      try {
        genbankDocumentEntry0 = new GenbankDocumentEntry(entryInformation0, uRLDocument0, logReadListener0);
        fail("Expecting exception: IOException");
      
      } catch(Throwable e) {
         //
         // Could not find: www.someFakeButWellFormedURL.org
         //
         verifyException("org.evosuite.runtime.mock.java.net.EvoHttpURLConnection", e);
      }
  }

  @Test(timeout = 4000)
  public void test02()  throws Throwable  {
      SimpleEntryInformation simpleEntryInformation0 = new SimpleEntryInformation();
      MockFile mockFile0 = new MockFile("");
      ZipFileDocument zipFileDocument0 = new ZipFileDocument(mockFile0, "");
      LogReadListener logReadListener0 = new LogReadListener("");
      GenbankDocumentEntry genbankDocumentEntry0 = null;
      try {
        genbankDocumentEntry0 = new GenbankDocumentEntry(simpleEntryInformation0, zipFileDocument0, logReadListener0);
        fail("Expecting exception: FileNotFoundException");
      
      } catch(Throwable e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("org.evosuite.runtime.mock.java.io.MockFileInputStream", e);
      }
  }

  @Test(timeout = 4000)
  public void test03()  throws Throwable  {
      SimpleEntryGroup simpleEntryGroup0 = new SimpleEntryGroup();
      uk.ac.sanger.artemis.Entry entry0 = simpleEntryGroup0.createEntry("");
      EntryInformation entryInformation0 = entry0.getEntryInformation();
      MSPcrunchDocumentEntry mSPcrunchDocumentEntry0 = new MSPcrunchDocumentEntry(entryInformation0);
      PublicDBDocumentEntry publicDBDocumentEntry0 = new PublicDBDocumentEntry(entryInformation0, mSPcrunchDocumentEntry0, false);
      GenbankDocumentEntry genbankDocumentEntry0 = new GenbankDocumentEntry(entryInformation0, mSPcrunchDocumentEntry0, true);
  }

  @Test(timeout = 4000)
  public void test04()  throws Throwable  {
      EntryInformation entryInformation0 = SimpleEntryInformation.getDefaultEntryInformation();
      GenbankDocumentEntry genbankDocumentEntry0 = null;
      try {
        genbankDocumentEntry0 = new GenbankDocumentEntry(entryInformation0, (Entry) null, true);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.io.SimpleDocumentEntry", e);
      }
  }

  @Test(timeout = 4000)
  public void test05()  throws Throwable  {
      SimpleEntryGroup simpleEntryGroup0 = new SimpleEntryGroup();
      uk.ac.sanger.artemis.Entry entry0 = simpleEntryGroup0.createEntry(";V-oi+Twl~~1");
      GenbankStreamSequence genbankStreamSequence0 = new GenbankStreamSequence("-12-");
      Bases bases0 = new Bases(genbankStreamSequence0);
      Range range0 = new Range(0);
      uk.ac.sanger.artemis.Entry entry1 = entry0.truncate(bases0, range0);
      EntryInformation entryInformation0 = entry1.getEntryInformation();
      GenbankDocumentEntry genbankDocumentEntry0 = new GenbankDocumentEntry(entryInformation0);
  }

  @Test(timeout = 4000)
  public void test06()  throws Throwable  {
      GenbankDocumentEntry genbankDocumentEntry0 = null;
      try {
        genbankDocumentEntry0 = new GenbankDocumentEntry((EntryInformation) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.io.SimpleEntryInformation", e);
      }
  }

  @Test(timeout = 4000)
  public void test07()  throws Throwable  {
      MSPcrunchEntryInformation mSPcrunchEntryInformation0 = new MSPcrunchEntryInformation();
      MSPcrunchDocumentEntry mSPcrunchDocumentEntry0 = new MSPcrunchDocumentEntry(mSPcrunchEntryInformation0);
      GenbankDocumentEntry genbankDocumentEntry0 = new GenbankDocumentEntry(mSPcrunchDocumentEntry0);
      assertTrue(PublicDBDocumentEntry.IGNORE_OBSOLETE_FEATURES);
  }

  @Test(timeout = 4000)
  public void test08()  throws Throwable  {
      SimpleEntryGroup simpleEntryGroup0 = new SimpleEntryGroup();
      uk.ac.sanger.artemis.FeatureVector featureVector0 = new uk.ac.sanger.artemis.FeatureVector();
      FilteredEntryGroup filteredEntryGroup0 = new FilteredEntryGroup(simpleEntryGroup0, featureVector0, " ");
      uk.ac.sanger.artemis.Entry entry0 = filteredEntryGroup0.createEntry("p}");
      EntryInformation entryInformation0 = entry0.getEntryInformation();
      GenbankDocumentEntry genbankDocumentEntry0 = new GenbankDocumentEntry(entryInformation0);
      JPasswordField jPasswordField0 = new JPasswordField();
      DatabaseDocument databaseDocument0 = new DatabaseDocument("PRODUCT STORED AS A CV (product_cv=yes) IN ", jPasswordField0, "", "", true);
      LogReadListener logReadListener0 = new LogReadListener("CDS");
      GenbankDocumentEntry genbankDocumentEntry1 = new GenbankDocumentEntry(entryInformation0, databaseDocument0, logReadListener0);
  }

  @Test(timeout = 4000)
  public void test09()  throws Throwable  {
      GenbankDocumentEntry genbankDocumentEntry0 = null;
      try {
        genbankDocumentEntry0 = new GenbankDocumentEntry((Entry) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.io.GenbankDocumentEntry", e);
      }
  }

  @Test(timeout = 4000)
  public void test10()  throws Throwable  {
      BlastEntryInformation blastEntryInformation0 = new BlastEntryInformation();
      GenbankDocumentEntry genbankDocumentEntry0 = new GenbankDocumentEntry(blastEntryInformation0);
      GenbankDocumentEntry genbankDocumentEntry1 = new GenbankDocumentEntry(blastEntryInformation0, genbankDocumentEntry0, false);
      assertFalse(genbankDocumentEntry1.isReadOnly());
  }
}