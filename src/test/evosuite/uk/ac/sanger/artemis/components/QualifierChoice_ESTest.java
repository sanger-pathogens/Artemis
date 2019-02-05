/*
 * This file was automatically generated by EvoSuite
 * Wed Sep 19 21:40:25 GMT 2018
 */

package uk.ac.sanger.artemis.components;

import org.junit.Test;
import static org.junit.Assert.*;
import static org.evosuite.shaded.org.mockito.Mockito.*;
import static org.evosuite.runtime.EvoAssertions.*;
import java.awt.HeadlessException;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.evosuite.runtime.ViolatedAssumptionAnswer;
import org.junit.runner.RunWith;
import uk.ac.sanger.artemis.components.QualifierChoice;
import uk.ac.sanger.artemis.io.BlastEntryInformation;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.MSPcrunchEntryInformation;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, separateClassLoader = false, useJEE = true) 
public class QualifierChoice_ESTest extends QualifierChoice_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test0()  throws Throwable  {
      QualifierChoice qualifierChoice0 = null;
      try {
        qualifierChoice0 = new QualifierChoice((EntryInformation) null, (Key) null, "QqCgb^]N'*#y=e76{}", false);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.components.QualifierChoice", e);
      }
  }

  @Test(timeout = 4000)
  public void test1()  throws Throwable  {
      MSPcrunchEntryInformation mSPcrunchEntryInformation0 = new MSPcrunchEntryInformation();
      Key key0 = Key.CDS;
      QualifierChoice qualifierChoice0 = null;
      try {
        qualifierChoice0 = new QualifierChoice(mSPcrunchEntryInformation0, key0, "locus_ag", true);
      
      } catch(Exception e) {
    	  	fail("Caught: Exception: " + e.getMessage());
      }
  }

  @Test(timeout = 4000)
  public void test2()  throws Throwable  {
      EntryInformation entryInformation0 = mock(EntryInformation.class, new ViolatedAssumptionAnswer());
      doReturn(false).when(entryInformation0).isValidQualifier(any(uk.ac.sanger.artemis.io.Key.class) , anyString());
      Key key0 = mock(Key.class, new ViolatedAssumptionAnswer());
      QualifierChoice qualifierChoice0 = null;
      try {
        qualifierChoice0 = new QualifierChoice(entryInformation0, key0, "9OC", false);
        fail("Expecting exception: HeadlessException");
      
      } catch(HeadlessException e) {
      }
  }

  @Test(timeout = 4000)
  public void test3()  throws Throwable  {
      BlastEntryInformation blastEntryInformation0 = new BlastEntryInformation();
      Key key0 = blastEntryInformation0.getDefaultKey();
      QualifierChoice qualifierChoice0 = null;
      try {
        qualifierChoice0 = new QualifierChoice(blastEntryInformation0, key0, (String) null, false);
      
      } catch(Exception e) {
    	  	fail("Caught: Exception: " + e.getMessage());
      }
  }
}