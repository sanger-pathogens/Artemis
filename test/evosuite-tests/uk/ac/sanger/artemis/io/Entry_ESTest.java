/*
 * This file was automatically generated by EvoSuite
 * Fri Jan 12 14:42:27 GMT 2018
 */

package uk.ac.sanger.artemis.io;

import org.junit.Test;
import static org.junit.Assert.*;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.junit.runner.RunWith;
import uk.ac.sanger.artemis.io.EmblDocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.PublicDBDocumentEntry;
import uk.ac.sanger.artemis.io.Sequence;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, separateClassLoader = true, useJEE = true) 
public class Entry_ESTest extends Entry_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test0()  throws Throwable  {
      SimpleEntryInformation simpleEntryInformation0 = new SimpleEntryInformation();
      PublicDBDocumentEntry publicDBDocumentEntry0 = new PublicDBDocumentEntry(simpleEntryInformation0);
      Sequence sequence0 = publicDBDocumentEntry0.getSequence();
      assertNull(sequence0);
  }

  @Test(timeout = 4000)
  public void test1()  throws Throwable  {
      SimpleEntryInformation simpleEntryInformation0 = new SimpleEntryInformation();
      EmblDocumentEntry emblDocumentEntry0 = new EmblDocumentEntry(simpleEntryInformation0);
      EntryInformation entryInformation0 = emblDocumentEntry0.getEntryInformation();
      assertNotSame(entryInformation0, simpleEntryInformation0);
  }
}