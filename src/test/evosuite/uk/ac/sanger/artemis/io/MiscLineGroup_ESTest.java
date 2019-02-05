/*
 * This file was automatically generated by EvoSuite
 * Thu Sep 20 13:50:50 GMT 2018
 */

package uk.ac.sanger.artemis.io;

import org.junit.Test;
import static org.junit.Assert.*;
import static org.evosuite.shaded.org.mockito.Mockito.*;
import static org.evosuite.runtime.EvoAssertions.*;
import java.io.Writer;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.evosuite.runtime.ViolatedAssumptionAnswer;
import org.junit.runner.RunWith;
import uk.ac.sanger.artemis.io.EmblMisc;
import uk.ac.sanger.artemis.io.FeatureHeader;
import uk.ac.sanger.artemis.io.GFFMisc;
import uk.ac.sanger.artemis.util.LinePushBackReader;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, useJEE = true) 
public class MiscLineGroup_ESTest extends MiscLineGroup_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test0()  throws Throwable  {
      LinePushBackReader linePushBackReader0 = mock(LinePushBackReader.class, new ViolatedAssumptionAnswer());
      doReturn((String) null).when(linePushBackReader0).readLine();
      GFFMisc gFFMisc0 = new GFFMisc(linePushBackReader0);
      String string0 = gFFMisc0.getLine();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test1()  throws Throwable  {
      LinePushBackReader linePushBackReader0 = mock(LinePushBackReader.class, new ViolatedAssumptionAnswer());
      doReturn((String) null).when(linePushBackReader0).readLine();
      FeatureHeader featureHeader0 = new FeatureHeader(linePushBackReader0);
      String string0 = featureHeader0.getLine();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test2()  throws Throwable  {
      EmblMisc emblMisc0 = new EmblMisc("*{U]O!ip/m2E#");
      // Undeclared exception!
      try { 
        emblMisc0.writeToStream((Writer) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.io.MiscLineGroup", e);
      }
  }

  @Test(timeout = 4000)
  public void test3()  throws Throwable  {
      EmblMisc emblMisc0 = new EmblMisc("d0q@fnF:m$j(%Nf)");
      String string0 = emblMisc0.getLine();
      assertNotNull(string0);
  }

  @Test(timeout = 4000)
  public void test4()  throws Throwable  {
      EmblMisc emblMisc0 = new EmblMisc("*{U]O!ip/m2E#");
      String string0 = emblMisc0.toString();
      assertEquals("*{U]O!ip/m2E#\n", string0);
  }

  @Test(timeout = 4000)
  public void test5()  throws Throwable  {
      LinePushBackReader linePushBackReader0 = mock(LinePushBackReader.class, new ViolatedAssumptionAnswer());
      doReturn("").when(linePushBackReader0).readLine();
      EmblMisc emblMisc0 = new EmblMisc(linePushBackReader0);
      emblMisc0.setLine("");
      assertEquals("", emblMisc0.getLine());
  }

  @Test(timeout = 4000)
  public void test6()  throws Throwable  {
      LinePushBackReader linePushBackReader0 = mock(LinePushBackReader.class, new ViolatedAssumptionAnswer());
      doReturn("").when(linePushBackReader0).readLine();
      EmblMisc emblMisc0 = new EmblMisc(linePushBackReader0);
      Writer writer0 = mock(Writer.class, new ViolatedAssumptionAnswer());
      emblMisc0.writeToStream(writer0);
      assertEquals("", emblMisc0.getLine());
  }
}