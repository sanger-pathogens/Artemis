/*
 * This file was automatically generated by EvoSuite
 * Fri Jan 12 16:20:13 GMT 2018
 */

package uk.ac.sanger.artemis;

import org.junit.Test;
import static org.junit.Assert.*;
import static org.evosuite.shaded.org.mockito.Mockito.*;
import static org.evosuite.runtime.MockitoExtension.*;
import static org.evosuite.runtime.EvoAssertions.*;
import java.io.InputStream;
import java.io.PipedInputStream;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.evosuite.runtime.ViolatedAssumptionAnswer;
import org.junit.runner.RunWith;
import uk.ac.sanger.artemis.ProcessMonitor;
import uk.ac.sanger.artemis.components.LogViewer;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, separateClassLoader = true, useJEE = true) 
public class ProcessMonitor_ESTest extends ProcessMonitor_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test0()  throws Throwable  {
      Process process0 = mock(Process.class, new ViolatedAssumptionAnswer());
      doReturn((InputStream) null).when(process0).getErrorStream();
      LogViewer logViewer0 = new LogViewer();
      ProcessMonitor processMonitor0 = new ProcessMonitor(process0, "yO", logViewer0);
      // Undeclared exception!
      try { 
        processMonitor0.run();
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("java.io.Reader", e);
      }
  }

  @Test(timeout = 4000)
  public void test1()  throws Throwable  {
      PipedInputStream pipedInputStream0 = new PipedInputStream();
      Process process0 = mock(Process.class, new ViolatedAssumptionAnswer());
      doReturn((InputStream) null).when(process0).getErrorStream();
      LogViewer logViewer0 = new LogViewer();
      ProcessMonitor processMonitor0 = new ProcessMonitor(process0, "", logViewer0);
      // Undeclared exception!
      try { 
        processMonitor0.run();
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("java.io.Reader", e);
      }
  }
}