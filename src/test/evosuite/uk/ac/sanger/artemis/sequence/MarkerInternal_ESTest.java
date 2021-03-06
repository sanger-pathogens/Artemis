/*
 * This file was automatically generated by EvoSuite
 * Thu Sep 20 12:47:09 GMT 2018
 */

package uk.ac.sanger.artemis.sequence;

import org.junit.Test;
import static org.junit.Assert.*;
import static org.evosuite.shaded.org.mockito.Mockito.*;
import static org.evosuite.runtime.EvoAssertions.*;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.evosuite.runtime.ViolatedAssumptionAnswer;
import org.junit.runner.RunWith;
import uk.ac.sanger.artemis.io.Sequence;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.MarkerInternal;
import uk.ac.sanger.artemis.sequence.SequenceChangeEvent;
import uk.ac.sanger.artemis.sequence.Strand;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, useJEE = true) 
public class MarkerInternal_ESTest extends MarkerInternal_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test00()  throws Throwable  {
      Bases bases0 = mock(Bases.class, new ViolatedAssumptionAnswer());
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(1).when(strand0).getSequenceLength();
      MarkerInternal markerInternal0 = new MarkerInternal(strand0, 1);
      SequenceChangeEvent sequenceChangeEvent0 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn((String) null).when(sequenceChangeEvent0).getSubSequence();
      doReturn(1, 102, 0).when(sequenceChangeEvent0).getType();
      // Undeclared exception!
      try { 
        markerInternal0.sequenceChanged(sequenceChangeEvent0);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.sequence.MarkerInternal", e);
      }
  }

  @Test(timeout = 4000)
  public void test01()  throws Throwable  {
      Bases bases0 = mock(Bases.class, new ViolatedAssumptionAnswer());
      doReturn((-97)).when(bases0).getLength();
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn("position: ", "").when(strand0).toString();
      doReturn(bases0).when(strand0).getBases();
      doReturn(1997).when(strand0).getSequenceLength();
      doReturn(false).when(strand0).isForwardStrand();
      MarkerInternal markerInternal0 = new MarkerInternal(strand0, 2);
      SequenceChangeEvent sequenceChangeEvent0 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn(5).when(sequenceChangeEvent0).getPosition();
      doReturn("\"n%_so%9.=w").when(sequenceChangeEvent0).getSubSequence();
      doReturn(2149, (-595), 69, (-5), 1).when(sequenceChangeEvent0).getType();
      markerInternal0.sequenceChanged(sequenceChangeEvent0);
      assertEquals((-9), markerInternal0.getPosition());
  }

  @Test(timeout = 4000)
  public void test02()  throws Throwable  {
      Bases bases0 = mock(Bases.class, new ViolatedAssumptionAnswer());
      doReturn(1).when(bases0).getLength();
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn("vH?6", "vH?6").when(strand0).toString();
      doReturn(bases0).when(strand0).getBases();
      doReturn(1).when(strand0).getSequenceLength();
      doReturn(false).when(strand0).isForwardStrand();
      MarkerInternal markerInternal0 = new MarkerInternal(strand0, 1);
      SequenceChangeEvent sequenceChangeEvent0 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn(2).when(sequenceChangeEvent0).getPosition();
      doReturn("/MBq-Ld,)").when(sequenceChangeEvent0).getSubSequence();
      doReturn(1302, 206, 1, 2, 1).when(sequenceChangeEvent0).getType();
      markerInternal0.sequenceChanged(sequenceChangeEvent0);
      assertEquals((-9), markerInternal0.getPosition());
  }

  @Test(timeout = 4000)
  public void test03()  throws Throwable  {
      Bases bases0 = mock(Bases.class, new ViolatedAssumptionAnswer());
      Sequence sequence0 = mock(Sequence.class, new ViolatedAssumptionAnswer());
      Bases bases1 = mock(Bases.class, new ViolatedAssumptionAnswer());
      Bases bases2 = mock(Bases.class, new ViolatedAssumptionAnswer());
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn("vH?6", "vH?6", "vH?6", "vH?6", (String) null).when(strand0).toString();
      doReturn(bases1, (Bases) null).when(strand0).getBases();
      doReturn(1, 5).when(strand0).getSequenceLength();
      doReturn(true, false).when(strand0).isForwardStrand();
      MarkerInternal markerInternal0 = new MarkerInternal(strand0, 1);
      markerInternal0.getStrand();
      markerInternal0.getStrand();
      SequenceChangeEvent sequenceChangeEvent0 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn((-1913)).when(sequenceChangeEvent0).getPosition();
      doReturn("vH?6").when(sequenceChangeEvent0).getSubSequence();
      doReturn((-3242), 2, (-1870), 1).when(sequenceChangeEvent0).getType();
      markerInternal0.sequenceChanged(sequenceChangeEvent0);
      markerInternal0.setPosition(2);
      SequenceChangeEvent sequenceChangeEvent1 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn((-3959), 1, 5).when(sequenceChangeEvent1).getType();
      markerInternal0.sequenceChanged(sequenceChangeEvent1);
      markerInternal0.getPosition();
      SequenceChangeEvent sequenceChangeEvent2 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn(2, 4).when(sequenceChangeEvent2).getType();
      markerInternal0.sequenceChanged(sequenceChangeEvent2);
      SequenceChangeEvent sequenceChangeEvent3 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn("vH?6").when(sequenceChangeEvent3).getSubSequence();
      doReturn(4165, (-1913), 1, 0).when(sequenceChangeEvent3).getType();
      // Undeclared exception!
      try { 
        markerInternal0.sequenceChanged(sequenceChangeEvent3);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.sequence.MarkerInternal", e);
      }
  }

  @Test(timeout = 4000)
  public void test04()  throws Throwable  {
      Bases bases0 = mock(Bases.class, new ViolatedAssumptionAnswer());
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn("vH?6", "vH?6").when(strand0).toString();
      doReturn(bases0, (Bases) null).when(strand0).getBases();
      doReturn(1).when(strand0).getSequenceLength();
      doReturn(true).when(strand0).isForwardStrand();
      MarkerInternal markerInternal0 = new MarkerInternal(strand0, 1);
      SequenceChangeEvent sequenceChangeEvent0 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn((-1913)).when(sequenceChangeEvent0).getPosition();
      doReturn("vH?6").when(sequenceChangeEvent0).getSubSequence();
      doReturn((-3242), 102, (-1870), 1).when(sequenceChangeEvent0).getType();
      markerInternal0.sequenceChanged(sequenceChangeEvent0);
      // Undeclared exception!
      try { 
        markerInternal0.setPosition(0);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.sequence.MarkerInternal", e);
      }
  }

  @Test(timeout = 4000)
  public void test05()  throws Throwable  {
      Bases bases0 = mock(Bases.class, new ViolatedAssumptionAnswer());
      doReturn((-5)).when(bases0).getLength();
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn("4O:3G:", "4O:3G:", "u4.ac.sanger.artemBs.Feature", "4O:3G:").when(strand0).toString();
      doReturn(bases0, (Bases) null).when(strand0).getBases();
      doReturn(1).when(strand0).getSequenceLength();
      doReturn(false, false).when(strand0).isForwardStrand();
      MarkerInternal markerInternal0 = new MarkerInternal(strand0, 1);
      SequenceChangeEvent sequenceChangeEvent0 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn(0).when(sequenceChangeEvent0).getPosition();
      doReturn("4O:3G:").when(sequenceChangeEvent0).getSubSequence();
      doReturn((-3378), (-28), (-3031), 1800, 1).when(sequenceChangeEvent0).getType();
      markerInternal0.sequenceChanged(sequenceChangeEvent0);
      SequenceChangeEvent sequenceChangeEvent1 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn(" > ").when(sequenceChangeEvent1).getSubSequence();
      doReturn(1, 113, 2, 0).when(sequenceChangeEvent1).getType();
      // Undeclared exception!
      try { 
        markerInternal0.sequenceChanged(sequenceChangeEvent1);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.sequence.MarkerInternal", e);
      }
  }

  @Test(timeout = 4000)
  public void test06()  throws Throwable  {
      Bases bases0 = mock(Bases.class, new ViolatedAssumptionAnswer());
      doReturn((Sequence) null).when(bases0).getSequence();
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases0).when(strand0).getBases();
      doReturn((-1275)).when(strand0).getSequenceLength();
      MarkerInternal markerInternal0 = null;
      try {
        markerInternal0 = new MarkerInternal(strand0, 2);
        fail("Expecting exception: Exception");
      
      } catch(Throwable e) {
         //
         // position: 2
         //
         verifyException("uk.ac.sanger.artemis.sequence.MarkerInternal", e);
      }
  }

  @Test(timeout = 4000)
  public void test07()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(0).when(strand0).getSequenceLength();
      MarkerInternal markerInternal0 = null;
      try {
        markerInternal0 = new MarkerInternal(strand0, 1);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.sequence.MarkerInternal", e);
      }
  }

  @Test(timeout = 4000)
  public void test08()  throws Throwable  {
      Bases bases0 = mock(Bases.class, new ViolatedAssumptionAnswer());
      doReturn((-5)).when(bases0).getLength();
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn("4O:3G:", "4O:3G:").when(strand0).toString();
      doReturn(bases0).when(strand0).getBases();
      doReturn(1).when(strand0).getSequenceLength();
      doReturn(false).when(strand0).isForwardStrand();
      MarkerInternal markerInternal0 = new MarkerInternal(strand0, 1);
      SequenceChangeEvent sequenceChangeEvent0 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn(0).when(sequenceChangeEvent0).getPosition();
      doReturn(" > ").when(sequenceChangeEvent0).getSubSequence();
      doReturn(1, 113, 2, 0, 0).when(sequenceChangeEvent0).getType();
      markerInternal0.sequenceChanged(sequenceChangeEvent0);
      assertEquals(4, markerInternal0.getPosition());
  }

  @Test(timeout = 4000)
  public void test09()  throws Throwable  {
      Bases bases0 = mock(Bases.class, new ViolatedAssumptionAnswer());
      doReturn(0).when(bases0).getLength();
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn("vH?6", "vH?6").when(strand0).toString();
      doReturn(bases0).when(strand0).getBases();
      doReturn(1).when(strand0).getSequenceLength();
      doReturn(false).when(strand0).isForwardStrand();
      MarkerInternal markerInternal0 = new MarkerInternal(strand0, 1);
      SequenceChangeEvent sequenceChangeEvent0 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn((-5)).when(sequenceChangeEvent0).getPosition();
      doReturn("vH?6").when(sequenceChangeEvent0).getSubSequence();
      doReturn((-3242), 102, (-1870), 1, 0).when(sequenceChangeEvent0).getType();
      markerInternal0.sequenceChanged(sequenceChangeEvent0);
      assertEquals(1, markerInternal0.getPosition());
  }

  @Test(timeout = 4000)
  public void test10()  throws Throwable  {
      Bases bases0 = mock(Bases.class, new ViolatedAssumptionAnswer());
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(5, 0).when(strand0).getSequenceLength();
      MarkerInternal markerInternal0 = new MarkerInternal(strand0, 2);
      // Undeclared exception!
      try { 
        markerInternal0.setPosition(2);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.sequence.MarkerInternal", e);
      }
  }

  @Test(timeout = 4000)
  public void test11()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(1).when(strand0).getSequenceLength();
      MarkerInternal markerInternal0 = new MarkerInternal(strand0, 1);
      SequenceChangeEvent sequenceChangeEvent0 = mock(SequenceChangeEvent.class, new ViolatedAssumptionAnswer());
      doReturn(3).when(sequenceChangeEvent0).getType();
      markerInternal0.sequenceChanged(sequenceChangeEvent0);
      assertEquals(1, markerInternal0.getPosition());
  }
}
