/*
 * This file was automatically generated by EvoSuite
 * Wed Sep 19 21:21:34 GMT 2018
 */

package uk.ac.sanger.artemis.plot;

import org.junit.Test;
import static org.junit.Assert.*;
import static org.evosuite.shaded.org.mockito.Mockito.*;
import static org.evosuite.runtime.EvoAssertions.*;
import java.awt.Color;
import java.awt.Graphics;
import java.io.StringReader;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.evosuite.runtime.ViolatedAssumptionAnswer;
import org.junit.runner.RunWith;
import uk.ac.sanger.artemis.io.FastaStreamSequence;
import uk.ac.sanger.artemis.io.GenbankStreamSequence;
import uk.ac.sanger.artemis.io.RawStreamSequence;
import uk.ac.sanger.artemis.io.Sequence;
import uk.ac.sanger.artemis.plot.ICDIAlgorithm;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.LinePushBackReader;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, separateClassLoader = false, useJEE = true) 
public class ICDIAlgorithm_ESTest extends ICDIAlgorithm_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test00()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false, false).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      iCDIAlgorithm0.getDefaultStepSize(10);
  }

  @Test(timeout = 4000)
  public void test01()  throws Throwable  {
      Bases bases0 = new Bases((Sequence) null);
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases0).when(strand0).getBases();
      doReturn(false, true).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      float[] floatArray0 = new float[0];
      // Undeclared exception!
      try { 
        iCDIAlgorithm0.getValues(0, (-149), floatArray0);
        fail("Expecting exception: Error");
      
      } catch(Error e) {
         //
         // internal error - unexpected exception: org.evosuite.runtime.mock.java.lang.MockThrowable: start: 0 > end: -148
         //
         verifyException("uk.ac.sanger.artemis.plot.ICDIAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test02()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false, false).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      float[] floatArray0 = new float[1];
      // Undeclared exception!
      try { 
        iCDIAlgorithm0.getValues(24, 2, floatArray0);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test03()  throws Throwable  {
      GenbankStreamSequence genbankStreamSequence0 = new GenbankStreamSequence("wwxIa8\"HE");
      Bases bases0 = new Bases(genbankStreamSequence0);
      Strand strand0 = bases0.getForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      float[] floatArray0 = new float[0];
      // Undeclared exception!
      try { 
        iCDIAlgorithm0.getValues(0, 1727, floatArray0);
        fail("Expecting exception: ArrayIndexOutOfBoundsException");
      
      } catch(ArrayIndexOutOfBoundsException e) {
         //
         // 0
         //
         verifyException("uk.ac.sanger.artemis.plot.ICDIAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test04()  throws Throwable  {
      ICDIAlgorithm iCDIAlgorithm0 = null;
      try {
        iCDIAlgorithm0 = new ICDIAlgorithm((Strand) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.ICDIAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test05()  throws Throwable  {
      StringReader stringReader0 = new StringReader("");
      LinePushBackReader linePushBackReader0 = new LinePushBackReader(stringReader0);
      RawStreamSequence rawStreamSequence0 = new RawStreamSequence(linePushBackReader0);
      Bases bases0 = new Bases(rawStreamSequence0);
      Strand strand0 = bases0.getReverseStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      iCDIAlgorithm0.getDefaultStepSize(28);
  }

  @Test(timeout = 4000)
  public void test06()  throws Throwable  {
      Bases bases0 = new Bases((Sequence) null);
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases0).when(strand0).getBases();
      doReturn(false, true).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      iCDIAlgorithm0.getDefaultStepSize(1);
  }

  @Test(timeout = 4000)
  public void test07()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false, false).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      iCDIAlgorithm0.getDefaultMinWindowSize();
  }

  @Test(timeout = 4000)
  public void test08()  throws Throwable  {
      GenbankStreamSequence genbankStreamSequence0 = new GenbankStreamSequence("wwxIa8\"HE");
      Bases bases0 = new Bases(genbankStreamSequence0);
      Strand strand0 = bases0.getForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      iCDIAlgorithm0.getDefaultMaxWindowSize();
  }

  @Test(timeout = 4000)
  public void test09()  throws Throwable  {
      FastaStreamSequence fastaStreamSequence0 = new FastaStreamSequence("invisible_qualifiers_gff");
      Bases bases0 = new Bases(fastaStreamSequence0);
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases0).when(strand0).getBases();
      doReturn(false, false).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      iCDIAlgorithm0.getDefaultWindowSize();
  }

  @Test(timeout = 4000)
  public void test10()  throws Throwable  {
      StringReader stringReader0 = new StringReader("");
      LinePushBackReader linePushBackReader0 = new LinePushBackReader(stringReader0);
      RawStreamSequence rawStreamSequence0 = new RawStreamSequence(linePushBackReader0);
      Bases bases0 = new Bases(rawStreamSequence0);
      Strand strand0 = bases0.getReverseStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      iCDIAlgorithm0.setRevCompDisplay(true);
      float[] floatArray0 = new float[2];
      // Undeclared exception!
      try { 
        iCDIAlgorithm0.getValues(656, (-2336), floatArray0);
        fail("Expecting exception: Error");
      
      } catch(Error e) {
         //
         // internal error - unexpected exception: org.evosuite.runtime.mock.java.lang.MockThrowable: start: 656 > end: -2336
         //
         verifyException("uk.ac.sanger.artemis.plot.ICDIAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test11()  throws Throwable  {
      FastaStreamSequence fastaStreamSequence0 = new FastaStreamSequence("invisible_qualifiers_gff");
      Bases bases0 = new Bases(fastaStreamSequence0);
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases0).when(strand0).getBases();
      doReturn(false, false).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      float[] floatArray0 = new float[7];
      iCDIAlgorithm0.getValues((-2489), 97, floatArray0);
      assertTrue(iCDIAlgorithm0.scalingFlag());
  }

  @Test(timeout = 4000)
  public void test12()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false, false).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      iCDIAlgorithm0.getAverage();
      assertEquals("Reverse Intrinsic Codon Deviation Index", iCDIAlgorithm0.getAlgorithmName());
      assertTrue(iCDIAlgorithm0.scalingFlag());
  }

  @Test(timeout = 4000)
  public void test13()  throws Throwable  {
      FastaStreamSequence fastaStreamSequence0 = new FastaStreamSequence("invisible_qualifiers_gff");
      Bases bases0 = new Bases(fastaStreamSequence0);
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases0).when(strand0).getBases();
      doReturn(false, false).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      Float float0 = iCDIAlgorithm0.getMinimumInternal();
      assertTrue(iCDIAlgorithm0.scalingFlag());
      assertEquals("Reverse Intrinsic Codon Deviation Index", iCDIAlgorithm0.getAlgorithmName());
      assertEquals(0.0F, (float)float0, 0.01F);
  }

  @Test(timeout = 4000)
  public void test14()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false, false).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      Float float0 = iCDIAlgorithm0.getMaximumInternal();
      assertEquals("Reverse Intrinsic Codon Deviation Index", iCDIAlgorithm0.getAlgorithmName());
      assertTrue(iCDIAlgorithm0.scalingFlag());
      assertEquals(1000.0F, (float)float0, 0.01F);
  }

  @Test(timeout = 4000)
  public void test15()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false, false).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      Graphics graphics0 = mock(Graphics.class, new ViolatedAssumptionAnswer());
      Color[] colorArray0 = new Color[6];
      iCDIAlgorithm0.drawLegend(graphics0, 0, 0, colorArray0);
      assertEquals("Reverse Intrinsic Codon Deviation Index", iCDIAlgorithm0.getAlgorithmName());
      assertTrue(iCDIAlgorithm0.scalingFlag());
  }

  @Test(timeout = 4000)
  public void test16()  throws Throwable  {
      FastaStreamSequence fastaStreamSequence0 = new FastaStreamSequence("invisible_qualifiers_gff");
      Bases bases0 = new Bases(fastaStreamSequence0);
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases0).when(strand0).getBases();
      doReturn(false, false).when(strand0).isForwardStrand();
      ICDIAlgorithm iCDIAlgorithm0 = new ICDIAlgorithm(strand0);
      int int0 = iCDIAlgorithm0.getValueCount();
      assertEquals(3, int0);
      assertEquals("Reverse Intrinsic Codon Deviation Index", iCDIAlgorithm0.getAlgorithmName());
      assertTrue(iCDIAlgorithm0.scalingFlag());
  }
}