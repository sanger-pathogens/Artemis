/*
 * This file was automatically generated by EvoSuite
 * Wed Sep 19 21:54:51 GMT 2018
 */

package uk.ac.sanger.artemis.plot;

import org.junit.Test;
import static org.junit.Assert.*;
import static org.evosuite.shaded.org.mockito.Mockito.*;
import static org.evosuite.runtime.EvoAssertions.*;
import java.awt.Graphics;
import java.awt.HeadlessException;
import javax.swing.DebugGraphics;
import javax.swing.JPasswordField;
import org.apache.batik.gvt.text.GVTAttributedCharacterIterator;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.evosuite.runtime.ViolatedAssumptionAnswer;
import org.evosuite.runtime.testdata.EvoSuiteFile;
import org.evosuite.runtime.testdata.FileSystemHandling;
import org.junit.runner.RunWith;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.FilteredEntryGroup;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.io.DatabaseStreamFeature;
import uk.ac.sanger.artemis.io.EmblStreamSequence;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.FastaStreamSequence;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.GenbankStreamSequence;
import uk.ac.sanger.artemis.io.PartialSequence;
import uk.ac.sanger.artemis.io.RawStreamSequence;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;
import uk.ac.sanger.artemis.plot.AGWindowAlgorithm;
import uk.ac.sanger.artemis.plot.GCDeviationAlgorithm;
import uk.ac.sanger.artemis.plot.LineAttributes;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.Marker;
import uk.ac.sanger.artemis.sequence.MarkerChangeEvent;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.DatabaseDocument;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, separateClassLoader = false, useJEE = true) 
public class GCDeviationAlgorithm_ESTest extends GCDeviationAlgorithm_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test00()  throws Throwable  {
      SimpleEntryInformation simpleEntryInformation0 = null;
      try {
        simpleEntryInformation0 = new SimpleEntryInformation((EntryInformation) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.io.SimpleEntryInformation", e);
      }
  }

  @Test(timeout = 4000)
  public void test01()  throws Throwable  {
      RawStreamSequence rawStreamSequence0 = new RawStreamSequence("");
      Bases bases0 = new Bases(rawStreamSequence0);
      SimpleEntryGroup simpleEntryGroup0 = new SimpleEntryGroup(bases0);
      Bases bases1 = simpleEntryGroup0.getBases();
      Strand strand0 = bases1.getForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      rawStreamSequence0.setFromChar(bases0.letter_index);
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getAverage();
      gCDeviationAlgorithm0.getMinimumInternal();
      gCDeviationAlgorithm0.getValueCount();
      float[] floatArray0 = new float[9];
      floatArray0[0] = 1894.4F;
      floatArray0[1] = (float) 5;
      floatArray0[0] = (float) 2;
      floatArray0[3] = 2840.7537F;
      floatArray0[4] = (float) (-102);
      floatArray0[5] = (float) 1;
      floatArray0[6] = (float) 5;
      floatArray0[7] = (float) (-102);
      floatArray0[8] = (float) (-5);
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getValues((-102), (-868), floatArray0);
        fail("Expecting exception: Error");
      
      } catch(Error e) {
         //
         // internal error - unexpected exception: org.evosuite.runtime.mock.java.lang.MockThrowable: start: -102 > end: -868
         //
         verifyException("uk.ac.sanger.artemis.plot.GCDeviationAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test02()  throws Throwable  {
      EmblStreamSequence emblStreamSequence0 = new EmblStreamSequence("GC Deviation (G-C)/(G+C)");
      Bases bases0 = new Bases(emblStreamSequence0);
      Strand strand0 = bases0.getReverseStrand();
      AGWindowAlgorithm aGWindowAlgorithm0 = new AGWindowAlgorithm(strand0);
      Strand strand1 = aGWindowAlgorithm0.getStrand();
      Bases bases1 = strand1.getBases();
      Strand strand2 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases1).when(strand2).getBases();
      doReturn(true).when(strand2).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand2);
      float[] floatArray0 = new float[0];
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getValues(1, 2, floatArray0);
        fail("Expecting exception: ArrayIndexOutOfBoundsException");
      
      } catch(ArrayIndexOutOfBoundsException e) {
         //
         // 0
         //
         verifyException("uk.ac.sanger.artemis.plot.GCDeviationAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test03()  throws Throwable  {
      FileSystemHandling.shouldThrowIOException((EvoSuiteFile) null);
      FastaStreamSequence fastaStreamSequence0 = new FastaStreamSequence("");
      Bases bases0 = new Bases(fastaStreamSequence0);
      Strand strand0 = bases0.getReverseStrand();
      AGWindowAlgorithm aGWindowAlgorithm0 = new AGWindowAlgorithm(strand0);
      Strand strand1 = aGWindowAlgorithm0.getStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand1);
      float[] floatArray0 = new float[8];
      floatArray0[0] = (float) 5;
      floatArray0[1] = (float) 1;
      floatArray0[2] = (float) 2;
      floatArray0[3] = (float) 5;
      floatArray0[4] = (float) 2;
      floatArray0[5] = (float) 2;
      floatArray0[6] = (float) 2;
      floatArray0[7] = (float) 8715;
      gCDeviationAlgorithm0.getValues(2, 8715, floatArray0);
  }

  @Test(timeout = 4000)
  public void test04()  throws Throwable  {
      char[] charArray0 = new char[4];
      charArray0[1] = 'Z';
      charArray0[3] = 'm';
      Short short0 = new Short((short)0);
      Integer integer0 = GVTAttributedCharacterIterator.TextAttribute.ARABIC_TERMINAL;
      PartialSequence partialSequence0 = new PartialSequence(charArray0, (-1017), 128, short0, integer0);
      Bases bases0 = new Bases(partialSequence0);
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases0).when(strand0).getBases();
      doReturn(true).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.max_min_disabled = false;
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.setUserMaxMin(false);
      gCDeviationAlgorithm0.getMinimumInternal();
      gCDeviationAlgorithm0.getValueCount();
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getAverage();
      gCDeviationAlgorithm0.getDefaultStepSize(1);
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.getDefaultMinWindowSize();
      gCDeviationAlgorithm0.setUserMaxMin(false);
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.getDefaultStepSize((-2069));
      float[] floatArray0 = new float[5];
      floatArray0[0] = (float) 1;
      floatArray0[1] = (float) 'Z';
      floatArray0[2] = (-905.2115F);
      floatArray0[3] = (float) '%';
      floatArray0[4] = (float) (-2069);
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getValues(1, 1, floatArray0);
        fail("Expecting exception: NegativeArraySizeException");
      
      } catch(NegativeArraySizeException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.io.PartialSequence", e);
      }
  }

  @Test(timeout = 4000)
  public void test05()  throws Throwable  {
      FastaStreamSequence fastaStreamSequence0 = new FastaStreamSequence("internal error - unexpected exception: ");
      Bases bases0 = new Bases(fastaStreamSequence0);
      Strand strand0 = bases0.getReverseStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getMaximumInternal();
      DebugGraphics debugGraphics0 = new DebugGraphics();
      LineAttributes[] lineAttributesArray0 = new LineAttributes[0];
      gCDeviationAlgorithm0.getAverage();
  }

  @Test(timeout = 4000)
  public void test06()  throws Throwable  {
      RawStreamSequence rawStreamSequence0 = new RawStreamSequence("");
      Bases bases0 = new Bases(rawStreamSequence0);
      SimpleEntryGroup simpleEntryGroup0 = new SimpleEntryGroup(bases0);
      Bases bases1 = simpleEntryGroup0.getBases();
      Strand strand0 = bases1.getForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getAverage();
      gCDeviationAlgorithm0.getMinimumInternal();
      gCDeviationAlgorithm0.getValueCount();
      float[] floatArray0 = new float[9];
      floatArray0[0] = 1894.4F;
      floatArray0[1] = (float) 5;
      floatArray0[2] = (float) 2;
      floatArray0[3] = 2840.7537F;
      floatArray0[4] = (float) (-102);
      floatArray0[5] = (float) 1;
      floatArray0[6] = (float) 5;
      floatArray0[7] = (float) (-102);
      floatArray0[8] = 0.0F;
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getValues((-102), (-868), floatArray0);
        fail("Expecting exception: Error");
      
      } catch(Error e) {
         //
         // internal error - unexpected exception: org.evosuite.runtime.mock.java.lang.MockThrowable: start: -102 > end: -868
         //
         verifyException("uk.ac.sanger.artemis.plot.GCDeviationAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test07()  throws Throwable  {
      SimpleEntryGroup simpleEntryGroup0 = new SimpleEntryGroup();
      FeatureVector featureVector0 = new FeatureVector();
      FilteredEntryGroup filteredEntryGroup0 = new FilteredEntryGroup(simpleEntryGroup0, featureVector0, "2N)?PCxH\"aRWc");
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getMaximum();
      LineAttributes[] lineAttributesArray0 = new LineAttributes[6];
      gCDeviationAlgorithm0.setAlgorithmName("].u\"N");
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getDefaultStepSize(1);
      gCDeviationAlgorithm0.getMaximum();
      gCDeviationAlgorithm0.getMaximum();
      LineAttributes lineAttributes0 = mock(LineAttributes.class, new ViolatedAssumptionAnswer());
      lineAttributesArray0[0] = lineAttributes0;
      gCDeviationAlgorithm0.getMinimum();
      gCDeviationAlgorithm0.getDefaultStepSize(663);
      gCDeviationAlgorithm0.setUserMin((-2529));
      lineAttributesArray0[1] = lineAttributes0;
      lineAttributesArray0[2] = lineAttributes0;
      gCDeviationAlgorithm0.getMinimum();
      gCDeviationAlgorithm0.setUserMin(663);
      lineAttributesArray0[3] = lineAttributes0;
      gCDeviationAlgorithm0.getDefaultMinWindowSize();
      lineAttributesArray0[4] = lineAttributes0;
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
  }

  @Test(timeout = 4000)
  public void test08()  throws Throwable  {
      SimpleEntryGroup simpleEntryGroup0 = new SimpleEntryGroup();
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.max_min_disabled = true;
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getMinimumInternal();
      gCDeviationAlgorithm0.setUserMax(1);
      gCDeviationAlgorithm0.setUserMaxMin(true);
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.setUserMin(2);
      gCDeviationAlgorithm0.getDefaultStepSize(877);
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.getDefaultMinWindowSize();
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.getDefaultMinWindowSize();
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getMaximumInternal();
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getAverage();
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test09()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.max_min_disabled = false;
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.setRevCompDisplay(false);
      gCDeviationAlgorithm0.setRevCompDisplay(true);
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.getMinimumInternal();
      gCDeviationAlgorithm0.getDefaultStepSize(10);
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.getMaximumInternal();
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getAverage();
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test10()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.max_min_disabled = false;
      gCDeviationAlgorithm0.setAlgorithmName("   ");
      gCDeviationAlgorithm0.setScalingFlag(false);
      gCDeviationAlgorithm0.setUserMax(2076.4F);
      gCDeviationAlgorithm0.getValueCount();
      gCDeviationAlgorithm0.getDefaultMinWindowSize();
      gCDeviationAlgorithm0.getDefaultMinWindowSize();
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.getValueCount();
      gCDeviationAlgorithm0.getDefaultStepSize(2959);
  }

  @Test(timeout = 4000)
  public void test11()  throws Throwable  {
      GCDeviationAlgorithm gCDeviationAlgorithm0 = null;
      try {
        gCDeviationAlgorithm0 = new GCDeviationAlgorithm((Strand) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test12()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getDefaultMinWindowSize();
      gCDeviationAlgorithm0.getMinimumInternal();
      gCDeviationAlgorithm0.setAlgorithmName("RC5:BHZ'My&bsqp");
      gCDeviationAlgorithm0.getDefaultStepSize(2443);
      gCDeviationAlgorithm0.getMinimumInternal();
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getAverage();
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test13()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getMinimumInternal();
  }

  @Test(timeout = 4000)
  public void test14()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getMaximumInternal();
  }

  @Test(timeout = 4000)
  public void test15()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getValueCount();
      float[] floatArray0 = new float[4];
      floatArray0[0] = (float) 1;
      floatArray0[1] = (float) 2;
      floatArray0[2] = (float) (-3849);
      floatArray0[3] = (float) (-3849);
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getValues((-3849), (-3849), floatArray0);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test16()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getMinimum();
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getAverage();
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test17()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getDefaultStepSize(5000);
      gCDeviationAlgorithm0.getDefaultMinWindowSize();
  }

  @Test(timeout = 4000)
  public void test18()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getDefaultStepSize(2);
  }

  @Test(timeout = 4000)
  public void test19()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getDefaultWindowSize();
  }

  @Test(timeout = 4000)
  public void test20()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      float[] floatArray0 = new float[0];
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getValues(1, 2, floatArray0);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test21()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getAverage();
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test22()  throws Throwable  {
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.setAlgorithmName("");
      gCDeviationAlgorithm0.getDefaultMinWindowSize();
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
  }

  @Test(timeout = 4000)
  public void test23()  throws Throwable  {
      FileSystemHandling.shouldThrowIOException((EvoSuiteFile) null);
      FastaStreamSequence fastaStreamSequence0 = new FastaStreamSequence("");
      Bases bases0 = new Bases(fastaStreamSequence0);
      Strand strand0 = bases0.getReverseStrand();
      AGWindowAlgorithm aGWindowAlgorithm0 = new AGWindowAlgorithm(strand0);
      Strand strand1 = aGWindowAlgorithm0.getStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand1);
      float[] floatArray0 = new float[8];
      floatArray0[1] = (float) 1;
      floatArray0[3] = (float) 5;
      floatArray0[4] = (float) 2;
      floatArray0[5] = (float) 2;
      floatArray0[6] = (float) 2;
      floatArray0[7] = (float) 8715;
      aGWindowAlgorithm0.setScalingFlag(false);
      gCDeviationAlgorithm0.getValues(2, 8715, floatArray0);
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.getAverage();
      gCDeviationAlgorithm0.getMaximumInternal();
      float[] floatArray1 = new float[4];
      floatArray1[0] = (float) 912;
      floatArray1[1] = (float) 1;
      floatArray1[2] = (float) 5;
      floatArray1[3] = (float) 0;
      // Undeclared exception!
      gCDeviationAlgorithm0.getValues(912, 1745, floatArray1);
  }

  @Test(timeout = 4000)
  public void test24()  throws Throwable  {
      FileSystemHandling.shouldAllThrowIOExceptions();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = null;
      try {
        gCDeviationAlgorithm0 = new GCDeviationAlgorithm((Strand) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test25()  throws Throwable  {
      RawStreamSequence rawStreamSequence0 = new RawStreamSequence("internal error - unexpected exception: ");
      Bases bases0 = new Bases(rawStreamSequence0);
      Strand strand0 = bases0.getReverseStrand();
      FileSystemHandling.shouldAllThrowIOExceptions();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getMaximum();
      gCDeviationAlgorithm0.getMaximumInternal();
      Graphics graphics0 = mock(Graphics.class, new ViolatedAssumptionAnswer());
      LineAttributes[] lineAttributesArray0 = new LineAttributes[0];
      gCDeviationAlgorithm0.getValueCount();
  }

  @Test(timeout = 4000)
  public void test26()  throws Throwable  {
      GenbankStreamSequence genbankStreamSequence0 = new GenbankStreamSequence("internal error - unexpected exception: ");
      Bases bases0 = new Bases(genbankStreamSequence0);
      Strand strand0 = bases0.getForwardStrand();
      AGWindowAlgorithm aGWindowAlgorithm0 = new AGWindowAlgorithm(strand0);
      Strand strand1 = aGWindowAlgorithm0.getStrand();
      Bases bases1 = strand1.getBases();
      Strand strand2 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases1).when(strand2).getBases();
      doReturn(false).when(strand2).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand2);
      gCDeviationAlgorithm0.setUserMax(2);
      float[] floatArray0 = new float[3];
      floatArray0[0] = (float) 1;
      floatArray0[1] = (float) 2;
      floatArray0[2] = (float) 1;
      gCDeviationAlgorithm0.getValues(1, 1, floatArray0);
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getDefaultStepSize(1);
  }

  @Test(timeout = 4000)
  public void test27()  throws Throwable  {
      FastaStreamSequence fastaStreamSequence0 = new FastaStreamSequence("");
      Bases bases0 = new Bases(fastaStreamSequence0);
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn(bases0).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.max_min_disabled = true;
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getMinimumInternal();
      gCDeviationAlgorithm0.setUserMax(1);
      gCDeviationAlgorithm0.setUserMaxMin(true);
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.setUserMin(2);
      gCDeviationAlgorithm0.getDefaultStepSize(877);
      fastaStreamSequence0.getCharSubSequence(1, 5);
      gCDeviationAlgorithm0.getDefaultWindowSize();
      gCDeviationAlgorithm0.getDefaultMinWindowSize();
      fastaStreamSequence0.setFromChar(bases0.letter_index);
      gCDeviationAlgorithm0.getMaximumInternal();
      bases0.getCCount();
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getDefaultWindowSize();
      Integer integer0 = GVTAttributedCharacterIterator.TextAttribute.UNDERLINE_ON;
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.getAverage();
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.getValueCount();
  }

  @Test(timeout = 4000)
  public void test28()  throws Throwable  {
      float[] floatArray0 = new float[5];
      floatArray0[0] = (float) (-731);
      floatArray0[1] = (float) 0;
      floatArray0[2] = (float) (-731);
      floatArray0[3] = (float) 0;
      floatArray0[4] = (float) (-731);
      Marker marker0 = null;
      DatabaseStreamFeature databaseStreamFeature0 = null;
      try {
        databaseStreamFeature0 = new DatabaseStreamFeature((Feature) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.io.GFFStreamFeature", e);
      }
  }

  @Test(timeout = 4000)
  public void test29()  throws Throwable  {
      char[] charArray0 = new char[4];
      charArray0[1] = 'Z';
      charArray0[3] = 'm';
      Short short0 = new Short((short)0);
      Integer integer0 = GVTAttributedCharacterIterator.TextAttribute.ARABIC_TERMINAL;
      PartialSequence partialSequence0 = new PartialSequence(charArray0, (-1017), 128, short0, integer0);
      Bases bases0 = new Bases(partialSequence0);
      Strand strand0 = mock(Strand.class, new ViolatedAssumptionAnswer());
      doReturn((Bases) null).when(strand0).getBases();
      doReturn(false).when(strand0).isForwardStrand();
      GCDeviationAlgorithm gCDeviationAlgorithm0 = new GCDeviationAlgorithm(strand0);
      gCDeviationAlgorithm0.max_min_disabled = false;
      gCDeviationAlgorithm0.getMaximumInternal();
      gCDeviationAlgorithm0.setUserMaxMin(false);
      gCDeviationAlgorithm0.getMinimumInternal();
      gCDeviationAlgorithm0.getValueCount();
      gCDeviationAlgorithm0.getDefaultMaxWindowSize();
      // Undeclared exception!
      try { 
        gCDeviationAlgorithm0.getAverage();
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.plot.BaseAlgorithm", e);
      }
  }

  @Test(timeout = 4000)
  public void test30()  throws Throwable  {
      JPasswordField jPasswordField0 = new JPasswordField();
      DatabaseDocument databaseDocument0 = new DatabaseDocument("", jPasswordField0, "", (String) null);
      // Undeclared exception!
      try { 
        databaseDocument0.getLinePushBackReader();
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         
      }
  }
}