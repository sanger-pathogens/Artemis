/*
 * This file was automatically generated by EvoSuite
 * Thu Sep 20 14:08:23 GMT 2018
 */

package org.gmod.schema.pub;

import org.junit.Test;
import static org.junit.Assert.*;
import static org.evosuite.shaded.org.mockito.Mockito.*;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.evosuite.runtime.ViolatedAssumptionAnswer;
import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.pub.Pub;
import org.gmod.schema.pub.PubProp;
import org.junit.runner.RunWith;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, useJEE = true) 
public class PubProp_ESTest extends PubProp_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test0()  throws Throwable  {
      Integer integer0 = new Integer(1);
      PubProp pubProp0 = new PubProp((CvTerm) null, (Pub) null, "g_Wm", integer0);
      Integer integer1 = pubProp0.getRank();
      assertEquals(1, (int)integer1);
  }

  @Test(timeout = 4000)
  public void test1()  throws Throwable  {
      CvTerm cvTerm0 = mock(CvTerm.class, new ViolatedAssumptionAnswer());
      Pub pub0 = mock(Pub.class, new ViolatedAssumptionAnswer());
      Integer integer0 = new Integer(0);
      PubProp pubProp0 = new PubProp(cvTerm0, pub0, "", integer0);
      Integer integer1 = pubProp0.getRank();
      assertEquals(0, (int)integer1);
  }

  @Test(timeout = 4000)
  public void test2()  throws Throwable  {
      Integer integer0 = new Integer((-1225));
      PubProp pubProp0 = new PubProp((CvTerm) null, (Pub) null, "g\"m@k]AHXYNuQ", integer0);
      Integer integer1 = pubProp0.getRank();
      assertEquals((-1225), (int)integer1);
  }

  @Test(timeout = 4000)
  public void test3()  throws Throwable  {
      CvTerm cvTerm0 = mock(CvTerm.class, new ViolatedAssumptionAnswer());
      Pub pub0 = mock(Pub.class, new ViolatedAssumptionAnswer());
      Integer integer0 = new Integer((-1096));
      PubProp pubProp0 = new PubProp(cvTerm0, pub0, "", integer0);
      Integer integer1 = pubProp0.getRank();
      assertEquals((-1096), (int)integer1);
  }

  @Test(timeout = 4000)
  public void test4()  throws Throwable  {
      Integer integer0 = new Integer((-1225));
      PubProp pubProp0 = new PubProp((CvTerm) null, (Pub) null, "g\"m@k]AHXYNuQ", integer0);
      CvTerm cvTerm0 = pubProp0.getCvTerm();
      assertNull(cvTerm0);
  }

  @Test(timeout = 4000)
  public void test5()  throws Throwable  {
      PubProp pubProp0 = new PubProp();
      Integer integer0 = pubProp0.getRank();
      assertNull(integer0);
  }

  @Test(timeout = 4000)
  public void test6()  throws Throwable  {
      CvTerm cvTerm0 = mock(CvTerm.class, new ViolatedAssumptionAnswer());
      doReturn("").when(cvTerm0).toString();
      Pub pub0 = mock(Pub.class, new ViolatedAssumptionAnswer());
      PubProp pubProp0 = new PubProp(cvTerm0, pub0, "");
      CvTerm cvTerm1 = pubProp0.getCvTerm();
      assertSame(cvTerm1, cvTerm0);
  }
}
