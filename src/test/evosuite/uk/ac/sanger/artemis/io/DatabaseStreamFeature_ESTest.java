/*
 * This file was automatically generated by EvoSuite
 * Thu Sep 20 14:09:53 GMT 2018
 */

package uk.ac.sanger.artemis.io;

import org.junit.Test;
import static org.junit.Assert.*;
import static org.evosuite.shaded.org.mockito.Mockito.*;
import static org.evosuite.runtime.EvoAssertions.*;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.evosuite.runtime.ViolatedAssumptionAnswer;
import org.junit.runner.RunWith;
import uk.ac.sanger.artemis.io.DatabaseStreamFeature;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.QualifierVector;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, useJEE = true) 
public class DatabaseStreamFeature_ESTest extends DatabaseStreamFeature_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test0()  throws Throwable  {
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
  public void test1()  throws Throwable  {
      Key key0 = mock(Key.class, new ViolatedAssumptionAnswer());
      doReturn((String) null).when(key0).getKeyString();
      Location location0 = mock(Location.class, new ViolatedAssumptionAnswer());
      doReturn((String) null).when(location0).toString();
      DatabaseStreamFeature databaseStreamFeature0 = new DatabaseStreamFeature(key0, location0, (QualifierVector) null);
      assertNotNull(databaseStreamFeature0);
      assertTrue(databaseStreamFeature0.isVisible());
      assertFalse(databaseStreamFeature0.isLazyLoaded());
      assertNull(databaseStreamFeature0.getGffSource());
      assertFalse(databaseStreamFeature0.isReadOnly());
      assertNull(databaseStreamFeature0.getGffSeqName());
      
      DatabaseStreamFeature databaseStreamFeature1 = new DatabaseStreamFeature(databaseStreamFeature0);
      assertNotNull(databaseStreamFeature1);
      assertFalse(databaseStreamFeature1.equals((Object)databaseStreamFeature0));
      assertTrue(databaseStreamFeature0.isVisible());
      assertFalse(databaseStreamFeature0.isLazyLoaded());
      assertNull(databaseStreamFeature0.getGffSource());
      assertFalse(databaseStreamFeature0.isReadOnly());
      assertNull(databaseStreamFeature0.getGffSeqName());
      assertNull(databaseStreamFeature1.getGffSeqName());
      assertFalse(databaseStreamFeature1.isLazyLoaded());
      assertNull(databaseStreamFeature1.getGffSource());
      assertFalse(databaseStreamFeature1.isReadOnly());
      assertTrue(databaseStreamFeature1.isVisible());
  }

  @Test(timeout = 4000)
  public void test2()  throws Throwable  {
      Key key0 = mock(Key.class, new ViolatedAssumptionAnswer());
      doReturn((String) null).when(key0).getKeyString();
      Location location0 = mock(Location.class, new ViolatedAssumptionAnswer());
      doReturn((String) null).when(location0).toString();
      DatabaseStreamFeature databaseStreamFeature0 = new DatabaseStreamFeature(key0, location0, (QualifierVector) null);
      assertNotNull(databaseStreamFeature0);
      assertTrue(databaseStreamFeature0.isVisible());
      assertNull(databaseStreamFeature0.getGffSource());
      assertFalse(databaseStreamFeature0.isLazyLoaded());
      assertNull(databaseStreamFeature0.getGffSeqName());
      assertFalse(databaseStreamFeature0.isReadOnly());
      
      DatabaseStreamFeature databaseStreamFeature1 = (DatabaseStreamFeature)databaseStreamFeature0.copy();
      assertNotNull(databaseStreamFeature1);
      assertNotSame(databaseStreamFeature0, databaseStreamFeature1);
      assertNotSame(databaseStreamFeature1, databaseStreamFeature0);
      assertFalse(databaseStreamFeature1.equals((Object)databaseStreamFeature0));
      assertTrue(databaseStreamFeature0.isVisible());
      assertNull(databaseStreamFeature0.getGffSource());
      assertFalse(databaseStreamFeature0.isLazyLoaded());
      assertNull(databaseStreamFeature0.getGffSeqName());
      assertFalse(databaseStreamFeature0.isReadOnly());
      assertNull(databaseStreamFeature1.getGffSeqName());
      assertFalse(databaseStreamFeature1.isLazyLoaded());
      assertFalse(databaseStreamFeature1.isReadOnly());
      assertNull(databaseStreamFeature1.getGffSource());
      assertTrue(databaseStreamFeature1.isVisible());
  }

  @Test(timeout = 4000)
  public void test3()  throws Throwable  {
      Key key0 = mock(Key.class, new ViolatedAssumptionAnswer());
      doReturn((String) null).when(key0).getKeyString();
      DatabaseStreamFeature databaseStreamFeature0 = null;
      try {
        databaseStreamFeature0 = new DatabaseStreamFeature(key0, (Location) null, (QualifierVector) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.io.GFFStreamFeature", e);
      }
  }
}