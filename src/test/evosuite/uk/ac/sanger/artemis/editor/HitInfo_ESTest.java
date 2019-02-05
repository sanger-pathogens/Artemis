/*
 * This file was automatically generated by EvoSuite
 * Thu Sep 20 12:18:31 GMT 2018
 */

package uk.ac.sanger.artemis.editor;

import org.junit.Test;
import static org.junit.Assert.*;
import static org.evosuite.runtime.EvoAssertions.*;
import java.util.Vector;
import org.evosuite.runtime.EvoRunner;
import org.evosuite.runtime.EvoRunnerParameters;
import org.junit.runner.RunWith;
import uk.ac.sanger.artemis.editor.HitInfo;

@RunWith(EvoRunner.class) @EvoRunnerParameters(mockJVMNonDeterminism = true, useVFS = true, useVNET = true, resetStaticState = true, useJEE = true) 
public class HitInfo_ESTest extends HitInfo_ESTest_scaffolding {

  @Test(timeout = 4000)
  public void test00()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "`");
      hitInfo0.setQueryPosition(2502, (-1));
  }

  @Test(timeout = 4000)
  public void test01()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", ",0=");
      hitInfo0.setUngapped(",0=");
      String string0 = hitInfo0.getUngapped();
      assertEquals(",0=", string0);
  }

  @Test(timeout = 4000)
  public void test02()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("0", "");
      hitInfo0.setUngapped("");
      String string0 = hitInfo0.getUngapped();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test03()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      hitInfo0.setSubjectRange("pt=}K2\"");
      String string0 = hitInfo0.getSubjectRange();
      assertEquals("pt=}K2\"", string0);
  }

  @Test(timeout = 4000)
  public void test04()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "UniProt");
      hitInfo0.setSubjectRange("");
      String string0 = hitInfo0.getSubjectRange();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test05()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("J!DE4=$8An", "J!DE4=$8An");
      hitInfo0.setStartPosition(415);
      int int0 = hitInfo0.getStartPosition();
      assertEquals(415, int0);
  }

  @Test(timeout = 4000)
  public void test06()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      hitInfo0.setStartPosition((-2257));
      int int0 = hitInfo0.getStartPosition();
      assertEquals((-2257), int0);
  }

  @Test(timeout = 4000)
  public void test07()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      hitInfo0.setScore("y<<`o^4GHZIQ$`GDd");
      String string0 = hitInfo0.getScore();
      assertEquals("y<<`o^4GHZIQ$`GDd", string0);
  }

  @Test(timeout = 4000)
  public void test08()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      hitInfo0.setScore("");
      String string0 = hitInfo0.getScore();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test09()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("D[$p{y", "D[$p{y");
      hitInfo0.setQueryRange("D[$p{y");
      String string0 = hitInfo0.getQueryRange();
      assertEquals("D[$p{y", string0);
  }

  @Test(timeout = 4000)
  public void test10()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("^B{c/kt.5UQ.2,`CD", "^B{c/kt.5UQ.2,`CD");
      hitInfo0.setQueryRange("");
      String string0 = hitInfo0.getQueryRange();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test11()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      hitInfo0.setOverlap("M");
      String string0 = hitInfo0.getOverlap();
      assertEquals("M", string0);
  }

  @Test(timeout = 4000)
  public void test12()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("f!U", "f!U");
      hitInfo0.setOverlap("");
      String string0 = hitInfo0.getOverlap();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test13()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", ";");
      hitInfo0.setOrganism(")");
      String string0 = hitInfo0.getOrganism();
      assertEquals(")", string0);
  }

  @Test(timeout = 4000)
  public void test14()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("0", "");
      hitInfo0.setOrganism("");
      String string0 = hitInfo0.getOrganism();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test15()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("Q.ss<P+]w", "Q.ss<P+]w");
      hitInfo0.setLength("`w_<Yy`D!JZxPX");
      String string0 = hitInfo0.getLength();
      assertEquals("`w_<Yy`D!JZxPX", string0);
  }

  @Test(timeout = 4000)
  public void test16()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      hitInfo0.setLength("");
      String string0 = hitInfo0.getLength();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test17()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "qZG~");
      hitInfo0.setIdentity("qZG~");
      String string0 = hitInfo0.getIdentity();
      assertEquals("qZG~", string0);
  }

  @Test(timeout = 4000)
  public void test18()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      hitInfo0.setIdentity("");
      String string0 = hitInfo0.getIdentity();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test19()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("_K:~XO4!d &5]", "vns[3{*<dMNe(1O");
      hitInfo0.setBLASTPInfo("_K:~XO4!d &5]");
      String string0 = hitInfo0.getID();
      assertEquals("~XO4!d", string0);
      assertNotNull(string0);
  }

  @Test(timeout = 4000)
  public void test20()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("Gq!f}t`}2c1Sv5P9>", " g%");
      hitInfo0.setBLASTPInfo(" g%");
      String string0 = hitInfo0.getID();
      assertNotNull(string0);
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test21()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      String string0 = hitInfo0.getHeader();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test22()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "`");
      hitInfo0.setGoAssociation("", "");
      String string0 = hitInfo0.getGoAssociation("");
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test23()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", ",0=");
      hitInfo0.setGeneName("tp*fu:|!,\"_~2q.a ");
      String string0 = hitInfo0.getGeneName();
      assertEquals("tp*fu:|!,\"_~2q.a ", string0);
  }

  @Test(timeout = 4000)
  public void test24()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("}N&BEQ0%7l`Q:_zmYq", "");
      hitInfo0.setGeneName("");
      String string0 = hitInfo0.getGeneName();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test25()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("^B{c/kt.5UQ.2,`CD", "^B{c/kt.5UQ.2,`CD");
      hitInfo0.setEndPosition(116);
      int int0 = hitInfo0.getEndPosition();
      assertEquals(116, int0);
  }

  @Test(timeout = 4000)
  public void test26()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      hitInfo0.setEValue(";");
      String string0 = hitInfo0.getEValue();
      assertEquals(";", string0);
  }

  @Test(timeout = 4000)
  public void test27()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "|");
      hitInfo0.setEValue("");
      String string0 = hitInfo0.getEValue();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test28()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("g-&@Sy", "g-&@Sy");
      hitInfo0.setEMBL("UNIPROT");
      String string0 = hitInfo0.getEMBL();
      assertEquals("UNIPROT", string0);
  }

  @Test(timeout = 4000)
  public void test29()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "|");
      hitInfo0.setEMBL("");
      String string0 = hitInfo0.getEMBL();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test30()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("Pj]~q", "");
      hitInfo0.setEC_number("Pj]~q");
      String string0 = hitInfo0.getEC_number();
      assertEquals("Pj]~q", string0);
  }

  @Test(timeout = 4000)
  public void test31()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("g-&@Sy", "g-&@Sy");
      hitInfo0.setEC_number("");
      String string0 = hitInfo0.getEC_number();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test32()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", ",0=");
      hitInfo0.setDescription("sgq'F\u0000i'n.");
      String string0 = hitInfo0.getDescription();
      assertEquals("sgq'F\u0000i'n.", string0);
  }

  @Test(timeout = 4000)
  public void test33()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", ";");
      hitInfo0.appendDescription("");
      String string0 = hitInfo0.getDescription();
      assertEquals("", string0);
  }

  @Test(timeout = 4000)
  public void test34()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("pD|Gz%]hv4 e0_ .Q", "pD|Gz%]hv4 e0_ .Q");
      hitInfo0.setBLASTPInfo("pD|Gz%]hv4 e0_ .Q");
      String string0 = hitInfo0.getAcc();
      assertEquals("e0_", string0);
      assertNotNull(string0);
  }

  @Test(timeout = 4000)
  public void test35()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", ";");
      // Undeclared exception!
      try { 
        hitInfo0.setGoAssociation("", (String) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("java.util.Hashtable", e);
      }
  }

  @Test(timeout = 4000)
  public void test36()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("R7N`Ac) 0{I", "R7N`Ac) 0{I");
      // Undeclared exception!
      try { 
        hitInfo0.setGO((String) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.editor.HitInfo", e);
      }
  }

  @Test(timeout = 4000)
  public void test37()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("lp%sT<B>;(0CX3be", "lp%sT<B>;(0CX3be");
      // Undeclared exception!
      try { 
        hitInfo0.setFastaHitInfo((String) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.editor.HitInfo", e);
      }
  }

  @Test(timeout = 4000)
  public void test38()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      // Undeclared exception!
      try { 
        hitInfo0.setEValue((String) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
      }
  }

  @Test(timeout = 4000)
  public void test39()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      // Undeclared exception!
      try { 
        hitInfo0.setEMBL((String) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
      }
  }

  @Test(timeout = 4000)
  public void test40()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("0", "");
      // Undeclared exception!
      try { 
        hitInfo0.setBLASTPInfo((String) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("uk.ac.sanger.artemis.editor.HitInfo", e);
      }
  }

  @Test(timeout = 4000)
  public void test41()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("D%&bPg]$L%n.", "|");
      hitInfo0.setGoAssociation(" bpSK!|qA$ri  6UwEw", "GrlT -9o;l3^.3FqHdr");
      // Undeclared exception!
      try { 
        hitInfo0.getGoAssociation((String) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
         //
         // no message in exception (getMessage() returned null)
         //
         verifyException("java.util.Hashtable", e);
      }
  }

  @Test(timeout = 4000)
  public void test42()  throws Throwable  {
      HitInfo hitInfo0 = null;
      try {
        hitInfo0 = new HitInfo("xiR}j|Ku8i<x.OM'equ", (String) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
      }
  }

  @Test(timeout = 4000)
  public void test43()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("k{q;sY| ]ga-#", "0mkFS(3");
      hitInfo0.setGoAssociation("k{q;sY| ]ga-#", "");
      String string0 = hitInfo0.getGoAssociation("");
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test44()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "|");
      hitInfo0.setGoAssociation("", "6");
      String string0 = hitInfo0.getGoAssociation("");
      assertNotNull(string0);
      assertEquals("6", string0);
  }

  @Test(timeout = 4000)
  public void test45()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      String string0 = hitInfo0.getGoAssociation("");
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test46()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("v", "v");
      hitInfo0.setGoAssociation("", "");
      hitInfo0.setGoAssociation("|", "0 ");
  }

  @Test(timeout = 4000)
  public void test47()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("D%&bPg]$L%n.", "|");
      hitInfo0.setGO("GrlT -9o;l3^.3FqHdr");
  }

  @Test(timeout = 4000)
  public void test48()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "I");
      hitInfo0.setGO("I");
      Vector vector0 = hitInfo0.getGO();
      assertFalse(vector0.isEmpty());
  }

  @Test(timeout = 4000)
  public void test49()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("/K+OBVxYD[&|Tk@7Bg;", " ");
      hitInfo0.setEValue("ei:@hw;");
  }

  @Test(timeout = 4000)
  public void test50()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      hitInfo0.setEMBL(".;");
  }

  @Test(timeout = 4000)
  public void test51()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("/K+OBVxYD[&|Tk@7Bg;", " ");
      hitInfo0.appendDescription("(EC ");
  }

  @Test(timeout = 4000)
  public void test52()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "I");
      hitInfo0.setDescription("p");
      hitInfo0.appendDescription("UNIPROT");
  }

  @Test(timeout = 4000)
  public void test53()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "I");
      hitInfo0.setDescription("");
      hitInfo0.appendDescription("UNIPROT");
  }

  @Test(timeout = 4000)
  public void test54()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("", "");
      hitInfo0.setOrganism("");
      hitInfo0.setOrganism("pt=}K2\"");
  }

  @Test(timeout = 4000)
  public void test55()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("pD|Gz%]hv4 e0_ .Q", "pD|Gz%]hv4 e0_ .Q");
      hitInfo0.setBLASTPInfo("pD|Gz%]hv4 e0_ .Q");
      hitInfo0.setFastaHitInfo("pD|Gz%]hv4 e0_ .Q");
  }

  @Test(timeout = 4000)
  public void test56()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("_K:~XO4!d &5]", "vns[3{*<dMNe(1O");
      hitInfo0.setFastaHitInfo("_K:~XO4!d &5]");
  }

  @Test(timeout = 4000)
  public void test57()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("pD|Gz%]hv4 e0_ .Q", "pD|Gz%]hv4 e0_ .Q");
      hitInfo0.setFastaHitInfo("pD|Gz%]hv4 e0_ .Q");
  }

  @Test(timeout = 4000)
  public void test58()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("n.$\"u", "n.$\"u");
      hitInfo0.setBLASTPInfo(" bW7s{,_zK$w~)Ze ");
  }

  @Test(timeout = 4000)
  public void test59()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("_K:~XO4!d &5]", "vns[3{*<dMNe(1O");
      hitInfo0.setBLASTPInfo("_K:~XO4!d &5]");
      String string0 = hitInfo0.getDB();
      assertNotNull(string0);
      assertEquals("_K", string0);
  }

  @Test(timeout = 4000)
  public void test60()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("2E", "2E");
      hitInfo0.setBLASTPInfo("<y~DDKeNj.7%j 4^;");
  }

  @Test(timeout = 4000)
  public void test61()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("pD|Gz%]hv4 e0_ .Q", "pD|Gz%]hv4 e0_ .Q");
      hitInfo0.setBLASTPInfo("pD|Gz%]hv4 e0_ .Q");
      hitInfo0.setBLASTPInfo("pD|Gz%]hv4 e0_ .Q");
  }

  @Test(timeout = 4000)
  public void test62()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("(EC ", "blastp");
      hitInfo0.setOrganism("F964[CM");
      hitInfo0.setOrganism("blastp");
  }

  @Test(timeout = 4000)
  public void test63()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("fasta", "fasta");
      // Undeclared exception!
      try { 
        hitInfo0.appendDescription((String) null);
        fail("Expecting exception: NullPointerException");
      
      } catch(NullPointerException e) {
      }
  }

  @Test(timeout = 4000)
  public void test64()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "`");
      String string0 = hitInfo0.getLength();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test65()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("jrjwrO", "");
      String string0 = hitInfo0.getEValue();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test66()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("0", "");
      String string0 = hitInfo0.getSubjectRange();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test67()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("g-&@Sy", "g-&@Sy");
      String string0 = hitInfo0.getScore();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test68()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("0", "");
      String string0 = hitInfo0.getOrganism();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test69()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("g-&@Sy", "g-&@Sy");
      String string0 = hitInfo0.getEMBL();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test70()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("/K+OBVxYD[&|Tk@7Bg;", " ");
      Vector vector0 = hitInfo0.getQueryPosition();
      assertTrue(vector0.isEmpty());
  }

  @Test(timeout = 4000)
  public void test71()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("/K+OBVxYD[&|Tk@7Bg;", " ");
      String string0 = hitInfo0.getDB();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test72()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "`");
      String string0 = hitInfo0.getIdentity();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test73()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("g-&@Sy", "g-&@Sy");
      String string0 = hitInfo0.getGeneName();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test74()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("/K+OBVxYD[&|Tk@7Bg;", " ");
      Vector vector0 = hitInfo0.getGO();
      assertNull(vector0);
  }

  @Test(timeout = 4000)
  public void test75()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("jrjwrO", "");
      String string0 = hitInfo0.getID();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test76()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "`");
      String string0 = hitInfo0.getDescription();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test77()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("/K+OBVxYD[&|Tk@7Bg;", " ");
      hitInfo0.setQueryPosition(2434, 2434);
      Vector vector0 = hitInfo0.getQueryPosition();
      assertEquals("[2434, 2434]", vector0.toString());
  }

  @Test(timeout = 4000)
  public void test78()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "`");
      int int0 = hitInfo0.getEndPosition();
      assertEquals(0, int0);
  }

  @Test(timeout = 4000)
  public void test79()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("J!DE4=$8An", "J!DE4=$8An");
      String string0 = hitInfo0.getQueryRange();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test80()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "`");
      int int0 = hitInfo0.getStartPosition();
      assertEquals(0, int0);
  }

  @Test(timeout = 4000)
  public void test81()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("R7N`Ac) 0{I", "R7N`Ac) 0{I");
      String string0 = hitInfo0.getEC_number();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test82()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("/K+OBVxYD[&|Tk@7Bg;", " ");
      String string0 = hitInfo0.getAcc();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test83()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("0", "");
      String string0 = hitInfo0.getUngapped();
      assertNull(string0);
  }

  @Test(timeout = 4000)
  public void test84()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo("/K+OBVxYD[&|Tk@7Bg;", " ");
      String string0 = hitInfo0.getHeader();
      assertNotNull(string0);
  }

  @Test(timeout = 4000)
  public void test85()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "`");
      hitInfo0.setEndPosition((-158));
      int int0 = hitInfo0.getEndPosition();
      assertEquals((-158), int0);
  }

  @Test(timeout = 4000)
  public void test86()  throws Throwable  {
      HitInfo hitInfo0 = new HitInfo(":", "`");
      String string0 = hitInfo0.getOverlap();
      assertNull(string0);
  }
}