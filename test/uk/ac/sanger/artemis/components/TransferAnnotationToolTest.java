/* TransferAnnotationToolTest.java
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2009  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

package uk.ac.sanger.artemis.components;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.net.URL;
import java.util.List;
import java.util.Vector;
import java.util.Enumeration;
import java.util.Hashtable;
import javax.swing.JCheckBox;

import junit.framework.Assert;

import org.junit.Before;
import org.junit.Test;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.SimpleEntryGroup;

import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.util.StringVector;

public class TransferAnnotationToolTest
{
  private QualifierPanel qPanel;
  private Vector genesToTransferTo = new Vector();
  private String productStr;
  private EntryGroup entryGroup;
  private Feature originatingFeature;
  
  /**
   * Open a flat file and create the components in the TransferAnnotaionTool
   * used to control the transfer of annotation.
   */
  @Before
  public void setup()
  {
    entryGroup = readFile();
    originatingFeature = getFeature("gene", "O", "CDS", entryGroup);

    qPanel = new QualifierPanel(originatingFeature, 
        originatingFeature.getKey().getKeyString());
    Hashtable qualifierHash = (Hashtable) qPanel.getQualifierCheckBoxes();

    Enumeration enumQ = qualifierHash.keys();
    productStr = null;
    while (enumQ.hasMoreElements())
    {
      JCheckBox cbQualifierName = (JCheckBox) enumQ.nextElement();

      if (cbQualifierName.getText().equals("product"))
      {
        cbQualifierName.setSelected(true);
        Vector cbQualifierValues = (Vector) qualifierHash.get(cbQualifierName);
        for (int i = 0; i < cbQualifierValues.size(); i++)
        {
          JCheckBox cbQualifierValue = (JCheckBox) cbQualifierValues.get(i);
          cbQualifierValue.setSelected(true);
          productStr = cbQualifierValue.getText();
        }
      }
      else
        cbQualifierName.setSelected(false);
    }

    JCheckBox cbNamesToTransferTo = new JCheckBox("N", true);
    genesToTransferTo.add(cbNamesToTransferTo);
  }
  

  /**
   * Test for appending qualifiers.
   */
  @Test
  public void testAppend()
  { 
    // append product
    TransferAnnotationTool.transferAnnotation(qPanel.getQualifierCheckBoxes(),
        genesToTransferTo, originatingFeature, entryGroup, true, false, false,
        new StringBuffer(), new StringBuffer());

    testTransferResult(entryGroup, 2, productStr);
  }
  
  /**
   * Test for appending qualifiers does not create duplicates when run
   * again.
   */
  @Test
  public void testAppend2()
  { 
    testAppend();
  }
   
  /**
   * Test for overwriting qualifiers.
   */
  @Test
  public void testOverwrite()
  { 
    // check it overwrites product
    TransferAnnotationTool.transferAnnotation(qPanel.getQualifierCheckBoxes(),
        genesToTransferTo, originatingFeature, entryGroup, true, true, false,
        new StringBuffer(), new StringBuffer());

    testTransferResult(entryGroup, 1, productStr);
  }
  
  /**
   * Check the result of the annotation transfer is correct.
   * @param entryGroup
   * @param nproducts
   * @param productStr
   */
  private void testTransferResult(final EntryGroup entryGroup,
                                  final int nproducts,
                                  final String productStr)
  {
    final Feature featureN = getFeature("gene", "N", "CDS", entryGroup);
    boolean hasTransferred = false;
    StringVector values = null;
    try
    {
      values = featureN.getQualifierByName("product").getValues();
      for (int j = 0; j < values.size(); j++)
      {
        if (((String) values.get(j)).equals(productStr))
          hasTransferred = true;
      }
    }
    catch (Exception e)
    {
      Assert.fail(e.getMessage());
    }

    assertEquals("Annotation transfer (product), number of values", values.size(), nproducts);
    assertTrue("Annotation transfer (product)", hasTransferred);
  }
  
  /**
   * Return an Artemis EntryGroup by reading in an example file
   * @return
   */
  private EntryGroup readFile()
  {
    URL url = TransferAnnotationToolTest.class.getResource("/etc/af063097.embl");

    //assertEquals("Number", 4, 4);
    final EntryInformation artemisEntryInformation =
      Options.getArtemisEntryInformation();
    final Document entryDocument =
        DocumentFactory.makeDocument(url.getFile());

    try
    {
      final uk.ac.sanger.artemis.io.Entry emblEntry =
                DocumentEntryFactory.makeDocumentEntry(artemisEntryInformation,
                		entryDocument, null);
      final Entry entry = new Entry(emblEntry);
      final EntryGroup entryGroup = new SimpleEntryGroup(entry.getBases());
      entryGroup.add(entry);
      return entryGroup;
    }
    catch(uk.ac.sanger.artemis.util.OutOfRangeException oore)
    {
      Assert.fail(oore.getMessage());
    }
    catch(uk.ac.sanger.artemis.sequence.NoSequenceException nse)
    {
      Assert.fail(nse.getMessage());
    }
    catch(Exception e)
    {
      Assert.fail(e.getMessage());
    }
    return null;
  }

  /**
   * Return a feature in an entry group by searching for the feature
   * key and a qualifier with a given value.
   * @param qualifierName
   * @param qualifierValue
   * @param key
   * @param entryGroup
   * @return
   */
  private Feature getFeature(final String qualifierName,
                             final String qualifierValue, 
                             final String key, 
                             final EntryGroup entryGroup)
  {
    FeatureVector features = entryGroup.getAllFeatures();
    FeaturePredicate newPredicate = new FeatureKeyQualifierPredicate(
        new Key(key), qualifierName, qualifierValue, false, false);

    for (int i = 0; i < features.size(); i++)
    {
      Feature f = features.elementAt(i);
      if (newPredicate.testPredicate(f))
        return f;
    }
    return null;
  }
}

