package uk.ac.sanger.artemis.editor;

import static org.junit.Assert.assertEquals;

import java.awt.GraphicsEnvironment;

import javax.swing.JDesktopPane;

import org.junit.Test;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.io.Utils;

/** 
 * 
 * Unit tests for the EvidenceViewer class
 * 
 * @author kp11
 *
 */
public class EvidenceViewerTest
{

	/**
	 * Test update of Pfam URL. 
	 */
	@Test
	public void testSetPfamUrl()
	{
		
		// ignore if in headless mode with no x11
	    if(GraphicsEnvironment.isHeadless())
	      return;
	    
	    // Given
	    
	    final EntryGroup egrp = Utils.getEntryGroup("/data/Pf3D7_01_02_v3.gff.gz");

	    FeatureVector features = egrp.getAllFeatures();
	    
	    Feature feature = Utils.getCDSFeatureByIdPrefix("PF3D7_0103500.1", features);
	    
	    FeatureVector featureVec = new FeatureVector();
	    
	    Feature overlapFeature = Utils.getCDSFeatureByIdPrefix("PF3D7_0103500.1", features);
	    
	    featureVec.add(overlapFeature);
	    
	    // When
	    
	    EvidenceViewer viewer = new EvidenceViewer(feature, featureVec, new JDesktopPane());
	    
	    // Then
	    
	    String pfamDb = Options.getOptions().getDatabaseHyperlinkProperty(Options.PFAM_HYPERLINK_PROPERTY_NAME);
		
	    assertEquals("https://pfam.xfam.org/family/", pfamDb); // Also exercise the new Options method.
	    assertEquals(pfamDb, viewer.pfamUrl);
	}
}
