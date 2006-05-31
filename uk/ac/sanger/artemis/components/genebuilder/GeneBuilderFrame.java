/* GeneBuilderFrame.java
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2006  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneBuilderFrame.java,v 1.1 2006-05-31 09:49:07 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import javax.swing.*;
import java.awt.*;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.components.QualifierTextArea;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.StreamQualifier;
import uk.ac.sanger.artemis.util.StringVector;

public class GeneBuilderFrame extends JFrame
{
  
  public GeneBuilderFrame(final Feature gene_feature)
  {
    super("Artemis Gene Builder: " + gene_feature.getIDString() +
          (gene_feature.isReadOnly() ?
          "  -  (read only)" :
          ""));
    
    QualifierTextArea qualifier_text_area = new QualifierTextArea();
    GFFStreamFeature gff_feature = (GFFStreamFeature)gene_feature.getEmblFeature();
    GeneComponentTree tree = new GeneComponentTree(gff_feature.getChadoGene(),
                                                   qualifier_text_area);
    getContentPane().add(new JScrollPane(tree), BorderLayout.WEST);
    
    GeneViewerPanel viewer = new GeneViewerPanel(gff_feature.getChadoGene());
    getContentPane().add(viewer, BorderLayout.CENTER);

    qualifier_text_area.setWrapStyleWord(true);
    qualifier_text_area.setText(getQualifierString(gene_feature));
    
    getContentPane().add(new JScrollPane(qualifier_text_area), BorderLayout.SOUTH);
    
    pack();
    setVisible(true);
  }

  /**
   *  Return a string containing one qualifier per line.  These are the
   *  original qualifiers, not the qualifiers from the qualifier_text_area.
   **/
  protected static String getQualifierString(final Feature feature) 
  {
    final StringBuffer buffer = new StringBuffer();
    final QualifierVector qualifiers = feature.getQualifiers();

    for(int qualifier_index = 0; qualifier_index < qualifiers.size();
        ++qualifier_index) 
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(qualifier_index);

      final QualifierInfo qualifier_info =
                       getEntryInformation(feature).getQualifierInfo(this_qualifier.getName());

      final StringVector qualifier_strings =
                       StreamQualifier.toStringVector(qualifier_info, this_qualifier);

      for(int value_index = 0; value_index < qualifier_strings.size();
          ++value_index)
      {
        final String qualifier_string = (String)qualifier_strings.elementAt(value_index);
        buffer.append(qualifier_string + "\n");
      }
    }

    return buffer.toString();
  }
  
  /**
   *  Return the EntryInformation object of the entry containing the feature.
   **/
  protected static EntryInformation getEntryInformation(final Feature feature) 
  {
    return feature.getEntry().getEntryInformation();
  }

}