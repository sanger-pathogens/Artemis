/* ArtemisUtil
 *
 * created: 2007
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2007  Genome Research Limited
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

package uk.ac.sanger.artemis.chado;

import java.util.List;
import java.util.Vector;

import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;
import org.gmod.schema.sequence.FeatureCvTermProp;

import uk.ac.sanger.artemis.components.Splash;

public class ArtemisUtils
{
  
  /**
   * Delete featureCvTerm and update associated feature_cvterm.rank's 
   * if appropriate
   * @param dao
   * @param featureCvTerm
   */
  public static void deleteFeatureCvTerm(GmodDAO dao, FeatureCvTerm featureCvTerm)
  {
    List featureCvTerms = dao.getFeatureCvTermsByFeature(featureCvTerm.getFeature());
       
    List featureCvTermDbXRefs = new Vector();
    
    if(featureCvTerm.getFeatureCvTermDbXRefs() != null &&
       featureCvTerm.getFeatureCvTermDbXRefs().size() > 0)
      featureCvTermDbXRefs = (List)featureCvTerm.getFeatureCvTermDbXRefs();
    
    List featureCvTermProps = new Vector();
    
    if(featureCvTerm.getFeatureCvTermProps()!= null &&
        featureCvTerm.getFeatureCvTermProps().size() > 0)
      featureCvTermProps = (List)featureCvTerm.getFeatureCvTermProps();
    
    // delete feature_cvterm and update ranks if appropriate
    FeatureCvTerm deleteme = null;
    Vector rankable = null;
    
    for(int i=0; i<featureCvTerms.size(); i++)
    {
      FeatureCvTerm this_feature_cvterm = (FeatureCvTerm)featureCvTerms.get(i);
      
      if(this_feature_cvterm.getCvTerm().getName().equals( 
         featureCvTerm.getCvTerm().getName() )  &&
         this_feature_cvterm.getCvTerm().getCv().getName().equals( 
             featureCvTerm.getCvTerm().getCv().getName() ))
      {     
         List this_featureCvTermDbXRefs = (List)this_feature_cvterm.getFeatureCvTermDbXRefs();
         List this_featureCvTermProps   = (List)this_feature_cvterm.getFeatureCvTermProps();
         
         if(this_featureCvTermDbXRefs == null)
           this_featureCvTermDbXRefs = new Vector();
         
         
         if(this_featureCvTermDbXRefs.size() != featureCvTermDbXRefs.size() ||
            featureCvTermProps.size() != this_featureCvTermProps.size())
         {
           if(rankable == null)
             rankable = new Vector();
           
           rankable.add(this_feature_cvterm);
           continue;
         }
         
         boolean found = true;
         for(int j=0; j<this_featureCvTermDbXRefs.size(); j++)
         {
           FeatureCvTermDbXRef fcd = (FeatureCvTermDbXRef)this_featureCvTermDbXRefs.get(j);
           FeatureCvTermProp fcp   = (FeatureCvTermProp)this_featureCvTermProps.get(j);
           if(!containsFeatureCvTermDbXRef(fcd, featureCvTermDbXRefs) ||
              !containsFeatureCvTermProp(fcp, featureCvTermProps))
           {
             Splash.logger4j.debug(fcp.getCvTerm().getName()+" "+fcp.getValue());
             
             found = false;
             break;
           }
         }
         if(!found)
         {
           if(rankable == null)
             rankable = new Vector();
           
           rankable.add(this_feature_cvterm);
           continue;
         }
         
         deleteme = this_feature_cvterm;
      }
    }
    dao.delete(deleteme);
    
    if(rankable != null)
    {
      // feature_cvterm.rank may need updating for those stored here
      for(int i=0; i<rankable.size(); i++)
      {
        FeatureCvTerm fc = (FeatureCvTerm)rankable.get(i);
        
        if(fc.getRank() == i)
          continue;
        
        Splash.logger4j.debug("UPDATE rank for "+ fc.getCvTerm().getCv().getName() + "   rank = " +
                              fc.getRank()+" -> "+i);
        fc.setRank(i);
        dao.merge(fc);
      }
    }
  }
  
  /**
   * Return true if the list contains a given feature_cvterm_dbxref
   * @param fcd
   * @param featureCvTermDbXRefs
   * @return
   */
  private static boolean containsFeatureCvTermDbXRef(FeatureCvTermDbXRef fcd, 
                                              List featureCvTermDbXRefs)
  {
    for(int i=0; i<featureCvTermDbXRefs.size(); i++)
    {
      FeatureCvTermDbXRef this_fcd = (FeatureCvTermDbXRef)featureCvTermDbXRefs.get(i);
      if( this_fcd.getDbXRef().getAccession().equals( fcd.getDbXRef().getAccession() ) &&
          this_fcd.getDbXRef().getDb().getName().equals( fcd.getDbXRef().getDb().getName() ))
        return true;
    }
    return false;
  }
  
  /**
   * Return true if the list contains a given feature_cvterm_prop
   * @param fcp
   * @param featureCvTermProps
   * @return
   */
  private static boolean containsFeatureCvTermProp(FeatureCvTermProp fcp, 
                                            List featureCvTermProps)
  {
    for(int i = 0; i < featureCvTermProps.size(); i++)
    {
      FeatureCvTermProp this_fcp = (FeatureCvTermProp)featureCvTermProps.get(i);
      if(this_fcp.getValue().equals(fcp.getValue()))
        return true;
    }
    return false;
  }
}