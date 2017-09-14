/* DatabaseTreeNode.java
 *
 * created: October 2006
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
*/

package uk.ac.sanger.artemis.components.database;

import java.awt.Frame;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.UnsupportedFlavorException;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.tree.DefaultMutableTreeNode;

import org.gmod.schema.organism.Organism;
import org.gmod.schema.organism.OrganismProp;
import org.gmod.schema.sequence.Feature;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.components.Splash;
import uk.ac.sanger.artemis.util.DatabaseDocument;

import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.io.IOException;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
* File node for local file tree manager
*/
public class DatabaseTreeNode extends DefaultMutableTreeNode 
                 implements Transferable, Serializable
{
  private static final long serialVersionUID = 1L;
  public static final DataFlavor DATABASETREENODE = 
          new DataFlavor(DatabaseTreeNode.class, "Work Package");
  public static final DataFlavor STRING_DATA_FLAVOUR = 
          new DataFlavor(String.class, "text/plain");
  private static final DataFlavor flavors[] = 
         { DATABASETREENODE, DataFlavor.stringFlavor };
  
  private String featureId;
  private String featureType;  // e.g. chromosome, mitochondrial_chromosome
  private transient Organism organism;
  private String organismCommonName;
  private boolean isLeaf = false;
  private String userName;
  private static transient DatabaseDocument dbDoc;
  private boolean explored = false;
   
  protected DatabaseTreeNode(final String name)
  { 
    super(name);
  }
  
  private DatabaseTreeNode(final String name, 
                          final boolean isLeaf)
  {
    super(name);
    this.isLeaf = isLeaf;
  }
    
  /**
   * Top node constructor
   * @param name
   * @param isLeaf
   * @param organism
   * @param userName
   * @param dbDoc
   */
  protected DatabaseTreeNode(final String name, 
                          final boolean isLeaf,
                          final Organism organism,
                          final String userName,
                          final DatabaseDocument dbDoc)
  { 
    super(name);
    this.isLeaf   = isLeaf;
    this.organism = organism;
    this.userName = userName;
    this.dbDoc    = dbDoc;
    setOrganismCommonName();
  }

  /**
   * Leaf constructor
   * @param name
   * @param organism
   * @param featureId
   * @param userName
   */
  private DatabaseTreeNode(final String name,
                          final Organism organism,
                          final String featureId,
                          final String featureType,
                          final String userName)
  { 
    super(name);
    this.organism  = organism;
    this.featureId = featureId;
    this.featureType = featureType;
    this.userName  = userName;
    if(getOrganism() != null)
      setOrganismCommonName();
  }
  
  public String getOrganismCommonName()
  {
	return organismCommonName;
  }
  
  private void setOrganismCommonName()
  {
	this.organismCommonName = getOrganism().getCommonName();
	if(organismCommonName == null || organismCommonName.equals(""))
	  organismCommonName = getOrganism().getGenus() + "." + getOrganism().getSpecies();
  }
  
  /**
   * Use the OrganismProps to set the translation table and
   * determine in this is a read only entry.
   * @param op
   * @return
   */
  public static boolean setOrganismProps(Set<OrganismProp> op, final boolean isMitochondrial)
  {
    Splash splash = getSplash();
    boolean readOnly = false;
    final Iterator<OrganismProp> it = op.iterator();
    while (it.hasNext()) 
    {
      OrganismProp organismProp = it.next();
      if(splash != null)
      {
        if( (isMitochondrial && 
             organismProp.getCvTerm().getName().equals("mitochondrialTranslationTable")) ||
            (!isMitochondrial && 
             organismProp.getCvTerm().getName().equals("translationTable")))
          splash.setTranslationTable(organismProp.getValue());
      }
      
     if(organismProp.getCvTerm().getName().equals("frozen") &&
        organismProp.getValue().equals("yes"))
       readOnly = true;
    } 
    return readOnly;
  }
  
  private static Splash getSplash()
  {
    Frame[] frames = JFrame.getFrames();
    for(int i=0;i<frames.length;i++)
    {
      if(frames[i] instanceof Splash)
        return (Splash)frames[i];
    }
    return null;
  }
  
  /** @return   true if node is a directory */
  public boolean getAllowsChildren() { return !isLeaf; }
  /** @return         true if node is a file */
  public boolean isLeaf() { return isLeaf; }
  /** @return         true if node is a directory */
  public boolean isDirectory() { return !isLeaf; }
  /** @return         true if explored */
  public boolean isExplored() { return explored; }

  /**
   * Explore the tree node if this is a node with child nodes.
   */
  public void explore()
  {
    if(isLeaf)
      return;
    
    List<Feature> sequenceList = dbDoc.getResidueFeatures(new Integer(getOrganism().getOrganismId()));
    Hashtable<String, DatabaseTreeNode> sequenceNode = new Hashtable<String, DatabaseTreeNode>();
    
    for(int i=0;i<sequenceList.size(); i++)
    {
      Feature f = sequenceList.get(i);
      DatabaseTreeNode typeNode;
      if(!sequenceNode.containsKey(f.getCvTerm().getName()))
      {
        typeNode = new DatabaseTreeNode(f.getCvTerm().getName(), false);
        add(typeNode);
        sequenceNode.put(f.getCvTerm().getName(), typeNode);
      }
      else
        typeNode = (DatabaseTreeNode) sequenceNode.get(f.getCvTerm().getName());
      
      DatabaseTreeNode seqNode = new DatabaseTreeNode(
          f.getUniqueName(), getOrganism(), 
          Integer.toString(f.getFeatureId()),
          f.getCvTerm().getName(), getUserName());

      seqNode.isLeaf = true;
      typeNode.add(seqNode);
      typeNode.explored = true;
    }
    explored = true;
    
    if(System.getProperty("database_manager_cache_off") == null)
      writeCache();
  }
  
  /**
   * Write out this node to the cache directory
   */
  private void writeCache()
  {
    try
    {
      File dir = new File(Options.CACHE_PATH);
      if(!dir.exists())
        dir.mkdirs();
      FileOutputStream fos = new FileOutputStream(Options.CACHE_PATH +
          ((String)dbDoc.getLocation()).replaceAll("[/:=\\?]", "_"));
      
      ObjectOutputStream out = new ObjectOutputStream(fos);
      out.writeObject(this);
      out.close();
    }
    catch(Exception ex)
    {
      JOptionPane.showMessageDialog(null, ex.getMessage());
    }
  }
  
    
// Transferable
  public DataFlavor[] getTransferDataFlavors()
  {
    return flavors;
  }

  public boolean isDataFlavorSupported(DataFlavor f)
  {
    if(f.equals(DATABASETREENODE) || f.equals(DataFlavor.stringFlavor))
      return true;
    return false;
  }

  public Object getTransferData(DataFlavor d)
        throws UnsupportedFlavorException, IOException
  {
    if(d.equals(DATABASETREENODE))
      return new DatabaseTreeNode((String)getUserObject(), 
                                  organism,
                                  getFeatureId(), getFeatureType(), getUserName());
    else if(d.equals(DataFlavor.stringFlavor))
    {
      String name = getOrganism().getCommonName();
      if(name == null || name.equals(""))
        name = getOrganism().getGenus() + "." + getOrganism().getSpecies();
      return name+":featureId="+getFeatureId();
    }
    else throw new UnsupportedFlavorException(d);
  }

  public String getFeatureId()
  {
    return featureId;
  }
  
  public String getFeatureType()
  {
    return featureType;
  }

  protected void setDbDoc(DatabaseDocument dbDoc)
  {
    DatabaseTreeNode.dbDoc = dbDoc;
  }
  
  public Organism getOrganism()
  {
    if(organism == null && organismCommonName != null)
      organism = dbDoc.getOrganismByCommonName(organismCommonName);

    return organism;
  }

  public String getUserName()
  {
    return userName;
  }
  
//Serializable
  private void writeObject(java.io.ObjectOutputStream out) throws IOException
  {
    out.defaultWriteObject();
  }

  private void readObject(java.io.ObjectInputStream in)
    throws IOException, ClassNotFoundException
  {
    in.defaultReadObject();
  }
}
