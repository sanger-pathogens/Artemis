/* 
 * This file is part of Artemis
 *
 * Copyright(C) 2013  Genome Research Limited
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
package uk.ac.sanger.artemis.io;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Vector;

import javax.swing.JOptionPane;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;

class GffToEMBL
{
  /**
   * Convert GFF to EMBL
   * @param inGff
   * @param outDir
   * @param emblSubmission
   * @param flatten
   * @param gzip
   */
  public GffToEMBL(final String inGff, 
                   final String outDir, 
                   final boolean emblSubmission,
                   final boolean flatten,
                   final boolean gzip)
  {
    final Entry entry = getEntry(inGff);
    if(entry == null || !(entry instanceof GFFDocumentEntry) )
    {
      JOptionPane.showMessageDialog(null, 
          "No GFF entry found.", "Error", JOptionPane.ERROR_MESSAGE);
      return;
    }
    
    EntryGroup egrp = null;
    try
    {
      uk.ac.sanger.artemis.Entry artEntry = new  uk.ac.sanger.artemis.Entry(entry);
      egrp = new SimpleEntryGroup(artEntry.getBases());
      egrp.add(artEntry);
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    }
    catch (NoSequenceException e) {} // no sequence


    if(egrp == null) // no sequence found
    {
      final FeatureVector features = entry.getAllFeatures();
      final HashSet<String> seqNames = new HashSet<String>();
      for(Feature f: features)
        seqNames.add(((GFFStreamFeature)f).getGffSeqName());
      
      for(String seq: seqNames)
      {
        Entry newEntry = getEntryForSeqName(seq, features);
        writeEMBL(newEntry, seq, outDir, emblSubmission, flatten, gzip);
      }
    }
    else
    {
      writeEMBL(entry, entry.getName(), outDir, emblSubmission, flatten, gzip);
    }
  }

  /**
   * Write EMBL file
   * @param entry
   * @param seqName
   * @param outDir
   * @param emblSubmission
   * @param gzip
   */
  private void writeEMBL(
      final Entry entry,
      final String seqName,
      final String outDir,
      final boolean emblSubmission,
      final boolean flatten,
      final boolean gzip)
  {
    final EntryInformation entryInfo;
    if(emblSubmission)
      entryInfo =Options.getDBEntryInformation();
    else
      entryInfo = Options.getArtemisEntryInformation();
    
    if(!flatten)
    {
      final FeatureVector features = entry.getAllFeatures();
      for(int i=0; i<features.size(); i++)
        ReadAndWriteEntry.addAllKeysQualifiers(entryInfo, features.elementAt(i));

      if(entry instanceof GFFDocumentEntry)
        ReadAndWriteEntry.addQualifierToEntryInfo(entryInfo, 
          (String)PublicDBDocumentEntry.getDatabaseQualifiersToRemove()[0]);
    }
      
    try
    {
      final EmblDocumentEntry emblEntry =
          new EmblDocumentEntry (entryInfo, entry, true);
      FileDocument out = new FileDocument(new File(outDir+File.separator+
          seqName+".embl"+(gzip ? ".gz" : "")));
      emblEntry.save(out);
      
      System.out.println("Written... "+out.getFile().getAbsolutePath());
    }
    catch (EntryInformationException e)
    {
      e.printStackTrace();
    }
    catch (IOException e)
    {
      JOptionPane.showMessageDialog(null, 
          e.getMessage(), "I/O Error", JOptionPane.ERROR_MESSAGE);
      e.printStackTrace();
    }
  }
  
  /**
   * Get a new entry for a sequence
   * @param seqName
   * @param features
   * @return
   */
  private Entry getEntryForSeqName(
                final String seqName, 
                final FeatureVector features)
  {
    final GFFDocumentEntry newEntry = new GFFDocumentEntry(null);
    for(Feature f: features)
    {
      final GFFStreamFeature gff = (GFFStreamFeature)f;
      if(gff.getGffSeqName().equals(seqName))
      try
      {
        newEntry.forcedAdd(new GFFStreamFeature(gff));
      }
      catch (ReadOnlyException e)
      {
        e.printStackTrace();
      }
    }
    return newEntry;
  }
  
  private Entry getEntry(final String fileName)
  {
    try
    {
      // move away index file
      File f = new File(fileName);
      if(IndexedGFFDocumentEntry.isIndexed( f ))
      {
        f = new File(f.getAbsolutePath() + ".tbi");
        File tmp = new File(f.getAbsolutePath() + ".old");
        f.renameTo(tmp);
        f = tmp;
      }
      else
        f = null;
        
      final Document doc = DocumentFactory.makeDocument(fileName);
      final Entry entry =  DocumentEntryFactory.makeDocumentEntry(
          Options.getArtemisEntryInformation(),doc,null);
      
      // move back index file
      if(f != null)
      {
        String name = f.getAbsolutePath();
        f.renameTo(new File(name.substring(0, name.length()-4)));
      }
      return entry;
    }
    catch(EntryInformationException e) 
    {
      JOptionPane.showMessageDialog(null, 
          e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
    }
    catch(IOException e) 
    {
      JOptionPane.showMessageDialog(null, 
          e.getMessage(), "I/O Error", JOptionPane.ERROR_MESSAGE);
      e.printStackTrace();
    }
    return null;
  }
  
  public static void main(String args[])
  {
    if( (args != null && args.length == 1 && args[0].startsWith("-h")) ||
        (args == null || args.length < 1))
    {
      System.out.println("-h\tshow help");
      System.out.println("-s\tspace separated list of GFF files to read and write out");
      System.out.println("-o\toutput directory");
      System.out.println("-f\t[y|n] flatten the gene model, default is y");
      System.out.println("-z\t[y|n] gzip output, default is y");
      System.out.println("-a\t[y|n] for EMBL submission format change to n, default is y");
      System.exit(0);
    }

    boolean gzip = true;
    boolean emblSubmission = true;
    boolean flatten = true;
    Vector<String> files = new java.util.Vector<String>();;
    String outDir = System.getProperty("user.dir"); // working directory
    for(int i=0; i<args.length; i++)
    {
      String s = args[i];
      if(s.equals("-z"))
      {
        if(i + 1 < args.length && args[i + 1].toLowerCase().equals("n"))
          gzip = false;
      }
      else if(s.equals("-a"))
      {
        if(i + 1 < args.length && args[i + 1].toLowerCase().equals("y"))
          emblSubmission = false;
      }
      else if(s.equals("-f"))
      {
        if(i + 1 < args.length && args[i + 1].toLowerCase().equals("n"))
          flatten = false;
      }
      else if(s.equals("-o"))
      {
        if(i + 1 < args.length)
          outDir= args[i + 1];
      }
      else if(args[i].toLowerCase().equals("-s"))
      {
        for(int j = i + 1; j < args.length; j++)
        {
          if(args[j].startsWith("-"))
            break;
          files.add(args[j]);
        }
      }
    }
    
    final File outDirFile = new File(outDir);
    if(!outDirFile.exists() && !outDirFile.mkdir())
    {
      JOptionPane.showMessageDialog(null, "Problems writing to "+
               outDirFile.getAbsolutePath(), "Error", JOptionPane.ERROR_MESSAGE);
      System.exit(0);
    }

    for(String fn: files)
      new GffToEMBL(fn, outDir, emblSubmission, flatten, gzip);
    System.exit(0);
  }
}
