/* 
 * created: 2010
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2010  Genome Research Limited
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
import java.io.Writer;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.JOptionPane;

import net.sf.picard.reference.FastaSequenceIndex;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFileFactory;

import uk.ac.sanger.artemis.io.Entry;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.util.CacheHashMap;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.URLDocument;

public class IndexFastaStream extends StreamSequence 
{
  private IndexedFastaSequenceFile indexSeqFile;
  private FastaSequenceIndex fastaIndex;
  private int len;
  private String contig;
  private CacheHashMap basesCache;  // used by charAt()
  
  public IndexFastaStream(Entry entry)
  {
    DocumentEntry doc = (DocumentEntry)entry;
    if(doc instanceof URLDocument)
    {
      //URL url = (URL)((URLDocument)doc).getLocation();
      // not supported yet
    }
    else
    {
      File fasta = ((FileDocument)doc.getDocument()).getFile();
      File parentDir = fasta.getParentFile();
      File fastaIndexFile;
      if(parentDir != null)
        fastaIndexFile = new File(parentDir.getAbsolutePath(), fasta.getName() + ".fai");
      else
        fastaIndexFile = new File(fasta.getName() + ".fai");
      
      fastaIndex = new FastaSequenceIndex(fastaIndexFile);
      
      try
      {
        indexSeqFile = new IndexedFastaSequenceFile(fasta, fastaIndex);
      }
      catch(IllegalArgumentException ie)
      {
        JOptionPane.showConfirmDialog(null, 
            "Expecting fasta extensions:\n"+
            ReferenceSequenceFileFactory.FASTA_EXTENSIONS.toString()+
            "\n"+ie.getMessage(), 
            "Error", JOptionPane.ERROR_MESSAGE);
      }
    }
    
    setContigByIndex(0);
  }
  
  /**
   * Test if the entry contains an indexed sequence
   * @param entry
   * @return
   */
  protected static boolean isIndexed(Entry entry)
  {
    try
    {
      if (entry instanceof DocumentEntry
          && ((DocumentEntry) entry).getDocument() instanceof FileDocument)
      {
        File fasta = ((FileDocument) ((DocumentEntry) entry).getDocument()).getFile();
        File parentDir = fasta.getParentFile();
        File fastaIndexFile;
        if(parentDir != null)
          fastaIndexFile = new File(parentDir.getAbsolutePath(), fasta.getName() + ".fai");
        else
          fastaIndexFile = new File(fasta.getName() + ".fai");
        if (fastaIndexFile.exists())
        {
          try
          {
            setExtensions();
          }
          catch(UnsupportedClassVersionError e)
          {
            System.err.println("Java version "+System.getProperty("java.version")+
                " does not support indexed fasta - use Java 1.6 or higher.");
            return false;
          }
          return true;
        }
      }
    }
    catch(Exception e)
    {
      e.printStackTrace();
    }
    return false;
  }
  
  protected boolean useIndex()
  {
    if(fastaIndex.size() == 1)
      return true;

    Object[] possibleValues = { "Use index", "Concatenate Sequences" };
    int sel = JOptionPane.showOptionDialog(null, 
        "Use the FASTA index file (to increase the performance)\n"+ 
        "or concatenate the sequences together?",
        "Indexed FASTA Found", 
        JOptionPane.DEFAULT_OPTION, JOptionPane.WARNING_MESSAGE,
        null, possibleValues, possibleValues[0]);
   
    return sel == 0;
  }
  
  /**
   * Add to supported FASTA allowed suffixes.
   */
  private static void setExtensions()
  {
    if(ReferenceSequenceFileFactory.FASTA_EXTENSIONS.contains(".dna"))
      return;
    ReferenceSequenceFileFactory.FASTA_EXTENSIONS.add(".dna");
    ReferenceSequenceFileFactory.FASTA_EXTENSIONS.add(".seq");
    ReferenceSequenceFileFactory.FASTA_EXTENSIONS.add(".fas");
    ReferenceSequenceFileFactory.FASTA_EXTENSIONS.add(".ffn");
  }
  
  public void setContigByIndex(int seqIndex) 
  {
    /*ReferenceSequence ref = getReferenceSequence(seqIndex);
    len = ref.length();
    contig = ref.getName();*/
    
    len = getLengthByIndex(seqIndex);
    contig = getContigByIndex(seqIndex);
    basesCache = null;
  }

  /**
   *  Return a the given range of bases as a String.  Returns an empty
   *  sequence if the end position is less than the start position.
   *  @param start The start base of the range.
   *  @param end The end base of the range.
   **/
  public String getSubSequence(int start, int end) 
  {
    byte b[] = indexSeqFile.getSubsequenceAt(contig, start, end).getBases();
    return new String(b).toLowerCase();
  }
  
  public char[] getCharSubSequence(int start, int end) 
  {
    return getSubSequence(start, end).toCharArray();
  }
  
  /**
   * Used by AddMenu.markAmbiguities() to retrieve the sequence character
   * at a specified position
   */
  public char charAt(final int i)
  {
    if(basesCache == null)
      basesCache = new CacheHashMap(250,50);
      
    if(!basesCache.containsKey(i))
    {
      int end = i+249;
      if(end > length())
        end = length();
      String seq = getSubSequence(i, end);
      for(int idx=0; idx<seq.length(); idx++)
        basesCache.put(i+idx, seq.charAt(idx));
    }

    return (Character)basesCache.get(i);
  }
  
  private int getLengthByIndex(int seqIndex)
  {
    Iterator it = fastaIndex.iterator();
    int i = 0;
    while(it.hasNext())
    {
      Object obj = it.next();
      if(i == seqIndex)
      {
        String parts[] = obj.toString().split(";");
        for(String part: parts)
        {
          if(part.trim().startsWith("size"))
            return Integer.parseInt(part.substring(5).trim());
        }
      }
      i++;
    }
    return -1;
  }
  
  private String getContigByIndex(int seqIndex)
  {
    Iterator it = fastaIndex.iterator();
    int i = 0;
    while(it.hasNext())
    {
      Object obj = it.next();
      if(i == seqIndex)
      {
        String c = obj.toString().split(" ")[1];
        if(c.endsWith(";"))
          c = c.substring(0, c.length()-1);
        return c;
      }
      i++;
    }
    return null;
  }
  
  public Vector<String> getContigs()
  {
    Iterator it = fastaIndex.iterator();
    Vector<String> contigs = new Vector<String>();
    while(it.hasNext())
    {
      String contig = it.next().toString().split(";")[0];
      if(contig.startsWith("contig "))
        contig = contig.substring(6).trim();
      contigs.add(  contig );
    }
    return contigs;
  }
  
  public ReferenceSequence getReferenceSequence(int seqIndex)
  {
    int i = 0;
    ReferenceSequence ref;
    indexSeqFile.reset();
    while( (ref=indexSeqFile.nextSequence()) != null ) 
    {  
      if(i == seqIndex)
        return ref;
       i++;
    }
    return null;
  }
  

  /**
   *  Returns the length of the sequence in bases.
   **/
  public int length() 
  {
    return len;
  }
  
  @Override
  public StreamSequence copy()
  {
    // TODO Auto-generated method stub
    return null;
  }

  @Override
  public int getFormatType()
  {
    return StreamSequenceFactory.INDEXED_FASTA_FORMAT;
  }
  
  public void setFromChar(final char dna[])
  {
    JOptionPane.showMessageDialog(null,"Read only sequence.", 
        "Warning", JOptionPane.WARNING_MESSAGE);
    throw new RuntimeException(new ReadOnlyException());
  }

  @Override
  public void writeToStream(Writer writer) throws IOException
  {
    // TODO Auto-generated method stub
  }
  
  public IndexedFastaSequenceFile getIndexSeqFile()
  {
    return indexSeqFile;
  }
  
  public FastaSequenceIndex getFastaIndex()
  {
    return fastaIndex;
  }
  
  public String getContig()
  {
    return contig;
  }
  
  public static void main(String args[])
  {
    EntryInformation new_entry_information =
      new SimpleEntryInformation(Options.getArtemisEntryInformation());

    try
    {
      Entry emblEntry = null;
      if(args[0].startsWith("http:"))
      {
        
        
      }
      else
      {
        final uk.ac.sanger.artemis.Entry entry =  
          new uk.ac.sanger.artemis.Entry(EntryFileDialog.getEntryFromFile(
                    null, new FileDocument(new File(args[0])),
                    new_entry_information, true));
        emblEntry = entry.getEMBLEntry();
      }
      
      IndexFastaStream istream = new IndexFastaStream(emblEntry);
      System.out.println(istream.getCharSubSequence(1, 8000));
    }
    catch (Exception e)
    {
      e.printStackTrace();
    }
  }
}