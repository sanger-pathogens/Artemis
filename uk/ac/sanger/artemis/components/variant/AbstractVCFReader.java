/* AbstractReader
 *
 * created: July 2011
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

package uk.ac.sanger.artemis.components.variant;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringReader;
import java.io.Writer;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.variant.BCFReader.BCFReaderIterator;

public abstract class AbstractVCFReader
{
  private boolean vcf_v4 = false;
  protected abstract String[] getSeqNames();
  protected abstract String getFileName();
  
  private BCFReaderIterator bcfIterator = null;
  private TabixReader.Iterator tabixIterator = null;
  private String header;
  
  /**
   * Read and return the next record.
   * @param chr     sequence name
   * @param sbeg    start base
   * @param send    end base
   * @return
   * @throws IOException
   */
  public VCFRecord getNextRecord(String chr, int sbeg, int send) throws IOException
  {
    VCFRecord record;
    if(this instanceof BCFReader)
    {
      if(bcfIterator == null)
        bcfIterator = ((BCFReader)this).query(chr, sbeg, send);

      record = bcfIterator.next();
      if(record == null)
        bcfIterator = null;
    }
    else
    {
      if(tabixIterator == null)
      {
        try
        {
          tabixIterator = ((TabixReader)this).query(chr+":"+sbeg+"-"+send);
        }
        catch(ArrayIndexOutOfBoundsException aob)
        {
          System.err.println(chr+":"+sbeg+"-"+send+" not found in "+((TabixReader)this).getFileName());
        }
      }
      if(tabixIterator == null)
        return null;

      String s = tabixIterator.next();
      if(s == null)
      {
        tabixIterator = null;
        return null;
      }
      record = VCFRecord.parse(s);
      
      if(record == null)
        tabixIterator = null;
    }
    return record;
  }
  
  
  protected static int readInt(final InputStream is) throws IOException {
    byte[] buf = new byte[4];
    is.read(buf);
    return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getInt();
  }
  
  
  protected static float readFloat(final InputStream is) throws IOException {
    byte[] buf = new byte[4];
    is.read(buf);
    return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getFloat();
  }

  protected static long readLong(final InputStream is) throws IOException {
    byte[] buf = new byte[8];
    is.read(buf);
    return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getLong();
  }
  
  protected String getName()
  {
    if(getFileName() == null)
      return null;
    File f = new File(getFileName());
    return f.getName();
  }
  
  /**
   * Export VCF file
   * @param writer
   * @param vcfView
   * @param features
   * @throws IOException
   */
  protected void write(Writer writer, VCFview vcfView, FeatureVector features) throws IOException
  {
    writer.write( getHeader()+"\n" );
    if(this instanceof BCFReader)
    {
      BCFReader reader = new BCFReader(getFileName());
      int sbeg = 0;
      int send = Integer.MAX_VALUE;
      VCFRecord record;
      
      while( (record = reader.nextRecord(null, sbeg, send)) != null)
      {
        int basePosition = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
        if( !vcfView.showVariant(record, features, basePosition, false) )
          continue;

        writer.write(record.toString()+"\n");
      }
      writer.close();
      reader.close();
      return;
    }
    
    TabixReader tr = new TabixReader(getFileName());
    String line;
    while ((line = tr.readLine()) != null)
    {
      if(line.startsWith("#"))
        continue;
      
      VCFRecord record = VCFRecord.parse(line);
      int basePosition = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
      if( !vcfView.showVariant(record, features, basePosition, tr.isVcf_v4()) )
        continue;
      writer.write(line+'\n');
    }
    writer.close();
  }
  
  /**
   * @return the vcf_v4
   */
  protected boolean isVcf_v4()
  {
    return vcf_v4;
  }
  
  /**
   * @param vcfV4 the vcf_v4 to set
   */
  protected void setVcf_v4(boolean vcfV4)
  {
    vcf_v4 = vcfV4;
  }
  
  protected void setHeader(String header)
  {
    this.header = header;
  }
  
  protected String getHeader()
  {
    return header;
  }
  
  protected List<Hashtable<String, String>> getINFO()
  {
    return getListOfLines("INFO");
  }
  
  private List<Hashtable<String, String>> getListOfLines(String lineType)
  {
    List<Hashtable<String, String>> listOfType = 
        new Vector<Hashtable<String, String>>();

    try
    {
      BufferedReader reader = new BufferedReader(new StringReader(getHeader()));
      String str;
      while ((str = reader.readLine()) != null)
      {
        if (str.startsWith("##"+lineType))
        {
          System.out.println(str);
          Hashtable<String, String> hash = new Hashtable<String, String>();
          str = str.substring(lineType.length() + 4, str.length()-1);

          String parts[] = str.split(","); 
          for(int i=0; i<parts.length; i++)
          {
            if(!parts[i].startsWith("Description"))
            {
              String thisPart[] = parts[i].split("=");
              if(thisPart.length == 2)
              {
                hash.put(thisPart[0], thisPart[1]);
              }
            }
            else if(parts[i].startsWith("Description"))
            {
              String thisPart[] = parts[i].split("=");
              
              if(thisPart.length == 2)
              {
                //search for closing quote
                while(thisPart[0].equals("Description") &&
                      thisPart[1].startsWith("\"") &&
                      thisPart[1].indexOf("\"",2) == -1 &&
                      i+1 < parts.length)
                {
                  thisPart[1] = thisPart[1] + "," + parts[++i];
                }

                if(thisPart[1].startsWith("\""))
                  thisPart[1] = thisPart[1].substring(1);
                if(thisPart[1].endsWith("\""))
                  thisPart[1] = thisPart[1].substring(0,thisPart[1].length()-1);
                hash.put(thisPart[0], thisPart[1]);
              }
            }
          }
          listOfType.add(hash);
        }
      }
    }
    catch (IOException ioe)
    {
      ioe.printStackTrace();
    }
    return listOfType;
  }
}