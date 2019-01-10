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
import java.util.Map;
import java.util.Vector;

import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.variant.BCFReader.BCFReaderIterator;

public abstract class AbstractVCFReader
{
  private boolean vcf_v4 = false;
  
  protected String[] sampleNames;
  protected abstract String[] getSeqNames();
  protected abstract String getFileName();
  protected int nsamples = -1;
  
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
      record = VCFRecord.parse(s, getNumberOfSamples());
      
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
   * @param manualHash
   * @param vcfFileName
   * @param writer
   * @param vcfView
   * @param features
   * @throws IOException
   */
  protected static void write(final Map<String, Boolean> manualHash, 
                              final String vcfFileName,
                              final int vcfIndex,
                              Writer writer, 
                              VCFview vcfView, 
                              FeatureVector features) throws IOException
  {
    // FILTER LINES
    if(IOUtils.isBCF(vcfFileName))
    {
      BCFReader reader = new BCFReader(vcfFileName);
      // FIX for old style BAM files
      AbstractVCFReader readers[] = vcfView.getVcfReaders();
      for(int i=0; i<readers.length; i++)
      {
        if(readers[i].getFileName().equals(vcfFileName))
          reader.newBCF = ((BCFReader)readers[i]).newBCF;
      }
      
      // add header
      String hdr = replaceFilterLines(reader.headerToString(), FilteredPanel.getHeader());
      writer.write( hdr );
      
      int sbeg = 0;
      int send = Integer.MAX_VALUE;
      VCFRecord record;
      
      while( (record = reader.nextRecord(null, sbeg, send)) != null)
      {
        int basePosition = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
        VCFFilter.setFilterString(manualHash, record, vcfView, basePosition, features, reader, vcfIndex);
        writer.write(record.toString()+"\n");
      }
      writer.close();
      reader.close();
      return;
    }
    
    final TabixReader tr = new TabixReader(vcfFileName);
    String line;
    boolean headerEnd = true;
    final StringBuffer buffHeader = new StringBuffer();

    while ((line = tr.readLine()) != null)
    {
      if(line.startsWith("##"))
      {
        if(!line.startsWith("##FILTER"))
          writer.write(line+'\n');
        buffHeader.append(line+'\n');
        continue;
      }
      else if(headerEnd)
      {
        writer.write(FilteredPanel.getHeader());
        headerEnd = false;
      }
      
      if(line.startsWith("#"))
      {
        writer.write(line+'\n');
        buffHeader.append(line+'\n');
        tr.setHeader(buffHeader.toString());
        continue;
      }
      
      VCFRecord record = VCFRecord.parse(line, tr.getNumberOfSamples());
      int basePosition = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
      VCFFilter.setFilterString(manualHash, record, vcfView, basePosition, features, tr, vcfIndex);
      writer.write(record.toString()+'\n');
    }
    writer.close();
  }
  
  /**
   * Return the header with the new Filter lines.
   * @param hdrLines
   * @param filterLines
   * @return
   */
  private static String replaceFilterLines(final String hdrLines, final String filterLines)
  {
    final StringBuffer buff = new StringBuffer();
    final BufferedReader readerStr = new BufferedReader(new StringReader(hdrLines));
    String str;
    try
    {
      while ((str = readerStr.readLine()) != null)
      {
        if (!str.startsWith("##Filter="))
        {
          if(str.startsWith("#CHROM"))
          {
            buff.append(filterLines);
            buff.append(str+"\n");
          }
          else
            buff.append(str+"\n");
        }
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }

    return buff.toString();
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
  
  protected int getNumberOfSamples()
  {
    if(nsamples < 1)
    {
      nsamples = 1;
      BufferedReader reader = new BufferedReader(new StringReader(header));
      try
      {
        String ln;
        while ((ln = reader.readLine()) != null)
        {
          if (ln.startsWith("#CHROM"))
          {
            int index = ln.indexOf("FORMAT");
            if(index > -1)
            {
              sampleNames = ln.substring(index+7).trim().split("[ \\t]");
              nsamples = sampleNames.length;
            }
          }
        }
      }
      catch (IOException e)
      {
        System.err.println("Problem calculating the number of samples.");
      }
    }
    return nsamples;
  }
  
  protected String[] getSampleNames()
  {
    return sampleNames;
  }

  protected List<HeaderLine> getFORMAT()
  {
    return getListOfLines("FORMAT");
  }
  
  protected List<HeaderLine> getFILTER()
  {
    return getListOfLines("FILTER");
  }
  
  protected List<HeaderLine> getINFO()
  {
    return getListOfLines("INFO");
  }
  
  /**
   * Return a list of the lines in the header. Each line is represented as a hash of
   * the key value pairs.
   * @param lineType
   * @return
   */
  private List<HeaderLine> getListOfLines(String lineType)
  {
    List<HeaderLine> listOfType = new Vector<HeaderLine>();

    try
    {
      BufferedReader reader = new BufferedReader(new StringReader(getHeader()));
      String str;
      while ((str = reader.readLine()) != null)
      {
        String origLine = new String(str);
        if (str.startsWith("##"+lineType))
        {
          listOfType.add(new HeaderLine(origLine, lineType, getLineHash(lineType, str)));
        }
      }
    }
    catch (IOException ioe)
    {
      ioe.printStackTrace();
    }
    return listOfType;
  }
  
  protected static Hashtable<String, String> getLineHash(String lineType, String str)
  {
    Hashtable<String, String> hash = new Hashtable<String, String>();
    hash.put("lineType", lineType);
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
    return hash;
  }
}