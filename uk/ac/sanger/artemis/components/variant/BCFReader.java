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

package uk.ac.sanger.artemis.components.variant;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.List;
import java.util.Vector;
import java.util.regex.Pattern;

import net.sf.samtools.util.BlockCompressedInputStream;


class BCFReader
{
  public static final int TAD_LIDX_SHIFT = 13; // linear index shift
  private static Pattern formatPattern = Pattern.compile("[^0-9]+");
  private BlockCompressedInputStream is;
  private FileInputStream indexFileStream;
  private List<BCFIndex> idx;
  
  // header information and names
  private List<String> seqNames;
  private List<String> sampleNames;
  private int nsamples;
  private String metaData;
  
  /**
   * @param bcf  BCF file
   * @throws IOException
   */
  public BCFReader(File bcf) throws IOException
  {
    is = new BlockCompressedInputStream(bcf);
    readHeader();
    
    File bcfIndex = new File(bcf.getAbsolutePath()+".bci");
    indexFileStream = new FileInputStream(bcfIndex);
    idx = loadIndex();
  }

  protected void seek(long off) throws IOException
  {
    is.seek(off);
  }
  
  protected VCFRecord next(int bid, int beg, int end) throws IOException
  {
    try
    {
      VCFRecord bcfRecord = readVCFRecord();
      if(bcfRecord.pos >= beg && bcfRecord.pos <= end)
        return bcfRecord;
      else if(bcfRecord.pos < beg)
      {
        while( (bcfRecord = readVCFRecord()).pos <= beg )
        {
          if(bcfRecord.pos >= beg && bcfRecord.pos <= end)
            return bcfRecord;
        }
      }
    } 
    catch(Exception e)
    {
      if(is.read() != -1)  // eof
        e.printStackTrace();
    }
    
    return null;
  }
  
  protected void close() throws IOException
  {
    is.close();
    indexFileStream.close();
  }
  
  private void readHeader() throws IOException
  {
    byte[] magic = new byte[4];
    is.read(magic);

    String line = new String(magic);
    if(!line.equals("BCF\4"))
    {
      throw new IOException("Not BCF format.");
    }

    // sequence names
    seqNames = getList(readInt(is));

    // sample names
    sampleNames = getList(readInt(is));
    nsamples = sampleNames.size();

    int len = readInt(is);   
    byte meta[] = new byte[len];
    is.read(meta);

    StringBuffer buff = new StringBuffer();
    for(int i=0; i<meta.length; i++)
      buff.append((char)meta[i]);

    metaData = buff.toString();
  }
  
  protected String headerToString()
  {
    StringBuffer head = new StringBuffer();
    head.append("##fileformat=VCFv4.0\n");
    head.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
    
    for(int i=0; i<sampleNames.size(); i++)
      head.append(sampleNames.get(i)+" ");
    return head.toString();
  }
  
  /**
   * Given the length of the concatenated names, that are NULL padded
   * construct a list of the Strings.
   * @param len
   * @return
   * @throws IOException
   */
  private List<String> getList(int len) throws IOException
  {
    byte names[] = new byte[len];
    is.read(names);
    
    List<String> list = new Vector<String>();
    StringBuffer buff = new StringBuffer();
    for(int i=0; i<names.length; i++)
    {
      if(names[i] != 0)
        buff.append((char)names[i]);
      else if(buff.length() > 0)
      {
        list.add(buff.toString());
        buff = new StringBuffer();
      }
    }
    return list;
  }
  
  private VCFRecord readVCFRecord() throws IOException
  {
    VCFRecord bcfRecord = new VCFRecord();
    bcfRecord.seqID = seqNames.get(readInt(is));    
    bcfRecord.pos = readInt(is)+1;
    bcfRecord.quality = readFloat(is);
    
    int slen = readInt(is);
    byte[] str = new byte[slen];
    is.read(str);

    getParts(str, bcfRecord);
    
    if(formatPattern.matcher(bcfRecord.format).matches())
    {     
      int nc  = 3;
      if(bcfRecord.alt.equals("."))
        nc = 1;

      String fmts[] = bcfRecord.format.split(":");
      for(int j=0; j<fmts.length; j++)
      {
        int nb = getByteSize(fmts[j],nc);
        str = new byte[nb];
        is.read(str);
        
        final String value;
        if(fmts[j].equals("GT"))
          value = getGTString(str[0]);
        else if(fmts[j].equals("PL"))
          value = getPLString(str, nc);
        else if(fmts[j].equals("DP")||fmts[j].equals("SP")||fmts[j].equals("GQ"))
          value = Integer.toString(byteToInt(str[0]));
        else
          value = "";
        bcfRecord.data.put(fmts[j], value);
      }
      
    }
    
    return bcfRecord;
  }
  
  /**
   * Make a string from the byte array (expanding NULL padding) to
   * determine the parts:- ID+REF+ALT+FILTER+INFO+FORMAT.
   * @param b
   * @return
   */
  private void getParts(byte[] b, VCFRecord bcfRecord)
  {
    StringBuffer buff = new StringBuffer();

    for(int i=0; i<b.length; i++)
    {
      if(b[i] == 0)
      {
        if(i == 0)
        {
          buff.append(". ");
          continue;
        }
        
        if(b[i-1] == 0 || (i < b.length-1 && b[i+1] == 0))
        {
          i++;
          buff.append(" . ");
        }
        else
          buff.append(" ");
        continue;
      }
      buff.append((char)b[i]);
    }
    
    String parts[] = buff.toString().replace("  ", " ").split(" ");

    bcfRecord.ID   = parts[0];
    bcfRecord.ref  = parts[1];
    bcfRecord.alt  = parts[2];
    bcfRecord.filter = parts[3];
    bcfRecord.info   = parts[4];
    bcfRecord.format = parts[5];
  }
  
  /**
   * DP uint16 t[n] Read depth
   * GL float[n*x] Log10 likelihood of data; x = m(m+1)/2 , m = #{alleles}
   * GT uint8 t[n] phase<<6 | allele1<<3 | allele2
   * GQ uint8 t[n] Genotype quality
   * HQ uint8 t[n*2] Haplotype quality
   * PL uint8 t[n*x] Phred-scaled likelihood of data
   * misc int32 t+char* NULL padded concatenated strings (integer equal to the length)
   * @param tag
   * @param nsamples
   * @param nc
   * @return
   */
  private int getByteSize(String tag, int nc)
  {
    if(tag.equals("DP"))        // Read depth
      return 2*nsamples;        // uint16_t[n]
    else if(tag.equals("GL"))   // Log10 likelihood
      return 4*nsamples*nc;     // float[nsamples*x]
    else if(tag.equals("GT"))   // phase<<6 | allele1<<3 | allele2
      return nsamples;          // uint8_t[n]
    else if(tag.equals("GQ"))   // Genotype quality
      return nsamples;          // uint8_t[n]
    else if(tag.equals("HQ"))   // Haplotype quality
      return 2*nsamples;        // uint8_t[n*2]
    else if(tag.equals("PL"))   // Phred-scaled likelihood
      return nsamples*nc;       // uint8_t[n*x]
    else if(tag.equals("SP"))   // 
      return nsamples;          // uint8_t[n]
    else                        // misc
      return 4*nsamples;        // uint32_t+char*
  }
  
  
  private String getPLString(byte[] b, int nc)
  {
    StringBuffer buff = new StringBuffer();
    for(int i=0;i<b.length; i++)
    {
      buff.append(byteToInt(b[i]));
      if(i<b.length-1)
        buff.append(",");
    }
    return buff.toString();
  }
  
  /**
   * GT genotype, allele values separated by Ó/Ó or Ò|Ó, i.e.
   * unphased or phased.
   * @param b
   * @return
   */
  private String getGTString(byte b)
  {
    return ((b >> 3) + ( (b >> 6 == 1) ? "|" : "/") + byteToInt(b));
  }
 
  private int byteToInt(byte b)
  {
    return (int)(b & 0xFF);
  }
  
  protected List<BCFIndex> loadIndex() throws IOException
  {
    BlockCompressedInputStream is = new BlockCompressedInputStream(indexFileStream);
    byte[] magic = new byte[4];
    is.read(magic);
    
    if(!new String(magic).equals("BCI\4"))
      System.err.println("Not a BCF index file:: "+new String(magic));
    
    int n = readInt(is);
    List<BCFIndex> idx = new Vector<BCFIndex>(n);
    
    for(int i=0; i<n; i++)
    {
      BCFIndex idx2 = new BCFIndex();
      idx2.n = readInt(is);
      idx2.index2_offset = new long[idx2.n];
      
      for(int j=0; j<idx2.n; j++)
        idx2.index2_offset[j] = readLong(is);

      idx.add(idx2);
    }
    return idx;
  }
  
  protected long queryIndex(int tid, int beg)
  {
    long min_off = -1;
    if (beg < 0) 
      beg = 0;
   
    long offset[] = idx.get(tid).index2_offset;
    int i;

    for(i = beg>>TAD_LIDX_SHIFT; i < idx.get(tid).n && offset[i] == 0; ++i);
    min_off = (i == idx.get(tid).n)? offset[idx.get(tid).n-1] : offset[i];
       
    return min_off;
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

  protected String getMetaData()
  {
    return metaData;
  }

  public static void main(String args[])
  {
    try
    {
      int sbeg = 0;
      int send = Integer.MAX_VALUE;
      if(args.length > 1)
      {
        sbeg = Integer.parseInt(args[1]);
        send = Integer.parseInt(args[2]);
      }
      
      BCFReader reader = new BCFReader(new File(args[0]));
      int bid = 0;
      
      long off = reader.queryIndex(bid, sbeg);
      reader.seek(off);

      System.out.println(reader.headerToString());
      VCFRecord bcfRecord;
      while( (bcfRecord = reader.next(bid, sbeg, send)) != null )
        System.out.println(bcfRecord.toString());
      
      reader.close();
    }
    catch (IOException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
}


class BCFIndex
{
  int n;
  long index2_offset[];
}