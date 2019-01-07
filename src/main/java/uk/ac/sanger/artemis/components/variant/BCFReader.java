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
import java.net.URL;
import java.util.List;
import java.util.Vector;
import java.util.regex.Pattern;

import uk.ac.sanger.artemis.util.FTPSeekableStream;

import htsjdk.samtools.util.BlockCompressedInputStream;

class BCFReader extends AbstractVCFReader
{
  public static final int TAD_LIDX_SHIFT = 13; // linear index shift
  private static Pattern formatPattern = Pattern.compile("[^0-9]+");
  private BlockCompressedInputStream is;
  private InputStream indexFileStream;
  private List<BCFIndex> idx;
  
  // header information and names
  private String[] seqNames;

  private String metaData;
  private String fileName;
  
  //
  // Nasty work around for backward compatibility with old BCF files.
  // Assume new BCF file unless ArrayIndexOutOfBoundsException thrown. 
  // The newer BCF files have:
  // SP type of long 
  protected boolean newBCF = true; 

  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(BCFReader.class);
  
  /**
   * @param bcf  BCF file or url
   * @throws IOException
   */
  public BCFReader(String bcf) throws IOException
  {
    if(bcf.startsWith("http"))
    {
      URL bcfURL = new URL(bcf);
      is = new BlockCompressedInputStream(bcfURL);
      indexFileStream = new URL(bcf+".bci").openStream();
      fileName = bcfURL.getFile();
    }
    else if(bcf.startsWith("ftp"))
    {
      URL bcfURL = new URL(bcf);
      FTPSeekableStream fss = new FTPSeekableStream(bcfURL);
      is = new BlockCompressedInputStream(fss);
      indexFileStream = new URL(bcf+".bci").openStream();
      fileName = bcfURL.getFile();
    }
    else 
    {
      File bcfFile = new File(bcf);
      is = new BlockCompressedInputStream(bcfFile);
      indexFileStream = new FileInputStream(new File(bcf+".bci"));
      fileName = bcfFile.getAbsolutePath();
    }
    idx = loadIndex();
    readHeader();
    logger4j.debug(bcf);
  }

  protected void seek(long off) throws IOException
  {
    is.seek(off);
  }
  
  protected VCFRecord nextRecord(String chr, int beg, int end) throws IOException
  {
    try
    {
      VCFRecord bcfRecord = readVCFRecord();
      if(chr != null && !bcfRecord.getChrom().equals(chr))
        return null;

      if(bcfRecord.getPos() >= beg && bcfRecord.getPos() <= end)
        return bcfRecord;
      else if(bcfRecord.getPos() < beg)
      {
        while( (bcfRecord = readVCFRecord()).getPos() <= beg )
        {
          if(chr != null && !bcfRecord.getChrom().equals(chr))
            return null;

          if(bcfRecord.getPos() >= beg && bcfRecord.getPos() <= end)
            return bcfRecord;
        }
        if(bcfRecord.getPos() >= beg && bcfRecord.getPos() <= end)
          return bcfRecord;
      }
    }
    catch(ArrayIndexOutOfBoundsException ae)
    {
      if(newBCF)
      {
        newBCF = false;
        logger4j.debug("This looks like an old style BCF: "+fileName);
      }
      else
        throw ae;
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
    seqNames = getArray(readInt(is));

    // sample names
    sampleNames = getArray(readInt(is));
    nsamples = sampleNames.length;

    int len = readInt(is);   
    byte meta[] = new byte[len];
    is.read(meta);

    StringBuffer buff = new StringBuffer();
    for(int i=0; i<meta.length; i++)
    {
      // expecting null terminated
      if(meta[i] == 0 && i == meta.length-1)
        continue;
      buff.append((char)meta[i]);
    }

    metaData = buff.toString();
  }
  
  protected String headerToString()
  {
    StringBuffer head = new StringBuffer();
    head.append("##fileformat=VCFv4.0\n");
    head.append(metaData);
    head.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
    
    for(int i=0; i<sampleNames.length; i++)
      head.append(sampleNames[i]+" ");
    return head.toString();
  }
  
  /**
   * Given the length of the concatenated names, that are NULL padded
   * construct a list of the Strings.
   * @param len
   * @return
   * @throws IOException
   */
  private String[] getArray(int len) throws IOException
  {
    byte b[] = new byte[len];
    is.read(b);
    
    List<String> names = new Vector<String>();
    StringBuffer buff = new StringBuffer();
    for(int i=0; i<b.length; i++)
    {
      if(b[i] != 0)
        buff.append((char)b[i]);
      else if(buff.length() > 0)
      {
        names.add(buff.toString());
        buff = new StringBuffer();
      }
    }
    
    String[] arr = new String[names.size()];
    for(int i=0; i< arr.length; i++)
      arr[i] = names.get(i);
    return arr;
  }
  
  private VCFRecord readVCFRecord() throws IOException
  {
    VCFRecord bcfRecord = new VCFRecord();
    bcfRecord.setChrom( seqNames[readInt(is)] );  
    bcfRecord.setPos ( readInt(is)+1 );
    bcfRecord.setQuality( readFloat(is) );

    int slen = readInt(is);
    if(slen > 5000)
      return bcfRecord;
    byte[] str = new byte[slen];
    is.read(str);

    getParts(str, bcfRecord);

    if(formatPattern.matcher(bcfRecord.getFormat()).matches())
    {
      int n_alleles;
      if(bcfRecord.getAlt().toString().equals(":")) // non-variant
      {
        n_alleles = 1;
        bcfRecord.setAlt(".");
      }
      else
        n_alleles = bcfRecord.getAlt().getNumAlleles();
      int nc  = (int) (n_alleles * ((float)(((float)n_alleles+1.f)/2.f)));

      //String fmts[] = VCFRecord.COLON_PATTERN.split( bcfRecord.getFormat() );
      
      final int nfmt = VCFRecord.countOccurrences(bcfRecord.getFormat(), ':')+1;
      final String fmts[] = VCFRecord.split(bcfRecord.getFormat() , ":", nfmt);
      
      bcfRecord.setGenoTypeData( new String[nsamples][fmts.length] );

      for(int i=0; i<fmts.length; i++)
      {
        int nb = getByteSize(fmts[i],nc);
        str = new byte[nb];
        is.read(str);

        if(fmts[i].equals("GT"))
        {
          for(int j=0; j<nsamples; j++)
            bcfRecord.getGenoTypeData()[j][i] = getGTString(str[j]);
        }
        else if(fmts[i].equals("PL"))
        {
          String pls[] = getPLString(str, nsamples);
          for(int j=0; j<nsamples; j++)
            bcfRecord.getGenoTypeData()[j][i] = pls[j];
        }
        else if(fmts[i].equals("DP")||fmts[i].equals("SP")||fmts[i].equals("GQ"))
        {
          for(int j=0; j<nsamples; j++)
            bcfRecord.getGenoTypeData()[j][i] = Integer.toString(byteToInt(str[j]));
        }
        else
        {
          bcfRecord.getGenoTypeData()[0][i] = new String(str);
        }
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
  private void getParts(final byte[] b, final VCFRecord bcfRecord)
  {
    int index = 0;
    final StringBuilder buff = new StringBuilder();
    
    for(int i=0; i<b.length; i++)
    {
      if(b[i] == 0 && (i == 0 || b[i-1] == 0) )
      {
        switch( index )
        {
          case 0:  bcfRecord.setID("."); break;
          case 1:  bcfRecord.setRef("."); break;
          case 2:  bcfRecord.setAlt(":"); break; // this looks like a non-variant site
          case 3:  bcfRecord.setFilter("."); break;
          case 4:  bcfRecord.setInfo("."); break;
          case 5:  bcfRecord.setFormat("."); break;
          default: return;
        }
        index++;
      }
      if(b[i] == 0)
        continue;
      
      buff.setLength(0);
      buff.append((char)b[i]);
      while(i < b.length-1 && b[i+1] != 0)
        buff.append((char)b[++i]);

      switch( index )
      {
        case 0:  bcfRecord.setID(buff.toString()); break;
        case 1:  bcfRecord.setRef(buff.toString()); break;
        case 2:  bcfRecord.setAlt(buff.toString()); break;
        case 3:  bcfRecord.setFilter(buff.toString()); break;
        case 4:  bcfRecord.setInfo(buff.toString()); break;
        case 5:  bcfRecord.setFormat(buff.toString()); break;
        default: return;
      }
      index++;
    }
  }
  

  /*private void getParts(final byte[] b, final VCFRecord bcfRecord)
  {
    final StringBuilder buff = new StringBuilder();
    for(int i=0; i<b.length; i++)
    {
      if(i == 0 && b[i] == 0)
        buff.append(". ");
      else if(b[i] == 0 && b[i-1] == 0)  // this looks like a non-variant site
        buff.append(" : ");
      else if(b[i] == 0)
        buff.append(" ");
      else
        buff.append((char)b[i]);
    }

    int ind = 0;
    int lastInd = 0;
    int i = 0;
    final int len = buff.length();
    
    while((ind = buff.indexOf(" ", ind)) > -1)
    {
      switch( i )
      {
        case 0:  bcfRecord.setID(buff.substring(lastInd, ind)); break;
        case 1:  bcfRecord.setRef(buff.substring(lastInd, ind)); break;
        case 2:  bcfRecord.setAlt(buff.substring(lastInd, ind)); break;
        case 3:  bcfRecord.setFilter(buff.substring(lastInd, ind)); break;
        case 4:  bcfRecord.setInfo(buff.substring(lastInd, ind)); break;
        case 5:  bcfRecord.setFormat(buff.substring(lastInd, ind)); break;
        default: return;
      }

      ind++;
      if(ind < len && buff.charAt(ind) == ' ')
        ind++;
      lastInd = ind;
      i++;
    }
  }*/
  
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
    {
      if(newBCF)                // type changed in bcftools
        return 4*nsamples;      // uint32_t[n]
      else
        return nsamples;        // uint8_t[n]
    }
    else                        // misc
      return 4*nsamples;        // uint32_t+char*
  }

  /**
   * Phred scaled likelihood of data.
   * @param b
   * @param nsamples
   * @return
   */
  private String[] getPLString(byte[] b, int nsamples)
  {
    String pls[] = new String[nsamples];
    int nb = b.length / nsamples;
    int cnt = 0;

    for(int i=0; i<nsamples; i++)
    {
      StringBuffer buff = new StringBuffer();
      for(int j=0;j<nb; j++)
      {
        buff.append(byteToInt(b[cnt]));
        if(j<nb-1)
          buff.append(",");
        cnt++;
      }
      pls[i] = buff.toString();
    }
    return pls;
  }
  
  /**
   * GT genotype, allele values separated by / or |, i.e.
   * unphased or phased.
   * @param b
   * @return
   */
  private String getGTString(byte b)
  {
    return ((b >> 3 & 7) + ( ((b >> 6 & 1 )== 1 ) ? "|" : "/") + (b & 7));
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
  
  protected int getSeqIndex(String chr)
  {
    for(int i=0; i<seqNames.length; i++)
    {
      if(seqNames[i].equals(chr))
        return i;
    }
        
    return -1;
  }
  
  protected long queryIndex(int tid, int beg)
  {
    long min_off = -1;
    if (beg < 0) 
      beg = 0;
   
    long offset[] = idx.get(tid).index2_offset;
    int i;

    try
    {
      for(i = beg>>TAD_LIDX_SHIFT; i < idx.get(tid).n && offset[i] == 0; ++i);
      min_off = (i == idx.get(tid).n)? offset[idx.get(tid).n-1] : offset[i];
    }
    catch(ArrayIndexOutOfBoundsException e)
    {
      return offset[offset.length-1];
    }
    return min_off;
  }

  protected String getMetaData()
  {
    return metaData;
  }

  protected String[] getSeqNames()
  {
    return seqNames;
  }
  
  protected String getFileName()
  {
    return fileName;
  }
  
  protected BCFReaderIterator query(String chr, int sbeg, int send) throws IOException
  {
    return new BCFReaderIterator(chr, sbeg, send);
  }
  
  protected class BCFReaderIterator
  {
    private String chr;
    private int sbeg;
    private int send;
    private int count = 0;
    
    protected BCFReaderIterator(String chr, int sbeg, int send)
    {
      this.chr = chr;
      this.sbeg = sbeg;
      this.send = send;
    }
    
    private boolean seekPosition() throws IOException
    {
      int bid = getSeqIndex(chr);
      if(bid < 0)
      {
        VCFview.logger4j.debug(chr+" NOT FOUND");
        return false;
      }
      long off = queryIndex(bid, sbeg);
      seek(off);
      return true;
    }
    
    public VCFRecord next() throws IOException
    {
      if(count == 0 && !seekPosition())
        return null;

      count+=1;
      return nextRecord(chr, sbeg, send);
    }
  }
  
  
  public static void main(String args[])
  {
    try
    {
      int sbeg = 0;
      int send = Integer.MAX_VALUE;
      String chr = null;
      if(args.length > 1)
      {
        String parts[] = args[1].split(":");
        chr = parts[0];
        
        String rgn[] = parts[1].split("-");
        sbeg = Integer.parseInt(rgn[0]);
        send = Integer.parseInt(rgn[1]);
      }
      
      BCFReader reader = new BCFReader(args[0]);
      int bid = 0;
      if(chr != null)
        bid = reader.getSeqIndex(chr);
      
      long off = reader.queryIndex(bid, sbeg);
      reader.seek(off);

      System.out.println(reader.headerToString());
      VCFRecord bcfRecord;
      while( (bcfRecord = reader.nextRecord(chr, sbeg, send)) != null )
      {
        if(chr != null && bcfRecord.getChrom().equals(chr))
          System.out.println(bcfRecord.toString());
        else
          break;
      }

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