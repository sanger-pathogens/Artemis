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
  
  public void query(File bcf, long offset, int beg, int end) throws IOException
  {
    BlockCompressedInputStream is = new BlockCompressedInputStream(bcf);
    is.seek(offset);
    
    for(int i=0; i<7299; i++)
    {
      int id = readInt(is);
      int pos = readInt(is)+1;
      float qual = readFloat(is);
      int slen = readInt(is);
      System.out.print("ID: "+id + " POS: "+pos + " QUAL: "+qual + " ");
      byte[] str = new byte[slen];
      is.read(str);

      String parts[] = getParts(str);
      
      String ref = parts[0];
      String alt = parts[1];
      String fmt = parts[parts.length-1];

      if(formatPattern.matcher(fmt).matches())
      {
        String info = parts[parts.length-2];
        System.out.println(info+"    "+fmt);
        
        String format = parts[parts.length-1];
        
        int nc  = 3;
        System.out.println("ALT:"+alt+" REF:"+ref);
        if(alt.equals("."))
          nc = 1;

        String fmts[] = format.split(":");
        for(int j=0; j<fmts.length; j++)
        {
          int nb = getByteSize(fmts[j],1,nc);
          str = new byte[nb];
          is.read(str);
          
          if(fmts[j].equals("GT"))
            System.out.println("nbytes = "+nb+" GT:"+getGTString(str[0]));
          else if(fmts[j].equals("PL"))
            System.out.println("nbytes = "+nb+" PL:"+getPLString(str, nc));
          else if(fmts[j].equals("DP")||fmts[j].equals("SP")||fmts[j].equals("GQ"))
          {
            System.out.println("nbytes = "+nb+" "+fmts[j]+":"+byteToInt(str[0]));
          } 

        }
      }
    }
  }
  
  /**
   * Make a string from the byte array. Expanding NULL padding.
   * @param b
   * @return
   */
  private String[] getParts(byte[] b)
  {
    StringBuffer buff = new StringBuffer();
    for(int i=0; i<b.length; i++)
    {
      if(i == 0 && b[i] == 0)
        continue;

      if(b[i] == 0)
      {
        if(i > 0 && b[i-1] == 0 && b[i+1] == 0)
          buff.append(".");
        else
          buff.append(" ");
      }
      else
        buff.append((char)b[i]);
    }
    
    return buff.toString().split(" ");
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
  private static int getByteSize(String tag, int nsamples, int nc)
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
  
  public static List<Index> load(File bcfIndex) throws IOException
  {
    FileInputStream fis = new FileInputStream(bcfIndex);
    BlockCompressedInputStream is = new BlockCompressedInputStream(fis);
    byte[] magic = new byte[4];
    is.read(magic);
    System.out.println(new String(magic));
    
    int n = readInt(is);
    
    List<Index> idx = new Vector<Index>(n);
    
    for(int i=0; i<n; i++)
    {
      Index idx2 = new Index();
      idx2.n = readInt(is);
      idx2.index2_offset = new long[idx2.n];
      
      for(int j=0; j<idx2.n; j++)
      {
        idx2.index2_offset[j] = readLong(is);
      }
      
      if(is.read() == -1)
        System.out.println("EOF");
      idx.add(idx2);
    }
    return idx;
  }
  
  public long queryIndex(List<Index> idx, int tid, int beg)
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
  
  public static int readInt(final InputStream is) throws IOException {
    byte[] buf = new byte[4];
    is.read(buf);
    return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getInt();
  }
  
  
  public static float readFloat(final InputStream is) throws IOException {
    byte[] buf = new byte[4];
    is.read(buf);
    return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getFloat();
  }

  public static long readLong(final InputStream is) throws IOException {
    byte[] buf = new byte[8];
    is.read(buf);
    return ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN).getLong();
  }

  public static void main(String args[])
  {
    try
    {
      List<Index> idx = load(new File(args[0]));
      
      int sbeg;
      int send;
      if(args.length < 3)
      {
        sbeg = 326758;
        send = sbeg+1;
      }
      else
      {
        sbeg = Integer.parseInt(args[2]);
        send = Integer.parseInt(args[3]);
      }
      
      BCFReader reader = new BCFReader();
      long off = reader.queryIndex(idx, 0, sbeg);
      System.out.println(off);
      reader.query(new File(args[1]), off, sbeg, send);
    }
    catch (IOException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
}

class BCFRecord
{
  int seqID;
  int pos;
  float quality;
  String info;
  String format;
}

class Index
{
  int n;
  long index2_offset[];
}