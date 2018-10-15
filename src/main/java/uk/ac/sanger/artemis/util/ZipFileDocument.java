/* ZipFileDocument
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2010  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.util;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

public class ZipFileDocument extends FileDocument
{
  private File zipFile;
  private String zipEntryName;
  private byte b[] = null;
  
  public ZipFileDocument(File zipFile, String zipEntryName)
  {
    super(zipFile);
    this.zipFile = zipFile;
    this.zipEntryName = zipEntryName;
  }
  
  /**
   *  Return the name of this Document.
   **/
  public String getName () 
  {
    return zipEntryName;
  }
  
  /**
   *  Create a new InputStream object from this Document.  The contents of the
   *  Document can be read from the InputStream.
   *  @exception IOException Thrown if the Document can't be read from
   *    (for example if it doesn't exist).
   **/
  public InputStream getInputStream() throws IOException
  {
    if(!zipFile.getName().endsWith(".zip"))
      return super.getInputStream();

    if(b == null)
      b = getEntryContent(zipFile, zipEntryName);
    
    if(b == null)
    {
      if(!zipEntryName.endsWith(".gz"))
      {
        zipEntryName = zipEntryName + ".gz";
        b = getEntryContent(zipFile, zipEntryName);
      }
      
      if(b == null)
        return null;
    }
    if(zipEntryName.endsWith(".gz"))
      return new WorkingGZIPInputStream(new ByteArrayInputStream(b));
    return new ByteArrayInputStream(b);
  }
  
  public String writeTmpFile(String txt) throws IOException
  {
    final File tmpFile = File.createTempFile(zipEntryName, "tmp");
    tmpFile.deleteOnExit();
  
    FileWriter out = new FileWriter(tmpFile);
    out.write(txt);
    out.close();
    return tmpFile.getAbsolutePath();
  }
  
  /**
   * Returns true if the zip file contains the entry name.
   * @return
   */
  private boolean containsEntry()
  {
    try
    {
      FileInputStream fis = new FileInputStream(zipFile);
      ZipInputStream zis = new ZipInputStream(
          new BufferedInputStream(fis));
      ZipEntry ze;
      
      while ((ze=zis.getNextEntry())!=null)
      {
        if( (ze.getName().equals(zipEntryName) || ze.getName().equals(zipEntryName+".gz")) &&
             !ze.isDirectory())
        {
          zipEntryName = ze.getName();
          
          ByteBuffer buff = new ByteBuffer();
          b = new byte[1];
          while(zis.read(b,0,1) != -1)
            buff.append(b);
          
          b = buff.getBytes();
          
          fis.close();
          zis.close();
          return true;
        }
      }
      fis.close();
      zis.close();
    }
    catch (IOException e){}
    
    return false;
  }
  
  /**
   *  Return true if and only if the Document referred to by this object exists
   *  and is readable.
   **/
  public boolean readable() 
  {
    if(getFile().exists() && 
       getFile ().canRead() && 
       zipFile.getName().endsWith(".zip")) 
    {
      if(containsEntry())
        return true;
    } 
    return false;
  }
  
  private static byte[] getEntryContent(File zipFile, String zipEntryName) throws IOException
  {
    ZipInputStream zis= new ZipInputStream(
                     new FileInputStream(zipFile));
      
    ZipEntry ze;
    while ((ze=zis.getNextEntry())!=null)
    {
      if(ze.isDirectory() || !ze.getName().equals(zipEntryName))
        continue;

      ByteBuffer buff = new ByteBuffer();
      byte[] b = new byte[1];
      while(zis.read(b,0,1) != -1)
        buff.append(b);
      return buff.getBytes();
    }
    return null;
  }

  public static void main(String[] args)
  {
    System.out.println(args[1]);

    try
    {
      ZipFileDocument document = new ZipFileDocument(new File(args[0]), args[1]);
      InputStream in = document.getInputStream();
      StringBuffer out = new StringBuffer();
      byte[] b = new byte[4096];
      for (int n; (n = in.read(b)) != -1;) {
          out.append(new String(b, 0, n));
      }
      System.out.println(new String(out));
      
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
}