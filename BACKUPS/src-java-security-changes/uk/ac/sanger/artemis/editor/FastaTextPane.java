/*
 *
 * created: Wed Aug 3 2004
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
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

package uk.ac.sanger.artemis.editor;

import java.awt.Point;
import javax.swing.JTextArea;
import javax.swing.JScrollPane;

import com.sshtools.sftp.SftpFileAttributes;

import java.awt.Font;
import java.awt.Dimension;
import java.io.InputStream;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.util.Vector;
import java.util.Enumeration;
import java.util.StringTokenizer;
import java.net.URL;
import java.net.MalformedURLException;

import uk.ac.sanger.artemis.components.filetree.FileList;
import uk.ac.sanger.artemis.util.Document;

public class FastaTextPane extends JScrollPane
{
  /***/
  private static final long serialVersionUID = 1L;
  private JTextArea textArea;
  private Vector<HitInfo> hitInfoCollection = null;
  private String format = null;
  private Document document;
  private int qlen;
  private Vector<FastaListener> listeners = new Vector<FastaListener>();
  private Vector<GetzThread> threads = new Vector<GetzThread>();
  protected static int MAX_HITS = 70;
  private static boolean remoteMfetch = false;
  private static boolean forceUrl = false;
  private boolean blastPlus = false;
  
  static 
  {
    if(System.getProperty("useSRS") != null)
      forceUrl = true;
  }

  public static HitInfo[] cacheHits = new HitInfo[BigPane.CACHE_SIZE];
  public static int nCacheHits = 0;
  private static String mfetchExecutablePath = "/nfs/disk100/pubseq/bin/mfetch";
  private static File mfetchExecutable = new File(mfetchExecutablePath);
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(FastaTextPane.class);
  
  public FastaTextPane(Document document)
  {
    super();

    this.document = document;
    format = getResultsFormat(); 

    if(format != null)
    {
      StringBuffer contents = null;
      if(format.equals("fasta"))
        contents = readFASTAFile(format);
      else if(format.equals("blastp"))
        contents = readBLASTPFile(format);
      textArea = new JTextArea(contents.toString());
      setTextAreaFont(BigPane.font);
      textArea.setEditable(false);

      setViewportView(textArea);
      setPreferredSize(new Dimension(500,300));
    }
  }

  protected void addFastaListener(FastaListener obj)
  {
    listeners.add(obj);
  }

  protected void reRead()
  {
    StringBuffer contents = null;

    format = getResultsFormat();
    if(format.equals("fasta"))
      contents = readFASTAFile(format);
    else if(format.equals("blastp"))
      contents = readBLASTPFile(format);

    textArea.setText(contents.toString());
    setViewportView(textArea);

    Enumeration<FastaListener> enumListeners = listeners.elements();
    while(enumListeners.hasMoreElements())
      enumListeners.nextElement().update();
  }

  /**
  *
  * Get the format of the results (e.g. FASTA or BLASTP).
  * @return format
  *
  */
  protected String getFormat()
  {
    return format;
  }

  private void setTextAreaFont(Font f)
  {
    textArea.setFont(f);
  }

  /**
  *
  * Get the format of the results in a file. FASTA and
  * BLASTP are supported.
  *
  */
  private String getResultsFormat()
  {
    InputStreamReader streamReader = null;
    BufferedReader buffReader = null;
    String line = null;
    String format = null;

    try
    {
      streamReader = new InputStreamReader(document.getInputStream());
      buffReader = new BufferedReader(streamReader);
      while( (line = buffReader.readLine()) != null)
      {
        if(line.startsWith("BLASTP"))
        {
          format = "blastp";
          if(line.indexOf("+") > -1)
            blastPlus = true;
          break;
        }
        else if(line.indexOf("FASTA") > -1)
        {
          format = "fasta";
          break;
        }
      }
      streamReader.close();
      buffReader.close();
    }
    catch(IOException ioe)
    {
      System.out.println("Cannot read file: " + 
          document.getName() +"\n"+ioe.getMessage());
    }
    
    return format;
  }


  private StringBuffer readBLASTPFile(final String format)
  {
    //File fn = dataFile;
    StringBuffer sbuff = new StringBuffer();

    InputStreamReader streamReader = null;
    BufferedReader buffReader = null;

    hitInfoCollection = new Vector<HitInfo>();
    try
    {
      streamReader = new InputStreamReader(document.getInputStream());
      buffReader = new BufferedReader(streamReader);

      String line = null;
      int textPosition = 0;
      int len     = 0;
      HitInfo hit = null;
      int ind1 = 0;
      //int endQuery = 0;

      while( (line = buffReader.readLine()) != null)
      {
        len = line.length()+1;
        sbuff.append(line+"\n");
        if(line.startsWith("Sequences producing significant alignments:"))
        {
          buffReader.readLine();
          while( !(line = buffReader.readLine()).equals("") )
          {
            textPosition += line.length()+1;
            sbuff.append(line+"\n");
  
            hit = new HitInfo(line,format);

            hitInfoCollection.add(hit);
          }
        }
        else if(line.indexOf(" ----") > -1 ||
                line.indexOf(" --- ") > -1 ||
                line.indexOf(" -- ") > -1 )
        {
        }
        else if(line.startsWith(">"))  // start of alignment
        {
          String currentID = line;

          int ind = line.indexOf(" ");
          if(ind > -1)
          {
            currentID = line.substring(1,ind);
            int indDot;
            if( (indDot = currentID.indexOf(".")) > 5)
            {
              String version = currentID.substring(indDot+1);
              try
              {
                int version_no = Integer.parseInt(version);
                // version number looks like uniprot so strip this
                if(version_no < 50)
                  currentID= currentID.substring(0,indDot);
              }
              catch(NumberFormatException nfe){}
            }
          }

          int ind2 = currentID.indexOf(":");
          if(ind2 > -1)
          {
            currentID = currentID.substring(ind2+1); 
          }

          if(hit != null)
            hit.setEndPosition(textPosition);
          
          hit = getHitInfo(currentID,hitInfoCollection);
          hit.setStartPosition(textPosition);

          String going = "";
          ind = line.indexOf("GO:");
          if(ind > -1)
            going = line.substring(ind+3);
          
          String nextLine = null;
//        buffReader.mark(210);

// get GO numbers
          while((nextLine = buffReader.readLine()).indexOf("Length") == -1)
          {
            len += nextLine.length()+1;
            sbuff.append(nextLine+"\n");
            if(going.equals("") && ((ind = nextLine.indexOf("GO:")) > -1))
              going = nextLine.substring(ind+3);
            else if(!going.equals(""))
              going = going.concat(nextLine);
          }
          
          if(!going.equals(""))
            hit.setGO(going);

          if(nextLine != null)
          {
            len += nextLine.length()+1;
            sbuff.append(nextLine+"\n");
            if( (ind1 = nextLine.indexOf("  Length = ")) > -1)
              hit.setLength(nextLine.substring(ind1+11));
          }

// get query start
          int start = 999999;
          int end   = 0;
          boolean seen = false;

          while((nextLine = buffReader.readLine()) != null &&
                !nextLine.startsWith(">"))
          {
            len += nextLine.length()+1;
            sbuff.append(nextLine+"\n");

            if(nextLine.startsWith(" Score ="))
            {
//            System.out.println("start:: "+start+"  end:: "+end);
              if(end != 0)
                hit.setQueryPosition(start,end);
              start = 999999;
              end   = 0;
              seen  = true;
            }
            else if(nextLine.startsWith("Query:"))
            {
              ind1 = nextLine.indexOf(" ",8);
              int nstart = Integer.parseInt(nextLine.substring(7,ind1).trim());
              if(nstart < start)
                start = nstart;
              end  = Integer.parseInt(nextLine.substring(nextLine.lastIndexOf(" ")).trim());
              seen = false;
            }
            buffReader.mark(100); 
          }
          if(!seen && end != 0)
            hit.setQueryPosition(start,end);

          buffReader.reset();
        }
        else if( (ind1 = line.indexOf("Identities = ")) > -1)
        {
          ind1 = line.indexOf("(",ind1)+1;
          if(ind1 > -1)
            hit.setIdentity(line.substring(ind1,line.indexOf(")",ind1)).trim());
        }
        else if( (ind1 = line.indexOf("  Length = ")) > -1)
          hit.setLength(line.substring(ind1+11));
//      else if(line.startsWith("Query: "))
//      {
//        hit.setQueryEnd(Integer.parseInt(line.substring(line.lastIndexOf(" ")).trim()));
//      }
        else if(line.startsWith("Query="))
        {
          int ind2 = 0;
          if(blastPlus)
          {
            ind2 = line.indexOf("Length=");
            if(ind2 == -1)
            {
              while((line = buffReader.readLine()).indexOf("Length=") < 0)
              {
                len += line.length()+1;
                sbuff.append(line+"\n");
              }
              ind2 = line.indexOf("Length=");
            }
            ind1 = line.length();
            qlen = Integer.parseInt(line.substring(ind2+7,ind1).trim());
          }
          else
          {
            ind1 = line.indexOf(" letters)");
            if(ind1 == -1)
            {
              String nextLine = null;
              while((nextLine = buffReader.readLine()).indexOf(" letters)") < 0)
              {
                len += nextLine.length()+1;
                sbuff.append(nextLine+"\n");
              }
              line = nextLine;
              ind1 = nextLine.indexOf(" letters)");
              ind2 = nextLine.indexOf("(");
            }
            else
              ind2 = line.indexOf("(");

            qlen = Integer.parseInt(line.substring(ind2+1,ind1).trim());
          }
        }
        textPosition += len;
      }

      if(hit != null)
        hit.setEndPosition(textPosition);

      streamReader.close();
      buffReader.close();

      GetzThread getz = new GetzThread(hitInfoCollection);
      getz.start();
      threads.add(getz);
    }
    catch (IOException ioe)
    {
      System.out.println("Cannot read file: " + document.getName());
    }
    return sbuff;
  }


  private StringBuffer readFASTAFile(final String format)
  {
    //File fn = dataFile;
    StringBuffer sbuff = new StringBuffer();

    InputStreamReader streamReader = null;
    BufferedReader buffReader = null;

    hitInfoCollection = new Vector<HitInfo>();
    try
    {
      streamReader = new InputStreamReader(document.getInputStream());
      buffReader = new BufferedReader(streamReader);

      String line = null;
      int textPosition = 0;
      int len    = 0;
      HitInfo hi = null;

      while( (line = buffReader.readLine()) != null)
      {
        len = line.length()+1;
        sbuff.append(line+"\n");  

        int ind1;

        if(line.endsWith(" aa"))
        {
          String tmp = line.substring(0,line.length()-3).trim();
          int in1 = tmp.lastIndexOf(" ");
          if(in1 > -1)
            qlen = Integer.parseInt(tmp.substring(in1).trim());
        }
        else if(line.startsWith("The best scores are:"))
        {
          while( !(line = buffReader.readLine()).equals("") )
          {
            textPosition += line.length()+1;
            sbuff.append(line+"\n");
            hitInfoCollection.add(new HitInfo(line,format)); 
          }
        }
        else if(line.startsWith(">>"))  // start of alignment
        {
          int ind = line.indexOf(" ");
          String currentID = line.substring(2,ind);
 
          int indDot;
          if( (indDot = currentID.indexOf(".")) > 5)
          {
            String version = currentID.substring(indDot+1);
            try
            {
              int version_no = Integer.parseInt(version);
              // version number looks like uniprot so strip this
              if(version_no < 50)
                currentID= currentID.substring(0,indDot);
            }
            catch(NumberFormatException nfe){}

            // HERE currentID= currentID.substring(0,indDot);
          }

          if(hi != null)
            hi.setEndPosition(textPosition);

          hi = getHitInfo(currentID,hitInfoCollection);
          hi.setStartPosition(textPosition);
        }
        else if(line.startsWith("Smith-Waterman")) // Smith-Waterman
        {
          ind1 = line.indexOf("score:");
          int ind2;
          if(ind1 > -1)
          {
            ind2 = line.indexOf(";",ind1);

            hi.setScore(line.substring(ind1+6,ind2));
     
            ind1 = ind2+1;
            ind2 = line.indexOf("identity");
            if(ind2 > -1)
              hi.setIdentity(line.substring(ind1,ind2).trim());
          
            ind1 = line.indexOf("(",ind2);
            if(ind1 > -1)
            {
              ind2 = line.indexOf("ungapped)",ind1);
              
              // fasta34 changed to similar
              if(ind2 == -1)
                ind2 = line.indexOf("similar)",ind1);
              
              if(ind2 > -1)
                hi.setUngapped(line.substring(ind1+1,ind2).trim());
            }

            ind1 = line.indexOf(" in ",ind2);
            ind2 = line.indexOf("(",ind1);
            if(ind1 > -1 && ind2 > -1)
              hi.setOverlap(line.substring(ind1+4,ind2).trim());
           
            ind1 = ind2+1;
            ind2 = line.indexOf(":",ind1);
            if(ind2 > -1)
            {
              String range = line.substring(ind1,ind2);
              hi.setQueryRange(range);
              int split = range.indexOf("-");
              if(split > -1)
              {   
//              hi.setQueryStart(Integer.parseInt(range.substring(0,split)));
//              hi.setQueryEnd(Integer.parseInt(range.substring(split+1)));
                hi.setQueryPosition(Integer.parseInt(range.substring(0,split)),
                                    Integer.parseInt(range.substring(split+1)));
              }
            }

            ind1 = ind2+1;
            ind2 = line.indexOf(")",ind1);
            if(ind2 > -1)
              hi.setSubjectRange(line.substring(ind1,ind2)); 
          }
        }
        else if( (ind1 = line.indexOf(" E():")) > -1)
        {
          StringTokenizer tok = new StringTokenizer(line.substring(ind1+5));
          hi.setEValue(tok.nextToken().trim());
        }
 
        textPosition += len;
      }
  
      if(hi != null)
        hi.setEndPosition(textPosition);
   
      streamReader.close();
      buffReader.close();

      GetzThread getz = new GetzThread(hitInfoCollection);
      getz.start();
      threads.add(getz);
    }
    catch (IOException ioe)
    {
      System.out.println("Cannot read file: " + document.getName());
    }
    return sbuff;
  }

  /**
  *
  * Get the query sequence length.
  * @return query sequence length.
  *
  */
  protected int getQueryLength()
  {
    return qlen;
  }

  protected Vector<HitInfo> getHitCollection()
  {
    return hitInfoCollection;
  }


  private HitInfo getHitInfo(String acc, final Vector<HitInfo> hits)
  {
    int ind = 0;
    acc     = acc.trim();
    
    if((ind = acc.indexOf(";")) > -1)
      acc = acc.substring(0,ind);

    HitInfo hit = null;
    for(int i=0;i<hits.size(); i++)
    {
      hit = (HitInfo)hits.get(i);
      if(hit.getAcc().equals(acc) ||
         hit.getID().equals(acc))
      {
        return hit;
      }
    }

    return null;
  }

   
  /**
  *
  * Stop all getz processes
  *
  */
  protected void stopGetz()
  {
    Enumeration<GetzThread> threadEnum = threads.elements();

    while(threadEnum.hasMoreElements())
    {
      GetzThread gthread = threadEnum.nextElement();
      if(gthread.isAlive())
        gthread.stopMe();
    }
  }

  class GetzThread extends Thread
  {
    private Vector<HitInfo> hitInfoCollection;
    private boolean keepRunning = true;

    protected GetzThread(Vector<HitInfo> hitInfoCollection)
    {
      this.hitInfoCollection = hitInfoCollection;
    }

    public void run()
    {
      getzCall(hitInfoCollection,hitInfoCollection.size());
    }

    protected void stopMe()
    {
      keepRunning = false;
    }

    /**
    *
    * Creates and executes an SRS query for all the hits in the
    * collection.
    * @param hit          HitInfo for a single hit.
    * @param ortholog     true if ortholog is selected.
    *
    */
    private void getzCall(final Vector<HitInfo> hits, final int nretrieve)
    {
      final String env[] = { "PATH=/usr/local/pubseq/bin/:/nfs/disk100/pubseq/bin/" };

      // split mfetch query up - max 70 hits per query
      int nhits = hits.size()/MAX_HITS + 1;
      StringBuffer querySRS = new StringBuffer();
      StringBuffer queryMfetch[] = new StringBuffer[nhits];
      
      int n = 0;
      int this_nhit = 0;
      nhits = 0;
      
      if(!getMfetchExecutable().exists() && !isRemoteMfetch() && !isForceUrl() &&
          System.getProperty("j2ssh") != null &&
          FileList.isConnected())
      {
        FileList fileList = new FileList();
        SftpFileAttributes attr = fileList.stat(mfetchExecutablePath);
        if(attr != null)
          remoteMfetch = attr.isFile();
      }
        
      Enumeration<HitInfo> ehits = hits.elements();
      while(ehits.hasMoreElements() && keepRunning)
      {
        if(n>nretrieve)
          break;
        HitInfo hit = ehits.nextElement();
        
        HitInfo cacheHit = checkCache(hit);
        if(cacheHit != null)
        {
          logger4j.debug("Retrieved early from cache "+cacheHit.getID()+
                         " cache size="+cacheHits.length);
          hit.setOrganism(cacheHit.getOrganism());
          hit.setDescription(cacheHit.getDescription());
          hit.setGeneName(cacheHit.getGeneName());
          if(cacheHit.getEMBL() != null)
            hit.setEMBL(cacheHit.getEMBL());
          continue;
        }
        
        if(n > 0)
        {
          querySRS.append("|");
          
          if(nhits > 0)
            queryMfetch[this_nhit].append("|acc:");
        }

        querySRS.append(hit.getAcc());
        
        if(queryMfetch[this_nhit] == null)
          queryMfetch[this_nhit] = new StringBuffer();
        
        queryMfetch[this_nhit].append(hit.getAcc());
        n++;
        nhits++;
        if(nhits > 70)
        {
          nhits = 0;
          this_nhit++;
        }
      }
      
      if(!keepRunning)
        return;
      
      boolean isLocalMfetchExists = getMfetchExecutable().exists();
      if(isForceUrl())
        isLocalMfetchExists = false;

      BufferedReader strbuff = null;
      final File fgetz = new File("/usr/local/pubseq/bin/getz");
      
      if(isLocalMfetchExists || remoteMfetch)
      {
        StringBuffer buff = new StringBuffer();
        for(int i=0; i<queryMfetch.length; i++)
        {
          if(queryMfetch[i] == null ||
             queryMfetch[i].toString().length() == 0 )
            continue;
          
          String mfetch = queryMfetch[i].toString();
          
          if(isLocalMfetchExists)  // local
          {
            String cmd[]   = { "mfetch", "-f", "acc org des gen",
              "-d", "uniprot", "-i", "acc:"+mfetch };
       
            ExternalApplication app = new ExternalApplication(cmd,
                                      env,null);
            buff.append(app.getProcessStdout());
          }
          else                 // remote
          {
            String cmd   = 
              "mfetch -f \"acc org des gen\" -d uniprot -i \"acc:"+mfetch+"\"" ;
            uk.ac.sanger.artemis.j2ssh.SshPSUClient ssh =
              new uk.ac.sanger.artemis.j2ssh.SshPSUClient(cmd);
            ssh.run();
            buff.append(ssh.getStdOutBuffer());
          }
        }
        
        final StringReader strread   = new StringReader(buff.toString());
        strbuff = new BufferedReader(strread);
      }
      else if(!fgetz.exists())
      {
        try
        {
          String queryURL = DataCollectionPane.srs_url+
                            "/wgetz?-noSession+-ascii+-f+acc%20org%20description%20gen+[uniprot-acc:"+
                            querySRS.toString()+"]+-lv+500";
          URL wgetz = new URL(queryURL);
          logger4j.debug(queryURL);
          
          InputStream in = wgetz.openStream();

          strbuff = new BufferedReader(new InputStreamReader(in));
          StringBuffer resBuff = new StringBuffer();
          String line;
          while((line = strbuff.readLine()) != null)
            resBuff.append(line);
   
          strbuff.close();
          in.close();

          String res = resBuff.toString();
          res= insertNewline(res, "OS ");
          res= insertNewline(res, "DE ");
          res= insertNewline(res, "GN ");
          res= insertNewline(res, "AC ");

          StringReader strread   = new StringReader(res);
          strbuff = new BufferedReader(strread);

//        System.out.println(DataCollectionPane.srs_url+
//                           "/wgetz?-f+acc%20org%20description%20gen+[uniprot-acc:"+
//                           query.toString()+"]+-lv+500");
//        System.out.println("HERE\n"+res);
        }
        catch(MalformedURLException e) {System.err.println(e);}
        catch(IOException e) {System.err.println(e);} 
      }
      else
      {
        String cmd[]   = { "getz", "-f", "acc org description gen",
                           "[uniprot-acc:"+querySRS.toString()+"]" };
                      
        ExternalApplication app = new ExternalApplication(cmd,
                                                   env,null);
        String res = app.getProcessStdout();
        StringReader strread   = new StringReader(res);
        strbuff = new BufferedReader(strread);
      }             

      HitInfo hit[] = new HitInfo[2];
      HitInfo this_hit = null;
      String line = null;
      String lineStrip = null;

      try             
      { 
        while((line = strbuff.readLine()) != null)
        { 
          line = line.trim();
          if(line.equals(""))
            continue; 
           
          if(line.length() < 3)  // empty description line
            continue;
        
          lineStrip = line.substring(3).trim();
          if(line.startsWith("AC "))
          {
            String acc1;
            String acc2 = null;
            
            if(lineStrip.endsWith(";"))
              lineStrip = lineStrip.substring(0, lineStrip.length()-1);
              
            int ind = lineStrip.indexOf(";");
            
            if(ind > -1)
            {
              acc1 = lineStrip.substring(0,ind);
              acc2 = lineStrip.substring(ind+1);
            }
            else 
              acc1 = lineStrip;
            
            hit[0] = getHitInfo(acc1,hits);
            if(acc2 != null)
              hit[1] = getHitInfo(acc2,hits);
            else
              hit[1] = null;
            
            if(hit[0] == null)
            {         
              logger4j.warn("HIT NOT FOUND "+lineStrip);
              continue;
            }         
                      
            hit[0].setGeneName("");
            
            if(hit[1] != null)
              hit[1].setGeneName("");
          }           
                
          for(int i = 0; i < hit.length; i++)
          {
            this_hit = hit[i];
            if(this_hit == null)
              continue;

            if(line.startsWith("OS "))
              this_hit.setOrganism(lineStrip);
            else if(line.startsWith("DE "))
              this_hit.appendDescription(lineStrip);
            else if(line.startsWith("GN "))
            {
              StringTokenizer tokGN = new StringTokenizer(lineStrip, ";");
              while(tokGN.hasMoreTokens())
              {
                line = tokGN.nextToken();
                if(line.startsWith("Name="))
                  this_hit.setGeneName(line.substring(5));
                // else
                // hit.appendDescription(line);
              }
            }
          }
        }

        strbuff.close();   
      }   
      catch(IOException ioe){}
      
      
      //
      //
      // mfetch
      if(isLocalMfetchExists || remoteMfetch)
      {
        getDbXRefWithMfetch(isLocalMfetchExists, queryMfetch, env, hits);
        return;
      }
      
      //
      //
      // SRS methods
      String res;
      ehits = hits.elements();
      while(ehits.hasMoreElements() && keepRunning)
      {               
        this_hit = (HitInfo)ehits.nextElement();
        
        HitInfo cacheHit = checkCache(this_hit);
        
        if(cacheHit != null && cacheHit.getEMBL() != null)
        {
          logger4j.debug("Retrieved from cache "+cacheHit.getID()+
                         " cache size="+cacheHits.length);
          this_hit.setEMBL(cacheHit.getEMBL());
          this_hit.setEC_number(cacheHit.getEC_number());
          this_hit.setGeneName(cacheHit.getGeneName());
          continue;
        }
        
        res = getUniprotLinkToDatabase(fgetz, isLocalMfetchExists, this_hit, env, "embl");
              
        int ind1 = res.indexOf("ID ");
        if(ind1 > -1) 
        {             
          StringTokenizer tok = new StringTokenizer(res);
          tok.nextToken();
          this_hit.setEMBL(tok.nextToken());
        }             
        else          
          this_hit.setEMBL("");
                      
        // EC_number  
        if(this_hit.getEC_number() != null)
          continue;   

        res = getUniprotLinkToDatabase(fgetz, isLocalMfetchExists, this_hit, env, "enzyme");

        ind1 = res.indexOf("ID ");
        if(ind1 > -1) 
        {             
          StringTokenizer tok = new StringTokenizer(res);
          tok.nextToken();
          this_hit.setEC_number(tok.nextToken());
        }
        
        addHitToCache(this_hit);
      }               
    }              
  }
  
  /**
   * Add a hit to the cache
   * @param thisHit
   */
  private void addHitToCache(final HitInfo thisHit)
  {
    if(nCacheHits >= cacheHits.length)
      nCacheHits = 0;
    cacheHits[nCacheHits] = thisHit;
    nCacheHits++;
  }
  
  /**
   * 
   * @param isLocalMfetchExists
   * @param queryMfetch
   * @param env
   * @param hits
   */
  private void getDbXRefWithMfetch(final boolean isLocalMfetchExists,
                                   final StringBuffer queryMfetch[],
                                   final String env[], 
                                   final Vector<HitInfo> hits)
  {
    String res = null;
    String line = null;

    for(int i = 0; i < queryMfetch.length; i++)
    {
      if(queryMfetch[i] == null || queryMfetch[i].toString().length() == 0)
        continue;

      final String mfetch = queryMfetch[i].toString();
      
      //
      // link to EMBL
      res = getUniprotLinkToDatabaseByMFetch(isLocalMfetchExists, mfetch, env,
          "embl");

      StringReader strread = new StringReader(res);
      BufferedReader strbuff = new BufferedReader(strread);

      try
      {
        while((line = strbuff.readLine()) != null)
        {
          String acc;
          if(line.startsWith("linked from:"))
          {
            acc = line.substring(12).trim();
            int ind = acc.indexOf('.');
            if(ind > -1)
              acc = acc.substring(0, ind);

            HitInfo thisHit = getHitInfo(acc, hits);            
            if(thisHit == null)
            {
              logger4j.warn(acc+" NOT FOUND");
              continue;
            }
            
            line = strbuff.readLine();

            int ind1 = line.indexOf("ID ");
            if(ind1 > -1)
            {
              int ind2 = line.indexOf(';');
              thisHit.setEMBL(line.substring(3, ind2).trim());
            }
            else
              thisHit.setEMBL("");
            addHitToCache(thisHit);
          }
        }
      }
      catch(IOException e)
      {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      
      final String accessions[] = mfetch.split("\\|");
      
      for(int j=0;j<accessions.length; j++)
      {
        String acc;
        if(accessions[j].startsWith("acc:"))
          acc = accessions[j].substring(4);
        else
          acc = accessions[j];
        HitInfo thisHit = getHitInfo(acc, hits);
        if(thisHit != null && thisHit.getEMBL() == null)
          thisHit.setEMBL("");
      }
      
      //
      // link to enzyme
      res = getUniprotLinkToDatabaseByMFetch(isLocalMfetchExists, mfetch, env, "enzyme");

      strread = new StringReader(res);
      strbuff = new BufferedReader(strread);

      try
      {
        while((line = strbuff.readLine()) != null)
        {
          String acc;
          if(line.startsWith("linked from:"))
          {
            acc = line.substring(12).trim();
            int ind = acc.indexOf('.');
            if(ind > -1)
              acc = acc.substring(0, ind);

            HitInfo thisHit = getHitInfo(acc, hits);
            if(thisHit == null)
            {
              logger4j.warn(acc + " NOT FOUND");
              continue;
            }

            line = strbuff.readLine();

            int ind1 = line.indexOf("ID ");
            if(ind1 > -1)
            {
              int ind2 = line.indexOf(';');
              if(ind2 < 0)
                ind2 = line.length();
              thisHit.setEC_number(line.substring(3, ind2).trim());
            }
          }
        }
        
        for(int j=0;j<accessions.length; j++)
        {
          String acc;
          if(accessions[j].startsWith("acc:"))
            acc = accessions[j].substring(4);
          else
            acc = accessions[j];
          HitInfo thisHit = getHitInfo(acc, hits);
          if(thisHit != null && thisHit.getEMBL() == null)
            thisHit.setEC_number("");
        }
      }
      catch(IOException e)
      {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
    
    return;
  }

  /**
   * 
   * @param hit
   * @return
   */
  protected static HitInfo checkCache(final HitInfo hit)
  {
    for(int i=0; i<cacheHits.length; i++)
    {
      if(cacheHits[i] != null &&
         cacheHits[i].getID().equals(hit.getID()))
        return cacheHits[i];
    }
    return null;
  }
  
  protected static String insertNewline(String s1, final String s2)
  {
    int index = 0;
    while((index = s1.indexOf(s2, index)) > -1)
    {
      if(index > 0)
        s1 = s1.substring(0,index-1)+"\n"+s1.substring(index);
      index++;
    }
    return s1;
  }

  /**
  *
  * Link Uniprot to the another database (e.g. EMBL or ENZYME)
  *
  */
  protected static String getUniprotLinkToDatabase(final File fgetz, 
                                                   final boolean isLocalMfetchExists, 
                                                   final HitInfo hit,
                                                   final String env[], final String DB)
  {
    String res = null;

    if(isLocalMfetchExists)
    {
      final String cmd[]   = { "mfetch", "-f", "id",
          "-d", "uniprot", "-i", "acc:"+hit.getID(), 
          "-l", DB };

      ExternalApplication app = new ExternalApplication(cmd,
                                  env,null);
      res = app.getProcessStdout();

    }
    else if(remoteMfetch)
    {
      final String cmd   = 
        "mfetch -f id -d uniprot -i \"acc:"+hit.getID()+"\" -l "+DB ;
      uk.ac.sanger.artemis.j2ssh.SshPSUClient ssh =
        new uk.ac.sanger.artemis.j2ssh.SshPSUClient(cmd);
      ssh.run();
      res = ssh.getStdOut();
    }
    else if(!fgetz.exists())
    {
      try
      {
        final String queryURL = DataCollectionPane.srs_url+
                                "/wgetz?-noSession+-ascii+-f+id+[uniprot-acc:"+hit.getAcc()+"]%3E"+DB;
        URL wgetz = new URL(queryURL);
        logger4j.debug(queryURL);
        
        InputStream in = wgetz.openStream();
        BufferedReader strbuff = new BufferedReader(new InputStreamReader(in));
        StringBuffer resBuff = new StringBuffer();
        String line;
        while((line = strbuff.readLine()) != null)
          resBuff.append(line);
   
        strbuff.close();
        in.close();

        res = resBuff.toString();

        if(res.indexOf("SRS error") > -1)
          return "";
 
//      System.out.println(DataCollectionPane.srs_url+
//                         "/wgetz?-f+id+[uniprot-acc:"+hit.getAcc()+"]%3E"+DB);
      }
      catch(MalformedURLException e) {System.err.println(e);}
      catch(IOException e) {System.err.println(e);}

    }
    else
    {
      String cmd3[] = { "getz", "-f", "id",
             "[libs={uniprot}-acc:"+hit.getID()+"]>"+DB };
      ExternalApplication app = new ExternalApplication(cmd3,env,null);
      res = app.getProcessStdout();
    }

    return res;
  }

  
  
  
  
  
  /**
  *
  * Link Uniprot to the another database (e.g. EMBL or ENZYME)
  *
  */
  protected static String getUniprotLinkToDatabaseByMFetch( 
                                                   final boolean isLocalMfetchExists, 
                                                   final String mfetchList,
                                                   final String env[], final String DB)
  {
    String res = null;

    if(isLocalMfetchExists)
    {
      final String cmd[]   = { "mfetch", "-f", "id",
          "-d", "uniprot", "-i", "acc:"+mfetchList, 
          "-L", DB };

      ExternalApplication app = new ExternalApplication(cmd,
                                  env,null);
      res = app.getProcessStdout();

    }
    else if(remoteMfetch)
    {
      final String cmd   = 
        "mfetch -f id -d uniprot -i \"acc:"+mfetchList+"\" -L "+DB ;
      uk.ac.sanger.artemis.j2ssh.SshPSUClient ssh =
        new uk.ac.sanger.artemis.j2ssh.SshPSUClient(cmd);
      ssh.run();
      res = ssh.getStdOut();
    }


    return res;
  }
  
  
  protected void show(Object obj)
  {
    if(obj instanceof HitInfo)
    {
      HitInfo hit = (HitInfo)obj;

      int start = hit.getStartPosition();
//    int end   = hit.getEndPosition();
//    textArea.moveCaretPosition(end);
      textArea.moveCaretPosition(start);

      Point pos  = getViewport().getViewPosition();
      Dimension rect = getViewport().getViewSize();
      double hgt = rect.getHeight()+pos.getY();
      pos.setLocation(pos.getX(),hgt);
      getViewport().setViewPosition(pos);
    }
  }

  public static boolean isRemoteMfetch()
  {
    return remoteMfetch;
  }

  public static void setRemoteMfetch(boolean remoteMfetch)
  {
    FastaTextPane.remoteMfetch = remoteMfetch;
  }

  public static boolean isForceUrl()
  {
    return forceUrl;
  }

  public static void setForceUrl(boolean forceUrl)
  {
    FastaTextPane.forceUrl = forceUrl;
  }

  public static File getMfetchExecutable()
  {
    return mfetchExecutable;
  }
 

}

