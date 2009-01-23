/*
 * Copyright (C) 2009  Genome Research Limited
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
 *
 */
package uk.ac.sanger.artemis.circular.digest;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.circular.Block;
import uk.ac.sanger.artemis.circular.DNADraw;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.OutOfRangeException;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;

import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.event.MouseMotionListener;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

/**
 * 
 */
public class CircularGenomeController
{

	public CircularGenomeController()
	{
	}

	/**
	 * Create in-silico PFGE from a restriction enzyme digest and draw
	 * alongside DNAPlotter
	 */
	protected void setup(CircularGenomeCommandBean cgcb, String fileName) throws Exception
	{
		try
		{
			EntryGroup entryGroup = getEntryGroupFromFile(fileName);		
			if(fileName == null)
			{
				fileName = 
		   		((File)entryGroup.getSequenceEntry().getRootDocument().getLocation()).getAbsolutePath()+
				  File.separator+entryGroup.getSequenceEntry().getName();
			}
			
			if(cgcb.getEnzymeName() == null)
			{
				// prompt for enzymes
				promptForEnzymes(cgcb);
			}
			
			File output = File.createTempFile("circular_genome", ".txt");
			String[] args =
			{
				System.getProperty("EMBOSS_ROOT") + "/bin/restrict", 
				fileName, "-auto", 
				"-limit", "y", 
				"-enzymes", cgcb.getEnzymeName(), 
				"-out", output.getCanonicalPath() 
			};
			
			ProcessBuilder pb = new ProcessBuilder(args);
			pb.redirectErrorStream(true);
			Process p = pb.start();
			System.err.print("**");
			try
			{
				InputStream is = p.getInputStream();
				int inchar;
				while ((inchar = is.read()) != -1)
				{
					char c = (char) inchar;
					System.err.print(c);
				}
				System.err.println("**");
				p.waitFor();
				System.err.println("Process exited with '" + p.exitValue() + "'");
			} 
			catch (InterruptedException exp)
			{
				exp.printStackTrace();
			}
			p.waitFor();
			
			final ReportDetails rd = Utils.findCutSitesFromEmbossReport(output.getCanonicalPath());

			if (rd.cutSites.size() == 0)
				JOptionPane.showMessageDialog(null,
						 "No cut site found.",
						 "RE Digest Results", JOptionPane.INFORMATION_MESSAGE);
		 
			final DNADraw dna = Utils.createDNADrawFromReportDetails(rd, entryGroup);

	    final InSilicoGelPanel inSilicoGelPanel = 
	    	new InSilicoGelPanel(rd.length, rd.cutSites, dna.getHeight(), output);
	    
			MouseMotionListener mouseMotionListener = new MouseMotionAdapter() 
	    {
	      public void mouseMoved(MouseEvent me)
	      {
	      	List<CutSite> cutSites = rd.cutSites;
	        final Block b = dna.getBlockAtLocation(me.getPoint());
	        if(b != null && b.isOverMe(me.getX(), me.getY()))
	        {
	          int bend   = b.getBend();
	          int bstart = b.getBstart();

	          if(bend == rd.length)
	          {
	          	((CutSite)cutSites.get(0)).setHighlighted(true);
	          	for (int i = 1; i < cutSites.size(); i++)
	          		((CutSite) cutSites.get(i)).setHighlighted(false);
	          }
	          else
						{
							for (int i = 0; i < cutSites.size(); i++)
							{
								CutSite cutSite = (CutSite) cutSites.get(i);

								if (bend == cutSite.getFivePrime())
									cutSite.setHighlighted(true);
								else
									cutSite.setHighlighted(false);
							}
						}
	        }
	        else
	        {
	        	for(int i=0; i<cutSites.size(); i++)
		        {
		          CutSite cutSite = (CutSite) cutSites.get(i);
		          cutSite.setHighlighted(false);
		        }
	        }
	        inSilicoGelPanel.repaint();
	      }
	    };
	    dna.addMouseMotionListener(mouseMotionListener);
			
			final JFrame f = new JFrame(cgcb.getEnzymeName()+" ::: "+fileName);
	    Dimension d = f.getToolkit().getScreenSize();

	    JPanel mainPanel = new JPanel(new FlowLayout(FlowLayout.CENTER, 0, 0));
	    mainPanel.add(inSilicoGelPanel);
	    mainPanel.add(dna);
	    
	    JScrollPane jsp = new JScrollPane(mainPanel);
	    jsp.getViewport().setBackground(Color.white);
	    f.getContentPane().add(jsp);
	    f.setJMenuBar(dna.createMenuBar());

	    f.pack();
	    f.setLocation(((int)d.getWidth()-f.getWidth())/4,
	                  ((int)d.getHeight()-f.getHeight())/2);
	    f.setVisible(true);
		} 
		catch (IOException exp)
		{
			exp.printStackTrace();
		}
	}
	
  /**
   * Create a DNADraw panel from a file
   * @param dna_current
   * @return
   */
  private static EntryGroup getEntryGroupFromFile(String fileName)
  {
    Options.getOptions();
    final EntryGroup entryGroup = new SimpleEntryGroup();
   
    try
    {
    	Entry entry = getEntry(fileName);
      entryGroup.add(entry);
      return entryGroup;
    }
    catch(OutOfRangeException e)
    {
      e.printStackTrace();
    }
    catch(NoSequenceException e)
    {
      JOptionPane.showMessageDialog(null, "No sequence found!", 
          "Sequence Missing", JOptionPane.WARNING_MESSAGE);
    }
    return null;
  }
  
  /**
   * Return an Artemis entry from a file 
   * @param entryFileName
   * @param entryGroup
   * @return
   * @throws NoSequenceException
   * @throws OutOfRangeException 
   */
  private static Entry getEntry(final String entryFileName) 
                   throws NoSequenceException, OutOfRangeException
  {

  	if(entryFileName == null)
  	{
    	// no file - prompt for a file
  		uk.ac.sanger.artemis.components.FileDialogEntrySource entrySource = 
        new uk.ac.sanger.artemis.components.FileDialogEntrySource(null, null);
      Entry entry = entrySource.getEntry(true);
      return entry;
  	}
  	
    final Document entry_document = DocumentFactory.makeDocument(entryFileName);
    final EntryInformation artemis_entry_information =
      Options.getArtemisEntryInformation();
    
    final uk.ac.sanger.artemis.io.Entry new_embl_entry =
      EntryFileDialog.getEntryFromFile(null, entry_document,
                                       artemis_entry_information,
                                       false);

    if(new_embl_entry == null)  // the read failed
      return null;

    Entry entry = null;
    try
    {
      entry = new Entry(new_embl_entry);
    } 
    catch(OutOfRangeException e) 
    {
      new MessageDialog(null, "read failed: one of the features in " +
          entryFileName + " has an out of range " +
                        "location: " + e.getMessage());
    }
    return entry;
  }
  
  private void promptForEnzymes(CircularGenomeCommandBean cgcb)
  {
  	String enzymes = JOptionPane.showInputDialog(null, 
  			"Enzymes", "HincII,hinfI,ppiI,hindiii");
  	cgcb.setEnzymeName(enzymes);
  }

	public static void main(String args[])
	{
		if(System.getProperty("EMBOSS_ROOT") == null)
		{
			String embossRoot = JOptionPane.showInputDialog(null,
					"Input the EMBOSS installation directory",
					"/usr/local/emboss");
			System.setProperty("EMBOSS_ROOT",embossRoot.trim());
		}
		
		CircularGenomeCommandBean command   = new CircularGenomeCommandBean();
		CircularGenomeController controller = new CircularGenomeController();
		try
		{
			String fileName = null;
			if(args != null && args.length > 0)
			{
				if(args.length == 1)
		  		fileName = args[0];
				
				for(int i=0; i<args.length; i++)
				{
					if(args[i].startsWith("-enz"))
						command.setEnzymeName(args[i+1]);
					else if(args[i].startsWith("-seq"))
						fileName = args[i+1];
				}
			}
			controller.setup(command, fileName);
		} 
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

};

class CircularGenomeCommandBean 
{
	private String enzymeName;

	public String getEnzymeName()
	{
		return this.enzymeName;
	}

	public void setEnzymeName(String enzymeName)
	{
		this.enzymeName = enzymeName;
	}
}
