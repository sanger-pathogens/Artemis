/* SearchResultViewer.java
 *
 * created: Thu Feb 24 2000
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/SearchResultViewer.java,v 1.5 2008-01-23 14:22:49 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.ExternalProgram;
import uk.ac.sanger.artemis.ExternalProgramException;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.ZipFileDocument;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

import javax.swing.JButton;



/**
 *  A component that displays the results of external searches, with the
 *  ability to send the results to a netscape process.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: SearchResultViewer.java,v 1.5 2008-01-23 14:22:49 tjc Exp $
 **/

public class SearchResultViewer extends FileViewer 
{
  /** */
  private static final long serialVersionUID = 1L;

  /**
   *  Create a new SearchResultViewer component.
   *  @param title The name to attach to the new JFrame.
   *  @param file_name The file to read into the new viewer.
   **/
  public SearchResultViewer(final String label,
                            final Document document)
      throws IOException 
  {
    super (label, false, false, false);

    try 
    {
      readFile(document.getReader());
    } 
    catch (IOException e) 
    {
      System.out.println("error while reading results: " + e);
      dispose();
      throw e;
    }
    setVisible(true);
    if(!Options.getOptions().getPropertyTruthValue("sanger_options")) 
      return;

    final JButton to_browser = new JButton("Send to browser");
    getButtonPanel().add(to_browser);

    to_browser.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        String fileName = document.toString();
        try 
        {
          if(document instanceof ZipFileDocument)
            fileName = ((ZipFileDocument)document).writeTmpFile(
                SearchResultViewer.this.getText());
          
          sendToBrowser(fileName);
        }
        catch (IOException e) 
        {
          System.out.println ("error while reading results: " + e);
          new MessageDialog(SearchResultViewer.this,
                            "Message",
                            "Send to browser failed: " + e);
        } 
        catch(ExternalProgramException e) 
        {
          System.out.println("error while reading results: " + e);
          new MessageDialog(SearchResultViewer.this,
                            "Message",
                            "Send to browser failed: " + e);
        }
      }
    });
  }

  /**
   *  Mark up the contents of the given file (which should contain blast or
   *  fasta output) and send it to a web browser (with netscape -remote).
   **/
  protected static void sendToBrowser(final String file_name)
      throws IOException, ExternalProgramException
  {
    final String[] arguments =
    {
      file_name
    };

    final Process process =
      ExternalProgram.startProgram("results_to_netscape", arguments);

    new ProcessWatcher(process, "results_to_netscape", false);
  }
}
