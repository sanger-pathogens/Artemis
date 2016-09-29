package uk.ac.sanger.artemis.components;

/* SelectionViewer.java
 *
 * created: Thu Mar  4 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 * $Header: /data/cvs/ArtemisV8/artemis/uk/ac/sanger/artemis/components/SelectionGAMOLAViewer.java,v 1.7 2013/07/25 00:15:28 $
 **/
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.ArrayList;
import java.util.ListIterator;
import java.util.logging.Level;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Qualifier;

/**
 * This component displays the current selection in a FileViewer component.
 *
 * @author Kim Rutherford
 * @version $Id: SelectionGAMOLAViewer.java,v 1.3 2006/05/09 02:33:20
 * @author warrenj
 *
 * @version $Id: SelectionGAMOLAViewer.java,v 1.3 2013/08/19 02:33:20
 * @author Luj Exp $
 * @Author Luj 19/08/2013 
 * added this java file, which shows annotation file when feature selected. It is very similar to SelectionViewer.java
 */
public class SelectionCustomViewer implements SelectionChangeListener,
        EntryGroupChangeListener {

    /**
     * Create a new SelectionViewer object from the given Selection.
     */
//String representing the GAMOLA sub directory that contains the extra info needed eg blast_database
    String extraInfoDir;
    //String representing the dir containing the source .gb or embl file that coordinates the
    //GAMOLA Annotation
    String entryDir;

    public SelectionCustomViewer(final Selection selection,
            final EntryGroup entry_group, String extraInfoDir) {
        System.out.println("GAMOLA View should show now");
        this.extraInfoDir = extraInfoDir;
        //file_viewer = new FileViewer("Artemis GAMOLA View " + extraInfoDir);
        Entry entryJon = entry_group.getSequenceEntry();
        entryDir = entryJon.getRootDocument().toString();
        System.out.println("entryDir= " + entryDir);

        this.selection = selection;
        this.entry_group = entry_group;
        // System.out.println("Selection= "+selection +"entry_group=
        // "+entry_group);

        // readSelection ();
        readDir();

        //selection.addSelectionChangeListener(this);

        //file_viewer.addWindowListener(new WindowAdapter() {
//			public void windowClosed(WindowEvent event) {
//				stopListening();
//			}
//		});

        //entry_group.addEntryGroupChangeListener(this);
    }

    /**
     * Remove this object as a selection change listener.
     */
    public void stopListening() {
        selection.removeSelectionChangeListener(this);
        entry_group.removeEntryGroupChangeListener(this);
    }

    /**
     * Implementation of the SelectionChangeListener interface. We listen to
     * SelectionChange events so that we can update the view to reflect the
     * current selection.
     */
    public void selectionChanged(SelectionChangeEvent event) {
        readSelection();
    }

    /**
     * Implementation of the EntryGroupChangeListener interface. We listen to
     * EntryGroupChange events so that we can get rid of the JFrame when the
     * EntryGroup is no longer in use (for example when the EntryEdit is
     * closed).
     */
    public void entryGroupChanged(final EntryGroupChangeEvent event) {
        switch (event.getType()) {
            case EntryGroupChangeEvent.DONE_GONE:
                entry_group.removeEntryGroupChangeListener(this);
                //file_viewer.dispose();
                break;
        }
    }

    /**
     * Read the selection into this SelectionViewer object.
     */
    public void readSelection() {
        final String base_selection_text = SelectionInfoDisplay
                .markerRangeText(selection, entry_group);
        // System.out.println("base selection= "+base_selection_text);

        final FeatureVector selection_features = selection.getAllFeatures();

        final StringBuffer buffer = new StringBuffer();

        final int MAX_FEATURES = 50;

        if (selection_features.size() > MAX_FEATURES) {
            buffer.append("first " + MAX_FEATURES + " features:\n\n");
        }

        for (int i = 0; i < selection_features.size() && i < 50; ++i) {
            buffer.append(selection_features.elementAt(i).toString());

            // System.out.println("selection.toString=
            // "+selection_features.elementAt (i).toString ());
        }

        if (selection_features.size() > 0) {
            buffer.append("\n");
        }

        if (base_selection_text != null && base_selection_text.length() > 0) {
            buffer.append(base_selection_text).append("\n\n");
        }

        // this is the maximum number of bases to show
        final int display_base_count;

        // this is the maximum number of residues to show
        final int display_residues_count;

        if (selection_features.size() == 0) {
            display_base_count = 2000;
            display_residues_count = 1000;
        } else {
            display_base_count = 600;
            display_residues_count = 300;
        }

        final String selection_bases = selection.getSelectedBases();

        if (selection_bases != null && selection_bases.length() > 0) {
            final StringVector base_summary = getBaseSummary(selection_bases);

            for (int i = 0; i < base_summary.size(); ++i) {
                buffer.append(base_summary.elementAt(i)).append("\n");
            }
        }

        if (selection_features.size() == 1) {
            double score = selection_features.elementAt(0)
                    .get12CorrelationScore();

            score = Math.round(score * 100) / 100.0;

            buffer.append("\nCorrelation score of the selected feature: "
                    + score + "\n");
        }

        buffer.append("\n");

        if (selection_features.size() > 1) {
            double correlation_score_total = 0;

            double max_gc_content = -999;
            double min_gc_content = 999;

            double max_correlation_score = -999999;
            double min_correlation_score = 999999;

            for (int i = 0; i < selection_features.size(); ++i) {
                final Feature this_feature = selection_features.elementAt(i);

                final double correlation_score = this_feature
                        .get12CorrelationScore();
                final double gc_content = this_feature.getPercentGC();

                correlation_score_total += correlation_score;

                if (min_correlation_score > correlation_score) {
                    min_correlation_score = correlation_score;
                }

                if (max_correlation_score < correlation_score) {
                    max_correlation_score = correlation_score;
                }

                if (min_gc_content > gc_content) {
                    min_gc_content = gc_content;
                }

                if (max_gc_content < gc_content) {
                    max_gc_content = gc_content;
                }
            }

            min_gc_content = Math.round(min_gc_content * 100) / 100.0;
            max_gc_content = Math.round(max_gc_content * 100) / 100.0;
            min_correlation_score = Math.round(min_correlation_score * 100) / 100.0;
            max_correlation_score = Math.round(max_correlation_score * 100) / 100.0;

            final double correlation_score_average = Math
                    .round(correlation_score_total / selection_features.size()
                    * 100) / 100.0;

            buffer
                    .append("Average correlation score of the selected features: "
                    + correlation_score_average + "\n\n");

            buffer.append("Min GC percentage of the selected features: "
                    + min_gc_content + "\n");
            buffer.append("Max GC percentage of the selected features: "
                    + max_gc_content + "\n");
            buffer.append("\n");
            buffer.append("Min correlation score of the selected features: "
                    + min_correlation_score + "\n");
            buffer.append("Max correlation score of the selected features: "
                    + max_correlation_score + "\n");
            buffer.append("\n");
        }

        if (selection_bases.length() > 0 && selection_features.size() <= 1) {
            if (selection_bases.length() > display_base_count) {
                // display the first display_base_count/2 bases and the last
                // display_base_count/2 bases
                buffer.append("first " + display_base_count / 2
                        + " bases of selection:\n");
                final String start_sub_string = selection_bases.substring(0,
                        display_base_count / 2);
                buffer.append(start_sub_string).append("\n\n");
                buffer.append("last " + display_base_count / 2
                        + " bases of selection:\n");
                final String sub_string = selection_bases
                        .substring(selection_bases.length()
                        - display_base_count / 2);
                buffer.append(sub_string).append("\n\n");
            } else {
                buffer.append("bases of selection:\n");
                buffer.append(selection_bases).append("\n\n");
            }

            if (selection_bases.length() >= 3) {
                final int residues_to_display;

                if (selection_bases.length() / 3 < display_residues_count) {
                    residues_to_display = selection_bases.length() / 3;

                    buffer.append("translation of the selected bases:\n");
                } else {
                    residues_to_display = display_residues_count;

                    buffer
                            .append("translation of the first "
                            + display_residues_count * 3
                            + " selected bases:\n");
                }

                final String residue_bases = selection_bases.substring(0,
                        residues_to_display * 3);

                final String residues = AminoAcidSequence.getTranslation(
                        residue_bases, true).toString();

                buffer.append(residues.toUpperCase()).append("\n\n");
            }
        }

        //file_viewer.setText(buffer.toString());
    }

    // jw written to replace the readSelection method when we want to search the
    // GAMOLA directories for
    // extra info
    public void readDir() {
        final String base_selection_text = SelectionInfoDisplay
                .markerRangeText(selection, entry_group);
        final FeatureVector selection_features = selection.getAllFeatures();


        final int MAX_FEATURES = 50;

        if (selection_features.size() == 0) {
            String warningMsg = "No feature selected - please select a feature for which you would like to view custom annotation.";
            JOptionPane.showMessageDialog(null, warningMsg, "Warning", JOptionPane.WARNING_MESSAGE);
        } else 
        {
            for (int i = 0; i < selection_features.size() && i < 50; ++i) {
                //buffer.append(selection_features.elementAt(i).toString());
                Feature feature = selection_features.elementAt(i);
                String geneName = feature.getGeneName();
                System.out.println("GAMOLA key= " + geneName + "key= "
                        + feature.getKey());

                Qualifier q;
                String evidence_id = "";
                try {
                    q = feature.getQualifierByName("evidence_id");
                    if (q != null) {
                        StringVector sv = q.getValues();
                        if (sv != null) {
                            evidence_id = sv.elementAt(0);
                        }
                    }
                } catch (InvalidRelationException ex) {
                    java.util.logging.Logger.getLogger(SelectionCustomViewer.class.getName()).log(Level.SEVERE, null, ex);
                }


                // need to get number from the end of the feature key (whatever I
                // use as the key)
                // then need to find directory (up and then COG or whatever option
                // has been selected
                // can get this from type attribute e.g. gene or COG etc?

                /*
                 if (evidence_id != "") {
                 //if (geneName != null) {

                 ArrayList files = getExtraInfoFile(evidence_id,
                 extraInfoDir);

                 ListIterator it = files.listIterator();
                 while (it.hasNext()) {
                 this.openInFileViewer((File) it.next());
                 //FileViewer.readFile((File) it.next());
                 }
                 }*/

                if (geneName != null) {

                    ArrayList files = getExtraInfoFile(geneName,
                            extraInfoDir);

                    ListIterator it = files.listIterator();
                    while (it.hasNext()) {
                        this.openInFileViewer((File) it.next(), geneName);
                    }
                }
            }
        }
    }

    private StringBuffer readFile(File extraGAMOLAInformationFile) {
        StringBuffer toView = new StringBuffer("");
        try {
            BufferedReader in = new BufferedReader(new FileReader(
                    extraGAMOLAInformationFile));
            String str;
            while ((str = in.readLine()) != null) {
                toView.append(str + "\n");
                //System.out.println(str);
            }
            in.close();
        } catch (IOException e) {
        }
        return toView;
    }

    /**
     * get the file that contains the extra info needed eg pfam or COG etc
     *
     * @param geneName name that contains the ORF number being refered to.
     * @param extraInfoDir - dir name needed corresponding to type of info eg
     * COG etc
     */
    /*
     private ArrayList getExtraInfoFile(String evidenceID, String extraInfoDir) {
     //Arraylist to hold files of interest
     ArrayList extraInfoFiles = new ArrayList();
     System.out.println("evidence_id= " + evidenceID + " extrainfoDir= "
     + extraInfoDir);
     File locationOfGenbankFile = new File(entryDir);

     File dir = locationOfGenbankFile;

     if (dir.isDirectory()) {
     File[] infoDirectories = dir.listFiles();
     if (infoDirectories == null) {
     // Either dir does not exist or is not a directory
     } else {
     for (int j = 0; j < infoDirectories.length; j++) {
     // Get filename of file or directory
     File infoDir = infoDirectories[j];
     String dirName = infoDir.getName();
     System.out.println(dirName);
     if (dirName.equalsIgnoreCase(extraInfoDir)
     && infoDir.isDirectory()) {
     // infoDir is the pfam one
     // so list the files under this dir and find the
     // file that ends with the
     // appropriate number signifying the ORF number.
     File[] infoFiles = infoDir.listFiles();
     for (int h = 0; h < infoFiles.length; h++) {
     File infoFile = infoFiles[h];
     //System.out.println("files in extraDir= "
     //+ infoFile.getName());
     if (infoFile.getName().substring(
     infoFile.getName().lastIndexOf("_") + 1,
     infoFile.getName().length())
     .equalsIgnoreCase(evidenceID)) {
     System.out.println("extraInfo File found "
     + infoFile.getName());
     * 
     extraInfoFiles.add(infoFile);
     }
     }

     }
     }
     }

     }
     return extraInfoFiles;
     }
    
     */
    private ArrayList getExtraInfoFile(String geneName, String extraInfoDir) {
        //Arraylist to hold files of interest
        ArrayList extraInfoFiles = new ArrayList();
        String afterDash = geneName.substring(geneName.lastIndexOf("_") + 1,
                geneName.length());
        System.out.println("afterDash= " + afterDash + " extrainfoDir= "
                + extraInfoDir);
        File locationOfGenbankFile = new File(entryDir);

        File dir = locationOfGenbankFile;

        if (dir.isDirectory()) {
            File[] infoDirectories = dir.listFiles();
            if (infoDirectories == null) {
                // Either dir does not exist or is not a directory
            } else {
                for (int j = 0; j < infoDirectories.length; j++) {
                    // Get filename of file or directory
                    File infoDir = infoDirectories[j];
                    String dirName = infoDir.getName();
                    System.out.println(dirName);
                    if (dirName.equalsIgnoreCase(extraInfoDir)
                            && infoDir.isDirectory()) {
                        // infoDir is the pfam one
                        // so list the files under this dir and find the
                        // file that ends with the
                        // appropriate number signifying the ORF number.
                        File[] infoFiles = infoDir.listFiles();
                        for (int h = 0; h < infoFiles.length; h++) {
                            File infoFile = infoFiles[h];
                            //System.out.println("files in extraDir= "
                            //+ infoFile.getName());
                            if (infoFile.getName().substring(
                                    infoFile.getName().lastIndexOf("_") + 1,
                                    infoFile.getName().length())
                                    .equalsIgnoreCase(afterDash)) {
                                System.out.println("extraInfo File found "
                                        + infoFile.getName());
                                extraInfoFiles.add(infoFile);
                            }
                        }

                    }
                }
            }

        }
        return extraInfoFiles;
    }

    /**
     * Return a summary of the base content of the given String. The return
     * value is a StringVector containing each line to be output. We return a
     * vector so that the caller can do things like indenting the lines.
     *
     * @return null if the bases String is empty or null
     */
    public static StringVector getBaseSummary(final String bases) {
        if (bases == null || bases.length() == 0) {
            return null;
        }

        final StringVector return_vector = new StringVector();

        long counts[] = new long[257];

        for (int i = 0; i < bases.length(); ++i) {
            final char this_char = bases.charAt(i);

            if (this_char > 255) {
                // don't know what this is, but it isn't a base
                ++counts[256];
            } else {
                ++counts[this_char];
            }
        }

        for (int i = 0; i < 256; ++i) {
            if (counts[i] > 0) {
                return_vector.add(new String(Character.toUpperCase((char) i)
                        + " content: " + counts[i] + "  (" + 10000 * counts[i]
                        / bases.length() / 100.0 + "%)"));
            }
        }

        final long non_ambiguous_base_count = counts['g'] + counts['c']
                + counts['a'] + counts['t'];

        if (non_ambiguous_base_count == bases.length()) {
            return_vector.add("(no ambiguous bases)");
        }

        return_vector.add("");

        return_vector.add(new String("GC percentage: "
                + (Math.round(100.0 * (counts['g'] + counts['c'])
                / bases.length() * 100) / 100.0)));

        if (non_ambiguous_base_count > 0
                && non_ambiguous_base_count < bases.length()) {
            return_vector
                    .add(new String("GC percentage of non-ambiguous bases: "
                    + (Math.round(100.0 * (counts['g'] + counts['c'])
                    / non_ambiguous_base_count * 100) / 100.0)));
        }

        return return_vector;
    }

    protected void openInFileViewer(File fileToOpen, String geneName) {
        BufferedReader br = null;
        try {
            String fileName = fileToOpen.getAbsolutePath();
            FileViewer file_viewer = new FileViewer(fileToOpen.getName());
            
            //Add gene name to the view custom annotation result. 
            file_viewer.appendString("gene=" + geneName + "\n"); 
            
            br = new BufferedReader(new FileReader(fileName));
            String line;
            while ((line = br.readLine()) != null) {
                file_viewer.appendString(line + "\n");
            }
            br.close();

        } catch (IOException ex) {
            java.util.logging.Logger.getLogger(SelectionCustomViewer.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                java.util.logging.Logger.getLogger(SelectionCustomViewer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    protected void openInWordPad(File fileToOpen) {
        String[] cmd = {"C:/Program Files/Windows NT/Accessories/wordpad.exe", fileToOpen.getAbsolutePath()};

        try {
            Process p = Runtime.getRuntime().exec(cmd);
            p.exitValue();
            //p.waitFor();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IllegalThreadStateException e) {
            //don't hold java up and don't worry about the process not finishing and
            //just carry on and let wordpad do what it likes!
            //e.printStackTrace();
        }
    }
    /**
     * This is the selection object that we are viewing.
     */
    final private Selection selection;
    /**
     * The FileViewer object that is displaying the selection.
     */
    //final private FileViewer file_viewer;
    /**
     * The EntryGroup that was passed to the constructor.
     */
    final EntryGroup entry_group;
}
