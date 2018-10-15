/* OrthoParalogTable.java
 * This file is part of Artemis
 *
 * Copyright (C) 2007  Genome Research Limited
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

package uk.ac.sanger.artemis.components.genebuilder.ortholog;

import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumn;
import javax.swing.table.TableModel;

import org.gmod.schema.sequence.FeatureLoc;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.chado.ArtemisUtils;
import uk.ac.sanger.artemis.chado.ClusterLazyQualifierValue;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.PartialSequence;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierLazyLoading;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

public class OrthoParalogTable extends AbstractMatchTable
{
  private static int NUMBER_COLUMNS = 10;
  private Vector rowData   = new Vector();
  private Vector tableData = new Vector(NUMBER_COLUMNS);
  private JTable table;
  private JButton infoLevelButton = new JButton("Details");
  private JPopupMenu popupMenu = new JPopupMenu();
  private boolean showCluster;
  
  //
  // column headings
  protected final static String CLUSTER_NAME_COL = "Cluster";
  protected final static String MATCH_NAME_COL = "Match";
  protected final static String ROW_TYPE_HIDE_COL = "Term";
  protected final static String ROW_TYPE_COL = "Type";
  protected final static String PROGRAM_COL = "Program";
  protected final static String ORGANISM_COL = "Organism";
  protected final static String PRODUCT_COL = "Product";
  protected final static String GENE_COL = "Gene";
  protected final static String LINK_COL = "Link";
  protected final static String VIEW_BUTTON_COL = "View";
  protected final static String REMOVE_BUTTON_COL = "";
  

  /**
   * Contruct a component for an ortholog or paralog line
   * @param doc
   * @param orthologQualifier
   * @param paralogQualifier
   * @param feature
   * @param showCluster
   */
  protected OrthoParalogTable(final DatabaseDocument doc,
                              final Qualifier orthologQualifier,
                              final Qualifier paralogQualifier,
                              final Feature feature,
                              final boolean showCluster)
  {
    this.origQualifiers = new QualifierVector();
    this.showCluster = showCluster;
    
    if(orthologQualifier != null)
      this.origQualifiers.add(orthologQualifier);
    if(paralogQualifier != null)
      this.origQualifiers.add(paralogQualifier);
    
    createPopupMenu(doc, feature);
    
    infoLevelButton.setOpaque(false);
    tableData.setSize(NUMBER_COLUMNS);
    
    tableData.setElementAt(CLUSTER_NAME_COL,0);
    tableData.setElementAt(MATCH_NAME_COL,1);
    tableData.setElementAt(ROW_TYPE_HIDE_COL,2);
    if(showCluster)
      tableData.setElementAt(PROGRAM_COL,3);
    else
      tableData.setElementAt(ROW_TYPE_COL,3);
    tableData.setElementAt(ORGANISM_COL,4);
    tableData.setElementAt(GENE_COL,5);
    tableData.setElementAt(LINK_COL,6);
    tableData.setElementAt(PRODUCT_COL,7);
    tableData.setElementAt(VIEW_BUTTON_COL,8);
    tableData.setElementAt(REMOVE_BUTTON_COL,9);
    
    // add row data
    int columnIndex;
  
    for(int i=0; i<origQualifiers.size(); i++)
    {
      final Vector qualifierValuesToDelete = new Vector();
      
      final Qualifier origQualifier = (Qualifier) origQualifiers.elementAt(i);
      
      // ensure the gene name is loaded as well
      if(origQualifier instanceof QualifierLazyLoading)
      {
        List lazyValues = ((QualifierLazyLoading)origQualifier).getLazyValues();
        for(int j=0; j<lazyValues.size(); j++)
        {
          ClusterLazyQualifierValue lazyValue = (ClusterLazyQualifierValue)lazyValues.get(j);
          
          if(!lazyValue.isLazyLoaded())
            lazyValue.setLoadGeneName(true);
        }
      }
      
      final StringVector values = origQualifier.getValues();

      // sort by their rank value
      Collections.sort(values, new OrthoParalogValueComparator());

      for(int j = 0; j < values.size(); j++)
      {
        StringVector rowStr = StringVector.getStrings((String) values.get(j),
            ";");

        if(rowStr.size() < 1)
        {
          qualifierValuesToDelete.add(values.get(j));
          continue;
        }
        if( (ArtemisUtils.getString(rowStr, "cluster_name=").equals("") && !showCluster) ||
            (!ArtemisUtils.getString(rowStr, "cluster_name=").equals("") && showCluster) )
          processRowData(rowStr, rowData, origQualifier.getName());
      }

      values.removeAll(qualifierValuesToDelete);
    }
    
    table = new JTable(rowData, tableData);
    setTable(table);
    
    // set hand cursor
    table.addMouseMotionListener( new MouseMotionAdapter() 
    {
      private Cursor handCursor = Cursor.getPredefinedCursor(Cursor.HAND_CURSOR);
      public void mouseMoved(MouseEvent e) 
      {
        int col = table.columnAtPoint(e.getPoint());
        final String colName = table.getColumnName(col);
     
        if(colName.equals(GENE_COL) || colName.equals(REMOVE_BUTTON_COL)) 
          table.setCursor(handCursor);
        else 
          table.setCursor(Cursor.getDefaultCursor());  
      }
    });
    
    table.addMouseListener(new MouseAdapter() 
    {
      public void mousePressed(MouseEvent e) 
      {
        showPopup(e);
      }

      public void mouseReleased(MouseEvent e) 
      {
        showPopup(e);
      }

      private void showPopup(MouseEvent e)
      {
        if(e.isPopupTrigger()) 
          popupMenu.show(e.getComponent(), e.getX(), e.getY());
      }
    });
    
    final TableColumn[] hideColumns = new TableColumn[4];
    
    hideColumns[0] = table.getColumn(ROW_TYPE_HIDE_COL);
    hideColumns[1] = table.getColumn(MATCH_NAME_COL);
    if(showCluster)
      hideColumns[2] = table.getColumn(REMOVE_BUTTON_COL);
    else
      hideColumns[2] = table.getColumn(CLUSTER_NAME_COL);
    hideColumns[3] = table.getColumn(PRODUCT_COL);

    for(int i=0; i<hideColumns.length; i++)
    {
      if(i == 3 && !showCluster)
        continue;
      hideColumns[i].setMinWidth(0);
      hideColumns[i].setMaxWidth(0);
    }

    table.setColumnSelectionAllowed(false);
    table.setRowSelectionAllowed(true);
    table.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
    table.setDragEnabled(true);
    table.setTransferHandler(new TableTransferHandler());
    
    TableModel tableModel = table.getModel();
    // remove button column
    TableColumn col;
    if(!showCluster)
    {
      col = table.getColumn(REMOVE_BUTTON_COL);
      col.setMinWidth(35);
      col.setMaxWidth(40);
      col.setPreferredWidth(40);
    }

    final OrthologRenderer renderer = new OrthologRenderer();

    for(columnIndex = 0; columnIndex <tableModel.getColumnCount();
        columnIndex++) 
    {
      col = table.getColumnModel().getColumn(columnIndex);
      col.setCellRenderer(renderer);
      col.setCellEditor(new CellEditing(new JTextField()));
    }

    packColumn(table, getColumnIndex(GENE_COL), 4);
    packColumn(table, getColumnIndex(LINK_COL), 4);
    packColumn(table, getColumnIndex(ORGANISM_COL), 4);
    packColumn(table, getColumnIndex(VIEW_BUTTON_COL), 4);

    if(showCluster)
    {
      packColumn(table, getColumnIndex(CLUSTER_NAME_COL), 4);
      packColumn(table, getColumnIndex(PROGRAM_COL), 4);
    }
    else
    {
      packColumn(table, getColumnIndex(ROW_TYPE_COL), 4);
      packColumn(table, getColumnIndex(PRODUCT_COL), 4);
    }
    
    // remove JButton column
    col = table.getColumn(REMOVE_BUTTON_COL);
    col.setCellEditor(new ButtonEditor(new JCheckBox(),
        (DefaultTableModel)table.getModel(), "X", doc));
    
    // remove JButton column
    col = table.getColumn(VIEW_BUTTON_COL);
    col.setCellEditor(new ButtonEditor(new JCheckBox(),
        (DefaultTableModel)table.getModel(), "VIEW", doc));
    
    // orthologue link
    col = table.getColumn(GENE_COL);
    col.setCellEditor(new LinkEditor(new JCheckBox(),
        (DefaultTableModel)table.getModel(), doc));
  }
  
  /**
   * Parse the qualifier string into row data
   * @param rowStr
   * @param rowData
   * @param qualifierName
   */
  private void processRowData(final StringVector rowStr,
                              final Vector rowData,
                              final String qualifierName)
  {
    final String orthoparalogs[] = ((String) rowStr.get(0)).split(",");

    String clusterName = "";
    if(rowStr.size() > 1)
    {
      clusterName = ArtemisUtils.getString(rowStr, "cluster_name=");
      if(!clusterName.equals(""))
        clusterName = clusterName.substring(13);
    }
    
    String matchName = "";
    if(rowStr.size() > 1)
    {
      matchName = ArtemisUtils.getString(rowStr, "match_name=");
      if(!matchName.equals(""))
        matchName = matchName.substring(11);
    }
    
    String program = "";
    if(rowStr.size() > 1)
    {
      program = ArtemisUtils.getString(rowStr, "program=");
      if(!program.equals(""))
        program = program.substring(8);
    }
    
    String product = "";
    if(rowStr.size() > 1)
    {
      product = ArtemisUtils.getString(rowStr, "product=");
      if(!product.equals(""))
        product = product.substring(8);
    }
    
    int columnIndex;
    for(int k = 0; k < orthoparalogs.length; k++)
    {
      Vector thisRowData = new Vector(NUMBER_COLUMNS);
      
      String geneNameAndLinkAndType[] = orthoparalogs[k].split("link=");
      String linkAndType[] = geneNameAndLinkAndType[1].split("type=");
      String gene[] = geneNameAndLinkAndType[0].trim().split(":");
      
      thisRowData.setSize(NUMBER_COLUMNS);
      
      columnIndex = tableData.indexOf(ORGANISM_COL);
      thisRowData.setElementAt(gene[0], columnIndex);
      
      columnIndex = tableData.indexOf(GENE_COL);
      thisRowData.setElementAt(geneNameAndLinkAndType[0].trim(), columnIndex);
      
      columnIndex = tableData.indexOf(LINK_COL);
      thisRowData.setElementAt(linkAndType[0].trim(), columnIndex);
      
      columnIndex = tableData.indexOf(CLUSTER_NAME_COL);
      thisRowData.setElementAt(clusterName, columnIndex);

      columnIndex = tableData.indexOf(MATCH_NAME_COL);
      thisRowData.setElementAt(matchName, columnIndex);
      
      columnIndex = tableData.indexOf(ROW_TYPE_HIDE_COL);
      thisRowData.setElementAt(qualifierName, columnIndex);
      
      columnIndex = tableData.indexOf(PRODUCT_COL);
      thisRowData.setElementAt(product, columnIndex);
      
      columnIndex = tableData.indexOf(ROW_TYPE_COL);
      
      if(columnIndex > -1)
      {
        final String symbol;
        if(linkAndType[1].trim().equals(MatchPanel.ORTHOLOG))
          symbol = "O";
        else
         symbol = "P";
        thisRowData.setElementAt(symbol, columnIndex);
      }
      
      columnIndex = tableData.indexOf(PROGRAM_COL);
      if(columnIndex > -1)
        thisRowData.setElementAt(program, columnIndex);
      
      rowData.add(thisRowData);
    }  
  }
  
  /**
   * Find if this contains any clusters
   * @return
   */
  protected static boolean hasCluster(final Qualifier orthoQualifier,
                                      final Qualifier paraQualifier,
                                      final GFFStreamFeature feature)
  {
    if(hasClusterOrOrthoParalog(true, orthoQualifier, feature))
      return true;
    return hasClusterOrOrthoParalog(true, paraQualifier, feature);
  }
  
  /**
   * Find if this contains any ortholog or paralog
   * @return
   */
  protected static boolean hasOrthoParlaog(final Qualifier orthoQualifier,
                                           final Qualifier paraQualifier,
                                           final GFFStreamFeature feature)
  {
    if(hasClusterOrOrthoParalog(false, orthoQualifier, feature))
    {
      forceBulkLoad(paraQualifier, feature);
      return true;
    }
    return hasClusterOrOrthoParalog(false, paraQualifier, feature);
  }
  
  private static boolean hasClusterOrOrthoParalog(final boolean lookForCluster, 
                                                  final Qualifier qualifier,
                                                  final GFFStreamFeature feature)
  {
    if(qualifier == null)
      return false;
    
    forceBulkLoad(qualifier, feature);
    
    StringVector values = qualifier.getValues();

    for(int j = 0; j < values.size(); j++)
    {
      StringVector rowStr = StringVector.getStrings((String) values.get(j),
          ";");

      if( (!ArtemisUtils.getString(rowStr, "cluster_name=").equals("") && lookForCluster) ||
          (ArtemisUtils.getString(rowStr, "cluster_name=").equals("")  && !lookForCluster) )
        return true;
    }
 
    return false;
  }
  
  /**
   * For long lists of ortho/paralogs this speeds-up their loading
   * @param qualifier
   * @param feature
   */
  private static void forceBulkLoad(final Qualifier qualifier,
                             final GFFStreamFeature feature)
  {
    if(qualifier instanceof QualifierLazyLoading && 
        !((QualifierLazyLoading)qualifier).isAllLazyValuesLoaded())
    {
      List values = ((QualifierLazyLoading)qualifier).getLazyValues();
      final DatabaseDocument document =
        (DatabaseDocument)((DocumentEntry)feature.getEntry()).getDocument();
      ClusterLazyQualifierValue.setClusterFromValueList(values, document);
    }    
  }
  
  /**
   * Find the parent feature
   * @param feature
   * @return
   */
  private Feature getParentFeature(final Feature feature)
  {
    final QualifierVector qualifiers = feature.getQualifiers();
    
 
    String featureId = null;
    for(int i=0; i<qualifiers.size(); i++)
    {
      Qualifier qualifier = (Qualifier)qualifiers.elementAt(i);
      if(qualifier.getName().equalsIgnoreCase("Derives_from") ||
         qualifier.getName().equalsIgnoreCase("Parent"))
      {
        featureId = (String) qualifier.getValues().get(0);
        break;
      }
    }
    FeatureNamePredicate predicate = new FeatureNamePredicate(featureId);
    final uk.ac.sanger.artemis.FeatureVector features = feature.getEntry().getAllFeatures();
    for(int i=0; i<features.size(); i++)
    {
      Feature thisFeature = features.elementAt(i);
      if(predicate.testPredicate(thisFeature))
        return thisFeature;
    }
    return null;
  }
  
  
  /**
   * Create the popup menu for the table
   *
   */
  private void createPopupMenu(final DatabaseDocument doc,
                               final Feature feature)
  {
    JMenuItem showSequenceMenu = new JMenuItem("Show selected sequences");
    popupMenu.add(showSequenceMenu);
    showSequenceMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Cursor orginalCursor = table.getCursor();
        table.setCursor(new Cursor(Cursor.WAIT_CURSOR)); 
        showAlignmentEditor(feature, doc, false);
        table.setCursor(orginalCursor); 
      }
    });


    JMenuItem showAASequenceMenu = new JMenuItem("Show selected amino acid sequences");
    popupMenu.add(showAASequenceMenu);
    showAASequenceMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Cursor orginalCursor = table.getCursor();
        table.setCursor(new Cursor(Cursor.WAIT_CURSOR)); 
        showAlignmentEditor(feature, doc, true);
        table.setCursor(orginalCursor); 
      }
    });
    
    JMenuItem openMenu = new JMenuItem("Open selected in Artemis");
    popupMenu.add(openMenu);
    openMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Cursor orginalCursor = table.getCursor();
        table.setCursor(new Cursor(Cursor.WAIT_CURSOR));
        int selectedRows[] = getTable().getSelectedRows();
        
        if(selectedRows.length > 1)
        {
          int select = JOptionPane.showConfirmDialog(null, 
              "Open all selected sequences in seperate Artemis windows?", 
              "Open Artemis x"+selectedRows.length, JOptionPane.OK_CANCEL_OPTION);
          if(select == JOptionPane.CANCEL_OPTION)
            return;
        }
        for(int i=0; i<selectedRows.length; i++)
          openArtemis(doc,selectedRows[i]);
        table.setCursor(orginalCursor); 
      }
    });
  }
  
  /**
   * Display this feature and selected features in the table in 
   * Jemboss alignment editor.
   * @param feature
   * @param doc
   * @param showPeptideSequence
   */
  private void showAlignmentEditor(final Feature feature, 
                                   final DatabaseDocument doc, 
                                   final boolean showPeptideSequence)
  {
    final int[] rows = table.getSelectedRows();
    final int orthoColumn = getColumnIndex(GENE_COL);
    final Vector seqs = new Vector();

    // find the gene feature
    Feature gene = feature;
    if(!feature.getKey().equals("gene"))
    {
      gene = getParentFeature(feature);
      if(!gene.getKey().equals("gene"))
        gene = getParentFeature(gene);
    }
    
    // find the exons for the gene
    if(gene != null)
    {
      ChadoCanonicalGene chadoGene = ((GFFStreamFeature)gene.getEmblFeature()).getChadoGene();
      StringBuffer buffer = new StringBuffer();
      try
      {
        String transcriptName =
          chadoGene.getTranscriptFromName((String) feature.getQualifierByName("ID").getValues().get(0));
        
        List exons = chadoGene.getSplicedFeaturesOfTranscript(transcriptName);
        
        for (int i = 0 ; i < exons.size () ; ++i) 
        {
          final uk.ac.sanger.artemis.io.Feature this_feature =
            (uk.ac.sanger.artemis.io.Feature)exons.get(i);
          buffer.append (((Feature)this_feature.getUserData()).getBases ());
        }
        
        final String seqStr;
        if(showPeptideSequence)
          seqStr = AminoAcidSequence.getTranslation (buffer.toString(), true).toString();
        else
          seqStr = buffer.toString();
        
        final String sysName = gene.getSystematicName();
        seqs.add(new org.emboss.jemboss.editor.Sequence(sysName, seqStr));
      }
      catch(InvalidRelationException e)
      {
        e.printStackTrace();
      }
    }
    
    for(int i=0; i<rows.length; i++)
    {
      String ortho = (String)table.getValueAt(rows[i], orthoColumn);
      final String reference[] = ortho.split(":");
      DatabaseDocument newdoc = new DatabaseDocument(doc, 
          reference[0], reference[1], true, stream_progress_listener);
      
      try
      {
        // gene sequence
        PartialSequence sequence = newdoc.getChadoSequence(reference[1]);
        
        // cds featureloc's
        List featureLocs = 
          newdoc.getCdsFeatureLocsByPeptideName(
              (String)table.getValueAt(rows[i], getColumnIndex(LINK_COL)));
        
        //
        int phase = 0;
        final StringBuffer sequenceBuffer = new StringBuffer();
        if(featureLocs != null)
        {
          for(int j = 0; j < featureLocs.size(); j++)
          {
            FeatureLoc featureLoc = (FeatureLoc) featureLocs.get(j);
            if(featureLoc.getPhase() != null)
              phase = featureLoc.getPhase().intValue();
            char[] subSeq = sequence.getCharSubSequence(
                (featureLoc.getFmin().intValue() + 1), 
                 featureLoc.getFmax().intValue());
            sequenceBuffer.append(subSeq);
          }
        }
        else
        {
          char[] subSeq = sequence.getSequence();
          sequenceBuffer.append(subSeq);
          if(sequence.getPhase() != null)
            phase = sequence.getPhase().intValue();
        }
        
        final String seqStr;
        if(showPeptideSequence)
          seqStr = AminoAcidSequence.getTranslation(
              sequenceBuffer.toString().substring(phase), true).toString();
        else
          seqStr = new String(sequenceBuffer.toString());
        
        seqs.add(new org.emboss.jemboss.editor.Sequence(ortho, seqStr));
      }
      catch(NullPointerException npe)
      {
        JOptionPane.showMessageDialog(null, 
            "Cannot get the sequence for "+ortho,
            "Warning", JOptionPane.WARNING_MESSAGE);
        npe.printStackTrace();
      }
    }
    
    org.emboss.jemboss.editor.AlignJFrame ajFrame =
          new org.emboss.jemboss.editor.AlignJFrame(seqs);
    ajFrame.setVisible(true);
  }
  
  /**
   * Called by AbstractMatchTable.updateQualifier()
   */
  protected String updateQualifierString(final int row)
  {
    final String type;
    if( ((String)getTable().getValueAt(row, getColumnIndex(ROW_TYPE_COL))).equals("O") )
      type = MatchPanel.ORTHOLOG;
    else
      type = MatchPanel.PARALOG;
    
    StringBuffer orthologStr = new StringBuffer(
        (String)getTable().getValueAt(row, getColumnIndex(GENE_COL))+
        " link="+
        (String)getTable().getValueAt(row, getColumnIndex(LINK_COL))+
        " type="+type);            // ortholog link
    orthologStr.append(";");
    
    String clusterName = (String)getTable().getValueAt(row, getColumnIndex(CLUSTER_NAME_COL));
    if(clusterName != null && !clusterName.equals(""))
      orthologStr.append("cluster_name="+clusterName+ ";" ); // cluster name

    String product = (String)getTable().getValueAt(row, getColumnIndex(PRODUCT_COL));
    if(product != null && !product.equals(""))
      orthologStr.append("product="+product+ ";" );
    
    
    String matchName = (String)getTable().getValueAt(row, getColumnIndex(MATCH_NAME_COL));
    if(matchName != null && !matchName.equals(""))
      orthologStr.append("match_name="+matchName+ ";" ); // match name
    
    orthologStr.append("rank="+row);
    
    return orthologStr.toString();
  }
  
  /**
   * Override AbstractMatchTables.getOtherValues(). 
   * Removes all values apart from those that are not displayed in
   * the table, i.e. depending on whether a cluster table is displayed
   * or the ortholog/paralog table.
   * @param origQualifier
   * @return
   */
  protected StringVector getOtherValues(final Qualifier origQualifier)
  {
    StringVector values = origQualifier.getValues();
    Vector removeValues = new Vector();
    for(int j=0; j<values.size(); j++)
    {
      String value = (String) values.elementAt(j);
      if( (value.indexOf("cluster_name=") > -1 && showCluster) ||
          (value.indexOf("cluster_name=") < 0  && !showCluster) )
        removeValues.add(value);
    }
    values.removeAll(removeValues);
    
    return values;
  }
  
  /**
   * Returns true if the qualifier name matches the row type, e.g.
   * orthologous_to / paralogous_to
   * @param qualifierName
   * @param row
   * @return
   */
  protected boolean isRowOfType(String qualifierName, int row)
  {
    String rowType = (String)getTable().getValueAt(row, getColumnIndex(ROW_TYPE_HIDE_COL));
    
    if(rowType.equals(qualifierName))
      return true;
    return false;
  }

  /**
   * Check whether ortholog/paralog qualifier string exists in a StringVector for that qualifier.
   * If the StringVector contains the hit, description return true.
   * @param qualStr
   * @param qualStringVector
   * @return
   */
  public static boolean containsStringInStringVector(final String qualStr, 
                                                     final StringVector qualStringVector)
  {
    final StringVector orth1 = StringVector.getStrings(qualStr, ";");
    final String clusterName1 = ArtemisUtils.getString(orth1, "cluster");
    //final String rank1 = ArtemisUtils.getString(orth1, "rank");
    String value1 = (String)orth1.get(0);
    int index;

    if(!clusterName1.equals("")) // if a part of a cluster
    {
      final String clusterElements1[] = value1.split(", ");
      final List findAll = Arrays.asList(clusterElements1);
      
      for(int i=0; i<findAll.size(); i++)
      {
        final String findMe = (String)findAll.get(i);
        boolean found = false;
        for(int j=0; j<qualStringVector.size(); j++)
        {
          String thisStr = (String)qualStringVector.get(j);
          StringVector orth2 = StringVector.getStrings(thisStr, ";");
          final String clusterName2 = ArtemisUtils.getString(orth2, "cluster_name");
          if(!clusterName1.equals(clusterName2))
            continue;
        
          String value2 = (String)orth2.get(0);
          if((index = value2.indexOf('=')) > -1)
            value2 = value2.substring(index+1);
          String clusterElements2[] = value2.split(", ");
        
          final List searchList = Arrays.asList(clusterElements2);

          if(searchList.contains(findMe))
            found = true;
        }
        if(!found)
          return false;
      }
      return true;
    }
    
    for(int i=0; i<qualStringVector.size(); i++)
    {
      String thisStr = (String)qualStringVector.get(i);
      
      StringVector orth2 = StringVector.getStrings(thisStr, ";");
      
      if(orth1.size() != orth2.size())
        continue;
      
      String value2 = (String)orth2.get(0);
      if((index = value2.indexOf('=')) > -1)
        value2 = value2.substring(index+1);
      
      if(!clusterName1.equals("") && !ArtemisUtils.getString(orth2, "cluster_name").equals(""))
        System.out.println(value1+"  ==>  "+value2);
      // ortholog/paralog/cluster
      if(value1.indexOf(value2) < 0 &&
         value2.indexOf(value1) < 0 )
        continue;
      
      // cluster name
      final String clusterName2 = ArtemisUtils.getString(orth2, "cluster_name");
      if(!clusterName1.equals(clusterName2))
        continue;
      
      // rank
      /*
      final String rank2 = ArtemisUtils.getString(orth2, "rank");
      if(!rank1.equals(rank2))
        continue;
      */
      
      // description
      /*
      if( orth1.size() > 1 && orth2.size() > 1 &&
          !((String)orth1.get(1)).equals((String)orth2.get(1)) )
        continue;
      */
      
      return true; 
    }
    return false;
  }
  
  /**
   * Renderer for the Ortholog cells
   */
  private class OrthologRenderer extends DefaultTableCellRenderer
  {  
    /** */
    private static final long serialVersionUID = 1L;
    private int minHeight = -1;
    
    private final JLabel gene = new JLabel();
    private final JLabel link = new JLabel();
    private final JLabel type = new JLabel();
    private final JLabel program = new JLabel();
    private final JLabel symbol = new JLabel();
    private final JLabel organism = new JLabel();
    private final JLabel product = new JLabel();
    private final JTextArea descriptionTextArea = new JTextArea();
    private final JLabel clusterName = new JLabel();
    private final JLabel matchName = new JLabel();
    private final JLabel buttRemove = new JLabel("X");
    private final JLabel buttView = new JLabel("VIEW");
    private Color fgColor = new Color(139,35,35);
    private Color fgLinkColor = Color.BLUE;
    
    public OrthologRenderer() 
    {
      gene.setForeground(Color.BLUE);
      gene.setOpaque(true);

      clusterName.setOpaque(true);
      organism.setOpaque(true);
      link.setOpaque(true);
      product.setOpaque(true);
      
      descriptionTextArea.setLineWrap(true);
      descriptionTextArea.setWrapStyleWord(true);

      buttRemove.setOpaque(true);
      buttRemove.setText("X");
      
      Font font = getFont().deriveFont(Font.BOLD);
      buttRemove.setFont(font);
      buttRemove.setToolTipText("REMOVE");
      buttRemove.setHorizontalAlignment(SwingConstants.CENTER);
      
      buttView.setOpaque(true);
      buttView.setFont(font);
      buttView.setHorizontalAlignment(SwingConstants.CENTER);
      
      symbol.setOpaque(true);
      symbol.setFont(font);
      symbol.setHorizontalAlignment(SwingConstants.CENTER);
    }
    

    public Component getTableCellRendererComponent(
        JTable table,
        Object value,
        boolean isSelected,
        boolean hasFocus,
        final int row,
        final int column) 
    {
      Component c = null;
      String text = null;
      if(value != null)
        text = (String)value;
      Dimension dim;

      TableColumn tableCol;
      if(column == getColumnIndex(GENE_COL))
      {
        String geneStr[] = text.split(":");
        if(geneStr.length > 1)
          gene.setText(geneStr[1]);
        
        if(isSelected) 
        {
          gene.setForeground(fgLinkColor);
          gene.setBackground(table.getSelectionBackground());
        } 
        else
        {
          gene.setForeground(fgLinkColor);
          gene.setBackground(UIManager.getColor("Button.background"));
        }
        
        c = gene;
      }
      else if(column == getColumnIndex(ORGANISM_COL))
      {
        organism.setText(text);
        c = organism;
      }
      else if(column == getColumnIndex(PRODUCT_COL))
      {
        product.setText(text);
        c = product;
      }
      else if(column == getColumnIndex(LINK_COL))
      {
        link.setText(text);
        c = link;
      }
      else if(column == getColumnIndex(CLUSTER_NAME_COL))
      {
        clusterName.setText(text);

        tableCol = table.getColumnModel().getColumn(column);
        clusterName.setSize(tableCol.getWidth(), table
            .getRowHeight(row));

        dim = clusterName.getPreferredSize();
        minHeight = Math.max(minHeight, dim.height);
        c = clusterName;
      }
      else if(column == getColumnIndex(MATCH_NAME_COL))
      {
        matchName.setText(text);
        c = matchName;
      }
      else if(column == getColumnIndex(ROW_TYPE_HIDE_COL))
      {
        type.setText(text);
        c = type;
      }
      else if(column == getColumnIndex(PROGRAM_COL))
      {
        program.setText(text);
        c = program;
      }
      else if(column == getColumnIndex(ROW_TYPE_COL))
      {
        symbol.setText(text);
        
        if(text.equals("O"))
          symbol.setForeground(fgColor);
        else
          symbol.setForeground(Color.GREEN);
        tableCol = table.getColumnModel().getColumn(column);
        symbol.setSize(tableCol.getWidth(), table
            .getRowHeight(row));

        dim = symbol.getPreferredSize();
        minHeight = Math.max(minHeight, dim.height);
        c = symbol;
      }
      else if(column == getColumnIndex(VIEW_BUTTON_COL))
      {
        if(isSelected) 
        {
          buttView.setForeground(fgColor);
          buttView.setBackground(table.getSelectionBackground());
        } 
        else
        {
          buttView.setForeground(fgColor);
          buttView.setBackground(UIManager.getColor("Button.background"));
        }
        c = buttView;
      }
      else if(column == getColumnIndex(REMOVE_BUTTON_COL))
      {
        if(isSelected) 
        {
          buttRemove.setForeground(fgColor);
          buttRemove.setBackground(table.getSelectionBackground());
        } 
        else
        {
          buttRemove.setForeground(fgColor);
          buttRemove.setBackground(UIManager.getColor("Button.background"));
        }
        c = buttRemove;
      }
      else
      {
        throw new RuntimeException("invalid column! " + column +
                                   " " + text);
      }

      // adjust row height for columns with multiple lines
      if(column < 3)
      {
        if(table.getRowHeight(row) < minHeight)
          table.setRowHeight(row, minHeight);

        minHeight = -1;
      }
      
      // highlight on selection
      if(isSelected) 
        c.setBackground(table.getSelectionBackground()); 
      else
        c.setBackground(Color.white);
      
      return c;
    }
  }

  public class OrthoParalogValueComparator implements Comparator
  {
    public int compare(Object o1, Object o2)
    {
      final String value1 = (String)o1;
      final String value2 = (String)o2;
      
      StringVector values1 = StringVector.getStrings((String)value1, ";");
      String rank1 = ArtemisUtils.getString(values1, "rank");
      if(!rank1.equals(""))
      {
        if(rank1.startsWith("rank=") || rank1.startsWith("rank "))
          rank1 = rank1.substring(5);
      }
      else
        rank1 = "0";
      
      StringVector values2 = StringVector.getStrings((String)value2, ";");
      String rank2 = ArtemisUtils.getString(values2, "rank");
      if(!rank2.equals(""))
      {
        if(rank2.startsWith("rank=") || rank2.startsWith("rank "))
          rank2 = rank2.substring(5);
      }
      else
        rank2 = "0";
      
      
      return (new Integer(rank1)).compareTo(new Integer(rank2));
    }   
  }
  
  public class FeatureNamePredicate implements FeaturePredicate
  {
    private String uniqueName;
    
    public FeatureNamePredicate(final String uniqueName)
    {
      this.uniqueName = uniqueName;
    }
    
    public boolean testPredicate(Feature feature)
    {
      try
      {
        String featureId = (String)feature.getQualifierByName("ID").getValues().get(0);
        if(featureId.equals(uniqueName))
          return true;
      }
      catch(InvalidRelationException e){}
      
      return false;
    }
    
  }
}