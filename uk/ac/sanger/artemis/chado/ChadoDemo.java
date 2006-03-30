/* 
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2006  Genome Research Limited
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

package uk.ac.sanger.artemis.chado;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.GridLayout;
import java.sql.SQLException;
import java.util.List;
import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPasswordField;
import javax.swing.JTextField;
import javax.swing.JList;
import javax.swing.JScrollPane;
import javax.swing.JPanel;

public class ChadoDemo
{
  /** JDBC DAO */
  private JdbcDAO jdbcDAO = null;
  /** iBatis DAO */
  private IBatisDAO connIB = null;
  private String location;
  private JPasswordField pfield;
  
  public ChadoDemo(final int feature_id)
  {
    try
    {
      setLocation(true);
      getGene(feature_id);
    }
    catch(java.net.ConnectException ce)
    {
      ce.printStackTrace();
    }
    catch(SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...\n"+
          sqlExp.getMessage(), 
          "SQL Error",
          JOptionPane.ERROR_MESSAGE);
      sqlExp.printStackTrace();
    }
    
  }
  
  public void getGene(final int feature_id)
              throws java.net.ConnectException, SQLException
  {

    int index = location.indexOf('=') + 1;
    String schema = location.substring(index);

    ChadoDAO dao = getDAO();

    List schemas = dao.getSchema();

    final String port_options[] = { "FIND", "CANCEL" };

    Vector v_schemas = new Vector(schemas);
    v_schemas.add(0, "All");

    JPanel panel = new JPanel(new BorderLayout());
    JList schema_list = new JList(v_schemas);
    schema_list.setSelectedValue("All", true);
    JScrollPane jsp = new JScrollPane(schema_list);
    panel.add(jsp, BorderLayout.EAST);

    JTextField gene_text = new JTextField(20);
    panel.add(gene_text, BorderLayout.SOUTH);

    int select = JOptionPane.showOptionDialog(null, panel, "Find Gene",
                                            JOptionPane.DEFAULT_OPTION, 
                                            JOptionPane.WARNING_MESSAGE, null,
                                            port_options, port_options[0]);

    if(select == 1)
      return;

    String search_gene = gene_text.getText();
    schema = (String)schema_list.getSelectedValue();
    
    List schema_search;
    if(schema.equalsIgnoreCase("All"))
      schema_search = schemas;
    else
    {
      schema_search = new Vector();
      schema_search.add(schema);
    }
    
    List featureList = dao.getFeature(search_gene, schema_search);
    
    for(int i = 0; i < featureList.size(); i++)
    {
      ChadoFeature feature = (ChadoFeature) featureList.get(i);
      int fmin = feature.getFmin() + 1;
      int fmax = feature.getFmax();

      System.out.print(feature.getSchema() + " " + fmin + " " + fmax);
      System.out.print(" \t" + feature.getType_id());
      System.out.print(" \t" + feature.getProp_type_id());
      System.out.print(" \t" + feature.getStrand());
      System.out.print(" \t" + feature.getUniquename());
      System.out.print(" \t" + feature.getTimelastmodified().toString());
      System.out.println(" " + Integer.toString(feature.getId()));
    }

  }
  
  /**
   * 
   * Get the data access object (DAO).
   * 
   * @return data access object
   * 
   */
 private ChadoDAO getDAO()
    throws java.net.ConnectException, SQLException
 {
   if(System.getProperty("ibatis") == null)
   {
     if(jdbcDAO == null)
      jdbcDAO = new JdbcDAO(location, pfield); 
     return jdbcDAO;
   }
   else
   {
     if(connIB == null)
       connIB = new IBatisDAO(pfield);
     return connIB;
   }
 }
 
 /**
  * 
  * Set the database location as:
  * jdbc:postgresql://localhost:13001/chadoCVS?user=es2
  * 
  */
 protected boolean setLocation(final boolean prompt_user)
 {
   Container bacross = new Container();
   bacross.setLayout(new GridLayout(6, 2, 5, 5));

   JLabel lServer = new JLabel("Server : ");
   bacross.add(lServer);
   JTextField inServer = new JTextField("localhost");
   bacross.add(inServer);

   JLabel lPort = new JLabel("Port : ");
   bacross.add(lPort);
   JTextField inPort = new JTextField("5432");
   bacross.add(inPort);

   JLabel lDB = new JLabel("Database : ");
   bacross.add(lDB);
   JTextField inDB = new JTextField("chado");
   bacross.add(inDB);

   JLabel lUser = new JLabel("User : ");
   bacross.add(lUser);
   JTextField inUser = new JTextField("afumigatus");
   bacross.add(inUser);

   JLabel lpasswd = new JLabel("Password : ");
   bacross.add(lpasswd);
   pfield = new JPasswordField(16);
   bacross.add(pfield);

   // given -Dchado=localhost:port/dbname?username
   if(System.getProperty("chado") != null)
   {
     String db_url = System.getProperty("chado").trim();
     int index;
     if((index = db_url.indexOf(":")) > -1)
     {
       inServer.setText(db_url.substring(0, index));
       int index2;
       if((index2 = db_url.indexOf("/")) > -1)
       {
         inPort.setText(db_url.substring(index + 1, index2));
         int index3;
         if((index3 = db_url.indexOf("?")) > -1)
         {
           inDB.setText(db_url.substring(index2 + 1, index3));
           inUser.setText(db_url.substring(index3 + 1));

           /*
            * if(!prompt_user) { location = "jdbc:postgresql://"
            * +inServer.getText().trim()+ ":" +inPort.getText().trim()+ "/"
            * +inDB.getText().trim()+ "?user=" +inUser.getText().trim(); return
            * true; }
            */
         }
       }
     }
   }

   int n = JOptionPane.showConfirmDialog(null, bacross,
                                         "Enter Database Address",
                                         JOptionPane.OK_CANCEL_OPTION,
                                         JOptionPane.QUESTION_MESSAGE);
   if(n == JOptionPane.CANCEL_OPTION)
     return false;

   location = "jdbc:postgresql://" + 
              inServer.getText().trim() + ":" +
              inPort.getText().trim() + "/" +
              inDB.getText().trim() + "?user=" +
              inUser.getText().trim();

   return true;
 }
 
  public static void main(String args[])
  { 
    new ChadoDemo(Integer.parseInt(args[0]));
  }
}