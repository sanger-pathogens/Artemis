/* DbSqlConfig.java
 *
 * created: 2005
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2005  Genome Research Limited
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

import com.ibatis.sqlmap.client.SqlMapClient;
import com.ibatis.sqlmap.client.SqlMapClientBuilder;
import com.ibatis.common.resources.Resources;

import uk.ac.sanger.artemis.util.DatabaseLocationParser;

import javax.swing.JPasswordField;
import java.util.Properties;
import java.io.Reader;

/**
 *
 * Initialises iBatis configuration.
 *
 */
public class DbSqlConfig
{

  private SqlMapClient sqlMap;

  /**
   *
   * Initialises iBatis configuration by reading the
   * artemis_sqlmap/chado_iBatis_config.xml file and the
   * database location properties defined by the system property, 
   * <i>e.g.</i> -Dchado=localhost:2997/chado?tjc
   * 
   * It is especially important that the -Dchado argument is not prefixed with
   * jdbc:postgresql:// as that is handled in iBatis
   *
   */
  public void init(JPasswordField fpasswd)
  {
    try
    {
     String resource = "artemis_sqlmap/chado_iBatis_config.xml";
     Reader reader = Resources.getResourceAsReader(resource);

     Properties properties = null;
     if(System.getProperty("chado") != null)
     {
       String url = System.getProperty("chado");
       DatabaseLocationParser dlp = new DatabaseLocationParser(url);
       properties = new Properties();

       properties.put("chado",dlp.getUnprefixedURL());

       properties.put("username",dlp.getUsername());
       
       if(fpasswd != null && fpasswd.getPassword().length > 0)
         properties.put("password", new String(fpasswd.getPassword()));
     }

     sqlMap =  SqlMapClientBuilder.buildSqlMapClient(reader, 
                                               properties);
    }
    catch(Exception e)
    {
       // If you get an error at this point, it doesnt matter what it was.  It is going to be
       // unrecoverable and we will want the app to blow up hard so we are aware of the
       // problem.  You should always log such errors and re-throw them in such a way that
       // you can be made immediately aware of the problem.
       e.printStackTrace();
       throw new RuntimeException("Error initializing DbSqlConfig class.  Cause: "  + e);
    }
  }

  
  public void init2(JPasswordField fpasswd)
  {
    try
    {
     String resource = "artemis_sqlmap/chado_iBatis_config.xml";
     Reader reader = Resources.getResourceAsReader(resource);
     sqlMap =  SqlMapClientBuilder.buildSqlMapClient(reader);
    }
    catch(Exception e)
    {
       e.printStackTrace();
       throw new RuntimeException("Error initializing DbSqlConfig class.  Cause: "  + e);
    }
  }
  
  public SqlMapClient getSqlMapInstance()
  {
    return sqlMap;
  }

}

