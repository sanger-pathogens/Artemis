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


package uk.ac.sanger.ibatis;

import com.ibatis.sqlmap.client.SqlMapClient;
import com.ibatis.sqlmap.client.SqlMapClientBuilder;
import com.ibatis.common.resources.Resources;

import javax.swing.JPasswordField;
import java.util.Properties;
import java.io.Reader;

public class DbSqlConfig
{

  private static SqlMapClient sqlMap;

  public static void init(JPasswordField fpasswd)
  {
    try
    {
     String resource = "artemis_sqlmap/chado_iBatis_config.xml";
     Reader reader = Resources.getResourceAsReader(resource);

     Properties properties = null;
     if(System.getProperty("chado") != null)
     {
       String url = System.getProperty("chado");
       int index  = url.indexOf("?");
       int index2 = url.indexOf("user=");
       properties = new Properties();

       int index3 = url.indexOf("://");
       if(index3 < 0)
         index3 = 0;
       else
         index3 = index3+3;

       properties.put("chado", url.substring(index3,index)); 

       if(index2 < 0)
         properties.put("username", url.substring(index+1));
       else
         properties.put("username", url.substring(index2+5));

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

  public static SqlMapClient getSqlMapInstance()
  {
    return sqlMap;
  }

}

