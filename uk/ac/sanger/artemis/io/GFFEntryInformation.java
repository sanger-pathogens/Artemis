/* GFFEntryInformation.java
 *
 * created: Thu Mar 30 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/GFFEntryInformation.java,v 1.6 2007-06-14 16:01:38 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.io.InputStream;
import java.io.IOException;
import java.util.Enumeration;
import java.util.Properties;

import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.Options;

/**
 *  An EntryInformation object for GFFDocumentEntry objects.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: GFFEntryInformation.java,v 1.6 2007-06-14 16:01:38 tjc Exp $
 **/

public class GFFEntryInformation extends SimpleEntryInformation 
{

  public GFFEntryInformation()
  {
    super(); 

    try
    {
      makeEntryInformation();
    }
    catch(QualifierInfoException e)
    {
      System.err.println("could not initialise the embl package: " +
                          e.getMessage());
      System.exit(1);
    }
    catch(IOException e)
    {
      System.err.println("could not initialise the embl package: " +
                          e.getMessage());
      System.exit(1);
    }
  }

 
  /**
   *  Return an EntryInformation object that is suitable for EMBL and GENBANK
   *  entries.
   **/
  private void makeEntryInformation()
      throws IOException, QualifierInfoException
  { 
    final InputStream feature_keys_stream =
      Options.class.getResourceAsStream("/etc/feature_keys_gff");

    final InputStream qualifier_types_stream =
      Options.class.getResourceAsStream("/etc/qualifier_types_gff");

    QualifierInfoVector qualifier_info_vector =
      readQualifierInfo(qualifier_types_stream, feature_keys_stream);

    for(int i = 0 ; i < qualifier_info_vector.size() ; ++i)
    {
      final QualifierInfo qualifier_info =
        qualifier_info_vector.elementAt(i);

      addQualifierInfo(qualifier_info);
    }
    
    //
    // include extra keys and qualifiers
    final QualifierInfoVector extra_qualifiers =
      Options.getOptions().getExtraGffQualifiers();
    final StringVector extra_keys =
      Options.getOptions().getOptionValues("extra_keys_gff");
    
    for(int i = 0 ; i < extra_keys.size() ; ++i)
    {
      final Key new_key = new Key((String)extra_keys.elementAt(i));
      //System.out.println(new_key.toString());
      addKey(new_key);
    }

    for(int i = 0 ; i < extra_qualifiers.size() ; ++i)
    {
      final QualifierInfo new_qualifier_info = extra_qualifiers.elementAt(i);
      //System.out.println(new_qualifier_info.getName());
      addQualifierInfo(new_qualifier_info);
    }

//  entry_information.setEMBLFormat(true);
  }

  /**
   *  Read the possible feature key and qualifier names and types from the two
   *  given streams (see etc/feature_keys and etc/qualifier_types for details
   *  on the formats).
   **/
  private static QualifierInfoVector
    readQualifierInfo(final InputStream qualifier_types_stream,
                      final InputStream feature_keys_stream)
      throws IOException 
  {

    final QualifierInfoVector return_vector = new QualifierInfoVector();

    Properties feature_properties = new Properties();
    final Properties qualifier_properties = new Properties();

    feature_properties.load(feature_keys_stream);
    qualifier_properties.load(qualifier_types_stream);

    // parse the feature_properties

    {
      final Properties new_feature_properties = new Properties();
      final Enumeration feature_enum = feature_properties.propertyNames();

      while(feature_enum.hasMoreElements()) 
      {
        String current_feature_name = (String) feature_enum.nextElement();

        final StringVector property_values =
          Options.getPropertyValues(feature_properties, current_feature_name);

        new_feature_properties.put(current_feature_name, property_values);
      }

      feature_properties = new_feature_properties;
    }

    final Enumeration qualifier_enum = qualifier_properties.propertyNames();

    while(qualifier_enum.hasMoreElements()) 
    {
      String current_qualifier_name = (String) qualifier_enum.nextElement();

      final StringVector current_qualifier_values =
        Options.getPropertyValues(qualifier_properties,
                                  current_qualifier_name);

      final boolean once_only =
        current_qualifier_values.elementAt(0).equals("yes");

      final String type_string = (String)current_qualifier_values.elementAt(1);

      // find the keys for which this qualifier name is valid or required

      final KeyVector valid_keys = new KeyVector();
      final KeyVector required_keys = new KeyVector();

      final Enumeration features_enum = feature_properties.propertyNames();

      while(features_enum.hasMoreElements()) 
      {
        final String current_key_string = (String)features_enum.nextElement();

        final Key current_key = new Key(current_key_string);

        final StringVector current_feature_qualifiers =
          (StringVector) feature_properties.get(current_key_string);

        if(current_feature_qualifiers.contains(current_qualifier_name)) 
          valid_keys.add(current_key);
        else
          if(current_feature_qualifiers.contains("@" +
                                                   current_qualifier_name)) 
          {
            valid_keys.add(current_key);
            required_keys.add(current_key);
          }
      }

      final int type = QualifierInfo.getQualifierTypeID(type_string);

      final QualifierInfo qualifier_info =
        new QualifierInfo(current_qualifier_name, type, valid_keys,
                           required_keys, once_only);

      return_vector.add(qualifier_info);
    }

    return return_vector;
  }
}
