/* GFF3AttributeBuilder.java
 * *
 * This file is part of Artemis
 *
 * Copyright (C) 2014 Genome Research Limited
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

package uk.ac.sanger.artemis.io;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.lang.StringBuilder;

import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.GFF3AttributeAggregator;
import uk.ac.sanger.artemis.util.StringVector;

public class GFF3AttributeBuilder {
  private HashMap<String, String>                  attrs        = new HashMap<String, String>();
  private HashMap<String, String>                  mappings     = new HashMap<String, String>();
  private HashMap<String, String>                  glue         = new HashMap<String, String>();
  private HashSet<String>                          ignores      = new HashSet<String>();
  private HashMap<String, String>                  clones       = new HashMap<String, String>();
  private HashMap<String, GFF3AttributeAggregator> aggs         = new HashMap<String, GFF3AttributeAggregator>();
  private HashSet<String>                          reserved     = new HashSet<String>(
                                                                    13);
  private static final GFF3AttributeAggregator defaultAgg = new GFF3AttributeAggregator() {
    @Override
    public String process(StringVector values) {
      StringBuilder buffer = new StringBuilder();
      if (values != null && values.size() > 0) {
        for (int value_index = 0; value_index < values.size(); ++value_index) {
          final String this_value = GFF3Encoder.encode(values.elementAt(value_index));
          if (value_index > 0 && value_index < (values.size())) {
            buffer.append(",");
          }
          buffer.append(this_value);
        }
      }
      return buffer.toString();
    }
  };

  final String  reserved_a[] = { "ID",
      "Name", "Alias", "Parent", "Derives_from", "Target", "Gap", "Note",
      "Dbxref", "Ontology_term", "Start_range", "End_range", "Is_circular" };

  public GFF3AttributeBuilder() {
    for (String s : reserved_a) {
      reserved.add(s);
    }
  }

  public void setMapping(String attr, String mapsto) {
    /* XXX: make sure only one level of mapping is used */
    mappings.put(attr, mapsto);
  }

  public void unsetMapping(String attr, String mapsto) {
    mappings.remove(attr);
  }

  public void setGlue(String attr, String glue_str) {
    glue.put(attr, glue_str);
  }

  public void unsetGlue(String attr) {
    glue.remove(attr);
  }

  public void setClone(String attr, String glue_str) {
    clones.put(attr, glue_str);
  }

  public void unsetClone(String attr) {
    clones.remove(attr);
  }

  public void setAggregator(String attr, GFF3AttributeAggregator agg) {
    aggs.put(attr, agg);
  }

  public void unsetAggregator(String attr, GFF3AttributeAggregator agg) {
    aggs.remove(attr);
  }

  public void ignore(String attr) {
    ignores.add(attr);
  }

  public void unignore(String attr) {
    ignores.remove(attr);
  }

  public void add(String attr, String val) {
    StringVector v = new StringVector(val);
    add(attr, v);
  }

  public void add(String attr, StringVector val) {
    String origAttr = attr;
    ArrayList<String> targetAttrs = new ArrayList<String>();
    // expand attributes
    if (clones.containsKey(attr))
      targetAttrs.add(clones.get(attr));
    if (mappings.containsKey(attr)) {
      attr = mappings.get(attr);
      targetAttrs.add(attr);
      if (clones.containsKey(attr))
        targetAttrs.add(clones.get(attr));
    } else {
      targetAttrs.add(attr);
    }
    // drop attributes with empty values
    if (val.size() == 1
        && val.elementAt(0).replaceAll("\\s+", "").equals(""))
      return;
    // process expanded list of attributes
    for (String this_attr : targetAttrs) {
      String aggregatedVal;
      // do we have an aggregator for this type?
      if (aggs.containsKey(origAttr)) {
        GFF3AttributeAggregator agg = aggs.get(origAttr);
        aggregatedVal = agg.process(val);
      } else if (aggs.containsKey(this_attr)) {
        GFF3AttributeAggregator agg = aggs.get(this_attr);
        aggregatedVal = agg.process(val);
      } else {
        aggregatedVal = defaultAgg.process(val);
      }
      // do not add empty values
      if (aggregatedVal == null)
        return;
      // append or set?
      if (attrs.containsKey(this_attr)) {
        String this_val = attrs.get(this_attr),
            this_glue = " ";
        if (glue.containsKey(this_attr))
          this_glue = glue.get(this_attr);
        this_val = this_val + this_glue + aggregatedVal;
        attrs.put(this_attr, this_val);
      } else {
        attrs.put(this_attr, aggregatedVal);
      }
    }
  }

  private String decapitalize(String line) {
    if (!reserved.contains(line)
        && Character.toUpperCase(line.charAt(0)) == line.charAt(0)) {
      return Character.toLowerCase(line.charAt(0)) + line.substring(1);
    } else {
      return line;
    }
  }

  public String get(String attr) throws EntryInformationException {
    if (mappings.containsKey(attr)) {
      attr = mappings.get(attr);
    }
    if (attrs.containsKey(attr)) {
      return attrs.get(attr);
    } else {
      throw new EntryInformationException("empty attribute value for " + attr);
    }
  }

  private Comparator<String> comparator = new Comparator<String>() {
    // make sure 'ID' is always at the beginning of the attribute list
    public int compare(String o1, String o2) {
      if (o1.equals("ID"))
        return -1;
      if (o2.equals("ID"))
        return 1;
      return o1.compareTo(o2);
    }
  };

  public String toString() {
    StringBuilder b = new StringBuilder();
    int i = 0;
    ArrayList<String> sortedAttrs = new ArrayList<String>(attrs.keySet());
    Collections.sort(sortedAttrs, comparator);
    for (String key : sortedAttrs) {
      if (!ignores.contains(key)) {
        String value = attrs.get(key); 
        if (i++ != 0)
          b.append(";");
        b.append(decapitalize(key) + "=" + value);
      }
    }
    return b.toString();
  }
}
