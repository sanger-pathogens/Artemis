/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2013  Genome Research Limited
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
package uk.ac.sanger.artemis.components;

import static org.junit.Assert.assertTrue;

import org.junit.Before;
import org.junit.Test;

import java.awt.GraphicsEnvironment;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Properties;

/**
 * Tests for loading and writing project properties for
 * Artemis.
 */
public class ProjectPropertyTest
{
  private HashMap<String, HashMap<String, String>> projects;
  @Before
  public void loadProperties()
  {
    InputStream ins = 
        this.getClass().getClassLoader().getResourceAsStream("data/project.properties");
    
    try
    {
      final Properties projectProps = new Properties();
      projectProps.load(ins);
      ins.close();
      projects = ProjectProperty.getProjectMap(projectProps);
    }
    catch (IOException e)
    {
      org.junit.Assert.fail(e.getMessage());
    }
  }
  
  @Test
  public void projectContent1()
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.getLocalGraphicsEnvironment().isHeadless())
      return;

    assertTrue("Styphi project",  projects.containsKey("Styphi"));
    if(projects.containsKey("Styphi"))
    {
      HashMap<String, String> thisProject = projects.get("Styphi");
      assertTrue("Styphi project sequence",  thisProject.containsKey("sequence"));
      assertTrue("Styphi project annotation",  thisProject.containsKey("annotation"));
    }
  }
  
  @Test
  public void projectContent2()
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.getLocalGraphicsEnvironment().isHeadless())
      return;

    assertTrue("PF3D7 project",  projects.containsKey("PF3D7"));
    if(projects.containsKey("PF3D7"))
    {
      HashMap<String, String> thisProject = projects.get("PF3D7");
      assertTrue("Styphi project sequence",  thisProject.containsKey("sequence"));
      assertTrue("Styphi project bam",  thisProject.containsKey("bam"));
    }
  }
  
  @Test
  public void launch()
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.getLocalGraphicsEnvironment().isHeadless())
      return;

    try
    {
      ProjectProperty pp = new ProjectProperty();
    }
    catch(Exception e)
    {
      org.junit.Assert.fail(e.getMessage());
    }
  }
  
  @Test
  public void writeAndReadProperties()
  {
    // ignore if in headless mode with no x11
    if(GraphicsEnvironment.getLocalGraphicsEnvironment().isHeadless())
      return;

    try
    {
      final File tmpFile = File.createTempFile("artemis.props", ".tmp");
      ProjectProperty.writeProperties(tmpFile, projects);
      
      final InputStream ins = new FileInputStream(tmpFile);
      final Properties projectProps = new Properties();
      projectProps.load(ins);
      ins.close();
      final HashMap<String, HashMap<String, String>>
         tmpProjects = ProjectProperty.getProjectMap(projectProps);
      
      for (String p: projects.keySet())
      {
        assertTrue("Contains "+p,  tmpProjects.containsKey(p));
        
        HashMap<String, String> props = projects.get(p);
        HashMap<String, String> tmpProps = tmpProjects.get(p);
        for(final String key: props.keySet())
        {
          // test property key (e.g. title, sequence, bam)
          assertTrue("Contains property "+key,  tmpProps.containsKey(key));
          
          // test values
          assertTrue("Contains property "+key+" = "+props.get(key), 
              props.get(key).equals(tmpProps.get(key)));
        }
      }
    }
    catch(IOException e)
    {
      org.junit.Assert.fail(e.getMessage());
    }

  }
  
}
