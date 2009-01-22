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

import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class EmbossDigestParser {
    
    private static final Logger logger = Logger.getLogger(EmbossDigestParser.class);
    
    String embossDir; 
    
    private List<String> digests = Collections.emptyList();
        
    public List<String> getDigests() {
        return digests;
    }

    public void afterPropertiesSet() {
        try {
            parseDigests();
        } catch (IOException e) {
            logger.warn(String.format("Emboss directory %s does't exist",embossDir));
        }
    }
    
    private void parseDigests() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(new File(embossDir+"REBASE/embossre.enz")));
        String line;
        digests = new ArrayList<String>();
        while ((line = br.readLine()) != null) {
            if (line.startsWith("#")) {
                continue;
            }
            int space = line.indexOf("\t");
            if (space>1) {
                digests.add(line.substring(0, space));
            } else {
                System.err.println("Couldn't get enzyme name from '"+line+"'");
            }
        }
    }

    public void setEmbossDir(String embossDir) {
        this.embossDir = embossDir;
    }

}
