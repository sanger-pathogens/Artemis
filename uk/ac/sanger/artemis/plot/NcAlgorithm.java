/* NcAlgorithm.java
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2004  Genome Research Limited
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

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.*;

/**
 *
 * Nc (Wright 1990)
 * Wright F (1990) The 'effective number of codons' used in a gene. Gene 87:23-9
 *
 *  @author Derek Gatherer
 *  original version 09-09-03
 *  revised 01-12-04
 *  division by zero bugremoved 08-12-04
 **/

public class NcAlgorithm extends BaseAlgorithm 
{
  /**
   *  Create a new NcAlgorithm object.
   *  @param strand The Strand to do the calculation on.
   **/
  public NcAlgorithm(final Strand strand) 
  {
    super(strand, makeName(strand), "Nc");
    setScalingFlag(true);
  }

  /**
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).  If the start/end pair
   *    doesn't give a multiple of three bases end is moved down so that it is
   *    a multiple of three.
   *  @param values The one return value for this algorithm is returned in
   *    this array.
   **/
  public void getValues(int start, int end, final float [] values) 
  {
    if(getStrand().isForwardStrand())  // rather than isRevCompDisplay()
    {
//    System.out.println("isRevCompDisplay does not activate here");
    }
    else
    {
      final int new_end =
        getStrand().getBases().getComplementPosition(start);
      final int new_start =
        getStrand().getBases().getComplementPosition(end);

      end = new_end;
      start = new_start;
    }
    final String sub_sequence;
    
    // add 1 or 2 if necessary to make the range a multiple of 3
    if(getStrand().isForwardStrand())
      end -= (end - start + 1) % 3;
    else
      start += (end - start + 1) % 3;

    try 
    {
      sub_sequence = getStrand().getSubSequence(new Range(start, end));
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    final char[] sequence_raw;
    sequence_raw = sub_sequence.toCharArray();

    float[][][][] p = new float [4][4][4][4]; // holds p
    float[][][][] obs_value = new float [4][4][4][4];
    float[] p2_phe = new float [3];  // hold p2 totals for each amino
    float[] p2_leu = new float [3];
    float[] p2_ile = new float [3];
    float[] p2_val = new float [3];
    float[] p2_ser = new float [3];
    float[] p2_pro = new float [3];
    float[] p2_thr = new float [3];
    float[] p2_ala = new float [3];
    float[] p2_tyr = new float [3];
    float[] p2_his = new float [3];
    float[] p2_gln = new float [3];
    float[] p2_asn = new float [3];
    float[] p2_lys = new float [3];
    float[] p2_asp = new float [3];
    float[] p2_glu = new float [3];
    float[] p2_cys = new float [3];
    float[] p2_arg = new float [3];
    float[] p2_gly = new float [3];
    float[] F2 = new float [3]; // hold Faa for each syn group
    float[] F3 = new float [3];
    float[] F4 = new float [3];
    float[] F6 = new float [3];
    float[] Nc = new float [3]; // hold Nc for each frame
    float NC2 = 0;
    float NC3 = 0;
    float NC4 = 0;
    float NC6 = 0;
    

// initialise all arrays to zero, or 1 to prevent infinity (check paper is this valid?)

    for(int x = 0 ; x < 4 ; ++x) {
      for(int y = 0 ; y < 4 ; ++y) {
        for(int z = 0 ; z < 4 ; ++z) {
          for(int a = 0 ; a < 4 ; ++a) {
            obs_value[x][y][z][a] = 0;
            p[x][y][z][a] = 0;
          }
        }
      }
    }

    char this_f_base = 0;
    char next_f_base = 0;
    char last_f_base = 0;
    int this_f_base_index = 0;
    int next_f_base_index = 0;
    int last_f_base_index = 0;

    for(int i = 0 ; i < sequence_raw.length - 5 ; i+=3) 
    {
      for(int frame = 0; frame < 3; frame++) 
      {
        this_f_base = sequence_raw[i + frame];
        next_f_base = sequence_raw[i + 1 + frame];
        last_f_base = sequence_raw[i + 2 + frame];

        this_f_base_index = Bases.getIndexOfBase(this_f_base);
        next_f_base_index = Bases.getIndexOfBase(next_f_base);
        last_f_base_index = Bases.getIndexOfBase(last_f_base);
        
        // ignore Ns
        if(this_f_base_index < 4 && next_f_base_index < 4 && last_f_base_index < 4)
          ++obs_value[this_f_base_index][next_f_base_index][last_f_base_index][(frame+start)%3];
      } 
    }
    
 // having collected the observed values in all 3 frames, formulate p and p2
 
    for(int frame = 0; frame < 3; ++frame) 
    {
 // Phe 001 000
      if(obs_value[0][0][0][frame]>0 || obs_value[0][0][1][frame]>0) 
      {
        p[0][0][0][frame] = obs_value[0][0][0][frame]/(obs_value[0][0][0][frame]+obs_value[0][0][1][frame]);
        p[0][0][1][frame] = obs_value[0][0][1][frame]/(obs_value[0][0][0][frame]+obs_value[0][0][1][frame]);
      }
      p2_phe[frame] = (float)Math.pow(p[0][0][0][frame],2)+(float)Math.pow(p[0][0][1][frame],2);
 // Leu = 002 003 100 101 102 103
      if(obs_value[0][0][2][frame]>0 || obs_value[0][0][3][frame]>0 || obs_value[1][0][0][frame]>0 ||
         obs_value[1][0][1][frame]>0 || obs_value[1][0][2][frame]>0 || obs_value[1][0][3][frame]>0) 
      {
        p[0][0][2][frame] = obs_value[0][0][2][frame]/(obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame]);
        p[0][0][3][frame] = obs_value[0][0][3][frame]/(obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame]);
        p[1][0][0][frame] = obs_value[1][0][0][frame]/(obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame]);
        p[1][0][1][frame] = obs_value[1][0][1][frame]/(obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame]);
        p[1][0][2][frame] = obs_value[1][0][2][frame]/(obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame]);
        p[1][0][3][frame] = obs_value[1][0][3][frame]/(obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame]);
      }
      p2_leu[frame] = (float)Math.pow(p[0][0][2][frame],2)+(float)Math.pow(p[0][0][3][frame],2)+(float)Math.pow(p[1][0][0][frame],2)+
          (float)Math.pow(p[1][0][1][frame],2)+(float)Math.pow(p[1][0][2][frame],2)+(float)Math.pow(p[1][0][3][frame],2);
// Ile = 200 201 202
      if(obs_value[2][0][0][frame]>0 || obs_value[2][0][1][frame]>0 || obs_value[2][0][2][frame]>0) 
      {
        p[2][0][0][frame] = obs_value[2][0][0][frame]/(obs_value[2][0][0][frame]+obs_value[2][0][1][frame]+obs_value[2][0][2][frame]);
        p[2][0][1][frame] = obs_value[2][0][1][frame]/(obs_value[2][0][0][frame]+obs_value[2][0][1][frame]+obs_value[2][0][2][frame]);
        p[2][0][2][frame] = obs_value[2][0][2][frame]/(obs_value[2][0][0][frame]+obs_value[2][0][1][frame]+obs_value[2][0][2][frame]);
      }
      p2_ile[frame] = (float)Math.pow(p[2][0][0][frame],2)+(float)Math.pow(p[2][0][1][frame],2)+(float)Math.pow(p[2][0][2][frame],2);
// Val = 300 301 302 303
      if(obs_value[3][0][0][frame]>0 || obs_value[3][0][1][frame]>0 || obs_value[3][0][2][frame]>0 ||
         obs_value[3][0][3][frame]>0) 
      {
        p[3][0][0][frame] = obs_value[3][0][0][frame]/(obs_value[3][0][0][frame]+obs_value[3][0][1][frame]+obs_value[3][0][2][frame]+obs_value[3][0][3][frame]);
        p[3][0][1][frame] = obs_value[3][0][1][frame]/(obs_value[3][0][0][frame]+obs_value[3][0][1][frame]+obs_value[3][0][2][frame]+obs_value[3][0][3][frame]);
        p[3][0][2][frame] = obs_value[3][0][2][frame]/(obs_value[3][0][0][frame]+obs_value[3][0][1][frame]+obs_value[3][0][2][frame]+obs_value[3][0][3][frame]);
        p[3][0][3][frame] = obs_value[3][0][3][frame]/(obs_value[3][0][0][frame]+obs_value[3][0][1][frame]+obs_value[3][0][2][frame]+obs_value[3][0][3][frame]);
      }
      p2_val[frame] = (float)Math.pow(p[3][0][0][frame],2)+(float)Math.pow(p[3][0][1][frame],2)+(float)Math.pow(p[3][0][2][frame],2)+(float)Math.pow(p[3][0][3][frame],2);
// Ser = 010 011 012 013 230 231
      if(obs_value[0][1][0][frame]>0 || obs_value[0][1][1][frame]>0 || obs_value[0][1][2][frame]>0 ||
         obs_value[0][1][3][frame]>0 || obs_value[2][3][0][frame]>0 || obs_value[2][3][1][frame]>0) 
      {
        p[0][1][0][frame] = obs_value[0][1][0][frame]/(obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
         +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame]);
        p[0][1][1][frame] = obs_value[0][1][1][frame]/(obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
         +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame]);
        p[0][1][2][frame] = obs_value[0][1][2][frame]/(obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
         +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame]);
        p[0][1][3][frame] = obs_value[0][1][3][frame]/(obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
        +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame]);
        p[2][3][0][frame] =  obs_value[2][3][0][frame]/(obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
        +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame]);
        p[2][3][1][frame] = obs_value[2][3][1][frame]/(obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
        +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame]);
      }
      p2_ser[frame] = (float)Math.pow(p[0][1][0][frame],2)+(float)Math.pow(p[0][1][1][frame],2)+(float)Math.pow(p[0][1][2][frame],2)+
          (float)Math.pow(p[0][1][3][frame],2)+(float)Math.pow(p[2][3][0][frame],2)+(float)Math.pow(p[2][3][1][frame],2);
// Pro = 110 111 112 113
      if(obs_value[1][1][0][frame]>0 || obs_value[1][1][1][frame]>0 || obs_value[1][1][2][frame]>0 ||
         obs_value[1][1][3][frame]>0) 
      {
        p[1][1][0][frame] = obs_value[1][1][0][frame]/(obs_value[1][1][0][frame]+obs_value[1][1][1][frame]+obs_value[1][1][2][frame]+obs_value[1][1][3][frame]);
        p[1][1][1][frame] = obs_value[1][1][1][frame]/(obs_value[1][1][0][frame]+obs_value[1][1][1][frame]+obs_value[1][1][2][frame]+obs_value[1][1][3][frame]);
        p[1][1][2][frame] = obs_value[1][1][2][frame]/(obs_value[1][1][0][frame]+obs_value[1][1][1][frame]+obs_value[1][1][2][frame]+obs_value[1][1][3][frame]);
        p[1][1][3][frame] = obs_value[1][1][3][frame]/(obs_value[1][1][0][frame]+obs_value[1][1][1][frame]+obs_value[1][1][2][frame]+obs_value[1][1][3][frame]);
      }
      p2_val[frame] = (float)Math.pow(p[1][1][0][frame],2)+(float)Math.pow(p[1][1][1][frame],2)+(float)Math.pow(p[1][1][2][frame],2)+(float)Math.pow(p[1][1][3][frame],2);
// Thr = 210 211 212 213
      if(obs_value[2][1][0][frame]>0 || obs_value[2][1][1][frame]>0 || obs_value[2][1][2][frame]>0 ||
         obs_value[2][1][3][frame]>0) 
      {
        p[2][1][0][frame] = obs_value[2][1][0][frame]/(obs_value[2][1][0][frame]+obs_value[2][1][1][frame]+obs_value[2][1][2][frame]+obs_value[2][1][3][frame]);
        p[2][1][1][frame] = obs_value[2][1][1][frame]/(obs_value[2][1][0][frame]+obs_value[2][1][1][frame]+obs_value[2][1][2][frame]+obs_value[2][1][3][frame]);
        p[2][1][2][frame] = obs_value[2][1][2][frame]/(obs_value[2][1][0][frame]+obs_value[2][1][1][frame]+obs_value[2][1][2][frame]+obs_value[2][1][3][frame]);
        p[2][1][3][frame] = obs_value[2][1][3][frame]/(obs_value[2][1][0][frame]+obs_value[2][1][1][frame]+obs_value[2][1][2][frame]+obs_value[2][1][3][frame]);
      }
      p2_thr[frame] = (float)Math.pow(p[2][1][0][frame],2)+(float)Math.pow(p[2][1][1][frame],2)+(float)Math.pow(p[2][1][2][frame],2)+(float)Math.pow(p[2][1][3][frame],2);
// Ala = 310 311 312 313
      if(obs_value[3][1][0][frame]>0 || obs_value[3][1][1][frame]>0 || obs_value[3][1][2][frame]>0 ||
         obs_value[3][1][3][frame]>0) 
      {
        p[3][1][0][frame] = obs_value[3][1][0][frame]/(obs_value[3][1][0][frame]+obs_value[3][1][1][frame]+obs_value[3][1][2][frame]+obs_value[3][1][3][frame]);
        p[3][1][1][frame] = obs_value[3][1][1][frame]/(obs_value[3][1][0][frame]+obs_value[3][1][1][frame]+obs_value[3][1][2][frame]+obs_value[3][1][3][frame]);
        p[3][1][2][frame] = obs_value[3][1][2][frame]/(obs_value[3][1][0][frame]+obs_value[3][1][1][frame]+obs_value[3][1][2][frame]+obs_value[3][1][3][frame]);
        p[3][1][3][frame] = obs_value[3][1][3][frame]/(obs_value[3][1][0][frame]+obs_value[3][1][1][frame]+obs_value[3][1][2][frame]+obs_value[3][1][3][frame]);
      }
      p2_ala[frame] = (float)Math.pow(p[3][1][0][frame],2)+(float)Math.pow(p[3][1][1][frame],2)+(float)Math.pow(p[3][1][2][frame],2)+(float)Math.pow(p[3][1][3][frame],2);
// Tyr = 020 021
      if(obs_value[0][2][0][frame]>0 || obs_value[0][2][1][frame]>0) 
      {
        p[0][2][0][frame] = obs_value[0][2][0][frame]/(obs_value[0][2][0][frame]+obs_value[0][2][1][frame]);
        p[0][2][1][frame] = obs_value[0][2][1][frame]/(obs_value[0][2][0][frame]+obs_value[0][2][1][frame]);
      }
      p2_tyr[frame] = (float)Math.pow(p[0][2][0][frame],2)+(float)Math.pow(p[0][2][1][frame],2);
// His = 120 121
      if(obs_value[1][2][0][frame]>0 || obs_value[1][2][1][frame]>0) 
      {
        p[1][2][0][frame] = obs_value[1][2][0][frame]/(obs_value[1][2][0][frame]+obs_value[1][2][1][frame]);
        p[1][2][1][frame] = obs_value[1][2][1][frame]/(obs_value[1][2][0][frame]+obs_value[1][2][1][frame]);
      }
      p2_his[frame] = (float)Math.pow(p[1][2][0][frame],2)+(float)Math.pow(p[1][2][1][frame],2);
// Gln = 122 123
      if(obs_value[1][2][2][frame]>0 || obs_value[1][2][3][frame]>0) 
      {
        p[1][2][2][frame] = obs_value[1][2][2][frame]/(obs_value[1][2][2][frame]+obs_value[1][2][3][frame]);
        p[1][2][3][frame] = obs_value[1][2][3][frame]/(obs_value[1][2][2][frame]+obs_value[1][2][3][frame]);
      }
      p2_gln[frame] = (float)Math.pow(p[1][2][2][frame],2)+(float)Math.pow(p[1][2][3][frame],2);
// Asn = 220 221
      if(obs_value[2][2][0][frame]>0 || obs_value[2][2][1][frame]>0) 
      {
        p[2][2][0][frame] = obs_value[2][2][0][frame]/(obs_value[2][2][0][frame]+obs_value[2][2][1][frame]);
        p[2][2][1][frame] = obs_value[2][2][1][frame]/(obs_value[2][2][0][frame]+obs_value[2][2][1][frame]);
      }
      p2_asn[frame] = (float)Math.pow(p[2][2][0][frame],2)+(float)Math.pow(p[2][2][1][frame],2);
// Lys = 222 223
      if(obs_value[2][2][2][frame]>0 || obs_value[2][2][3][frame]>0) 
      {
        p[2][2][2][frame] = obs_value[2][2][2][frame]/(obs_value[2][2][2][frame]+obs_value[2][2][3][frame]);
        p[2][2][3][frame] = obs_value[2][2][3][frame]/(obs_value[2][2][2][frame]+obs_value[2][2][3][frame]);
      }
      p2_lys[frame] = (float)Math.pow(p[2][2][2][frame],2)+(float)Math.pow(p[2][2][3][frame],2);
// Asp = 320 321
      if(obs_value[3][2][0][frame]>0 || obs_value[3][2][1][frame]>0) 
      {
        p[3][2][0][frame] = obs_value[3][2][0][frame]/(obs_value[3][2][0][frame]+obs_value[3][2][1][frame]);
        p[3][2][1][frame] = obs_value[3][2][1][frame]/(obs_value[3][2][0][frame]+obs_value[3][2][1][frame]);
      }
      p2_asp[frame] = (float)Math.pow(p[3][2][0][frame],2)+(float)Math.pow(p[3][2][1][frame],2);
// Glu = 322 323   
      if(obs_value[3][2][2][frame]>0 || obs_value[3][2][3][frame]>0) 
      {
        p[3][2][2][frame] = obs_value[3][2][2][frame]/(obs_value[3][2][2][frame]+obs_value[3][2][3][frame]);
        p[3][2][3][frame] = obs_value[3][2][3][frame]/(obs_value[3][2][2][frame]+obs_value[3][2][3][frame]);
      }
      p2_glu[frame] = (float)Math.pow(p[3][2][2][frame],2)+(float)Math.pow(p[3][2][3][frame],2);
//Cys = 030 031
      if(obs_value[0][3][0][frame]>0 || obs_value[0][3][1][frame]>0) 
      {
        p[0][3][0][frame] = obs_value[0][3][0][frame]/(obs_value[0][3][0][frame]+obs_value[0][3][1][frame]);
        p[0][3][1][frame] = obs_value[0][3][1][frame]/(obs_value[0][3][0][frame]+obs_value[0][3][1][frame]);
      }
      p2_cys[frame] = (float)Math.pow(p[0][3][0][frame],2)+(float)Math.pow(p[0][3][1][frame],2);
// Arg = 130 131 132 133 232 233
      if(obs_value[1][3][0][frame]>0 || obs_value[1][3][1][frame]>0 || obs_value[1][3][2][frame]>0 ||
         obs_value[1][3][3][frame]>0 || obs_value[2][3][2][frame]>0 || obs_value[2][3][3][frame]>0) 
      {
        p[1][3][0][frame] = obs_value[1][3][0][frame]/(obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame]);
        p[1][3][1][frame] = obs_value[1][3][1][frame]/(obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame]);
        p[1][3][2][frame] = obs_value[1][3][2][frame]/(obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame]);
        p[1][3][3][frame] = obs_value[1][3][3][frame]/(obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame]);
        p[2][3][2][frame] = obs_value[2][3][2][frame]/(obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame]);
        p[2][3][3][frame] = obs_value[2][3][3][frame]/(obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame]);
      }
      p2_arg[frame] = (float)Math.pow(p[1][3][0][frame],2)+(float)Math.pow(p[1][3][1][frame],2)+(float)Math.pow(p[1][3][2][frame],2)+
          (float)Math.pow(p[1][3][3][frame],2)+(float)Math.pow(p[2][3][2][frame],2)+(float)Math.pow(p[2][3][3][frame],2);
// Gly = 330 331 332 333
      if(obs_value[3][3][0][frame]>0 || obs_value[3][3][1][frame]>0 || obs_value[3][3][2][frame]>0 ||
         obs_value[3][3][3][frame]>0) 
      {
        p[3][3][0][frame] = obs_value[3][3][0][frame]/(obs_value[3][3][0][frame]+obs_value[3][3][1][frame]+obs_value[3][3][2][frame]+obs_value[3][3][3][frame]);
        p[3][3][1][frame] = obs_value[3][3][1][frame]/(obs_value[3][3][0][frame]+obs_value[3][3][1][frame]+obs_value[3][3][2][frame]+obs_value[3][3][3][frame]);
        p[3][3][2][frame] = obs_value[3][3][2][frame]/(obs_value[3][3][0][frame]+obs_value[3][3][1][frame]+obs_value[3][3][2][frame]+obs_value[3][3][3][frame]);
        p[3][3][3][frame] = obs_value[3][3][3][frame]/(obs_value[3][3][0][frame]+obs_value[3][3][1][frame]+obs_value[3][3][2][frame]+obs_value[3][3][3][frame]);
      }
      p2_gly[frame] = (float)Math.pow(p[3][1][0][frame],2)+(float)Math.pow(p[3][1][1][frame],2)+(float)Math.pow(p[3][1][2][frame],2)+(float)Math.pow(p[3][1][3][frame],2);

      p[2][0][3][frame] = 0;  // don't count Met
      p[0][3][3][frame] = 0;  // don't count Trp
      p[0][2][2][frame] = 0;  // don't count STOP
      p[0][2][3][frame] = 0;  // don't count STOP
      p[0][3][2][frame] = 0;  // don't count STOP

// having calculated p2 value, now do Faa, ie. F for each amino acid.
// various tedious routines to catch abset amino acids

      int F2_aa = 0;
      if(p2_cys[frame]>0) { F2_aa++; }
      if(p2_asp[frame]>0) { F2_aa++; }
      if(p2_glu[frame]>0) { F2_aa++; }
      if(p2_phe[frame]>0) { F2_aa++; }
      if(p2_his[frame]>0) { F2_aa++; }
      if(p2_lys[frame]>0) { F2_aa++; }
      if(p2_asn[frame]>0) { F2_aa++; }
      if(p2_gln[frame]>0) { F2_aa++; }
      if(p2_tyr[frame]>0) { F2_aa++; }
      if(F2_aa>0)
      {
        F2[frame] = (p2_cys[frame]+p2_asp[frame]+p2_glu[frame]+p2_phe[frame]+p2_his[frame]+p2_lys[frame]+
                    p2_asn[frame]+p2_gln[frame]+p2_tyr[frame])/F2_aa;
      }
      else { F2[frame] = 0; }
//      System.out.println("F2["+frame+"]="+F2[frame]+" cys: "+p2_cys[frame]+" asp: "+p2_asp[frame]+" glu: "+p2_glu[frame]
//      +" phe: "+p2_phe[frame]+" his: "+p2_his[frame]+" lys: "+p2_lys[frame]+" asn: "+p2_asn[frame]+" gln: "+p2_gln[frame]+" tyr: "+p2_tyr[frame]);

      int F6_aa = 0;
      if(p2_leu[frame]>0) { F6_aa++; }
      if(p2_arg[frame]>0) { F6_aa++; }
      if(p2_ser[frame]>0) { F6_aa++; }
      if(F6_aa>0)
      {
        F6[frame] = (p2_leu[frame]+p2_arg[frame]+p2_ser[frame])/F6_aa;
      }
      else { F6[frame] = 0; }
//      System.out.println("F6["+frame+"]="+F6[frame]);
      
      int F4_aa = 0;
      if(p2_gly[frame]>0) { F4_aa++; }
      if(p2_pro[frame]>0) { F4_aa++; }
      if(p2_thr[frame]>0) { F4_aa++; }
      if(p2_val[frame]>0) { F4_aa++; }
      if(p2_ala[frame]>0) { F4_aa++; }
      if(F4_aa>0)
      {
        F4[frame] = (p2_gly[frame]+p2_pro[frame]+p2_thr[frame]+p2_val[frame]+p2_ala[frame])/F4_aa;
      }
      else { F4[frame] = 0; }
//      System.out.println("F4["+frame+"]="+F4[frame]);
      
      if(p2_ile[frame]>0) 
      {
        F3[frame] = p2_ile[frame];
      } else {
           F3[frame]= (F2[frame] + F4[frame])/2;
      }
//      System.out.println("F3["+frame+"]="+F3[frame]);

// the line below gives the occasional div by zero error
//      Nc[frame] = (9/F2[frame])+(5/F4[frame])+(3/F6[frame])+(1/F3[frame])+2;
      if(F2[frame]>0)
      {
        NC2 = F2_aa/F2[frame];
      }
      else {  NC2 = 0; }
      if(F4[frame]>0)
      {
        NC4 = F4_aa/F4[frame];
      }
      else {  NC4 = 0; }
      if(F6[frame]>0)
      {
        NC6 = F6_aa/F6[frame];
      }
      else {  NC6 = 0; }
      if(F3[frame]>0)
      {
        NC3 = 1/F3[frame];
      }
      else {  NC3 = 0; }
      
      Nc[frame] = NC2+NC4+NC6+NC3+2;
//      System.out.println("Nc["+frame+"] is "+"NC2: "+NC2+" NC4: "+NC4+" NC6: "+NC6+" NC3: "+NC3+" + 2 = "+Nc[frame]);
      values[frame] = Nc[frame];
    }
  } 

  /**
   *  Return the number of values a call to getValues() will return - three
   *  in this case.
   **/
  public int getValueCount() 
  {
    return 3;
  }

  /**
   *  Return the default or optimal window size.
   *  @return null is returned if this algorithm doesn't have optimal window
   *    size.
   **/
  public Integer getDefaultWindowSize() 
  {
    final Integer super_window_size = super.getDefaultWindowSize();
    if(super_window_size != null) 
    {
      // the superclass version of getDefaultWindowSize() returns non-null
      // iff the user has set the window size in the options file
      return super_window_size;
    }
    return new Integer(500);
  }

  /**
   *  Return the default maximum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have maximum window
   *    size.
   **/
  public Integer getDefaultMaxWindowSize() 
  {
    final Integer super_max_window_size = super.getDefaultMaxWindowSize();
    if(super_max_window_size != null) 
    {
      // the superclass version of getDefaultMaxWindowSize() returns non-null
      // iff the user has set the max window size in the options file
      return super_max_window_size;
    }
    return new Integer(5000);
  }

  /**
   *  Return the default minimum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have minimum window
   *    size.
   **/
  public Integer getDefaultMinWindowSize() 
  {
    final Integer super_min_window_size = super.getDefaultMinWindowSize();
    if(super_min_window_size != null) 
    {
      // the superclass version of getDefaultMinWindowSize() returns non-null
      // iff the user has set the minimum window size in the options file
      return super_min_window_size;
    }
    return new Integer(24);
  }

  /**
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize(int window_size) 
  {
    if(window_size > 10) 
    {
      Integer step = new Integer(24);
//      Integer step = new Integer(window_size/10);
//      int step_int = step.intValue();
//      step_int+=(step_int % 3);  // step is multiple of 3
//      Integer step_out = new Integer(step_int);
//      return step_out;
      return step; 

    }
    else
      return null;
  }

  /**
   *  Return the maximum value of this algorithm.
   *  @return The maximum is 2.
   **/
  protected Float getMaximumInternal() 
  {
    return new Float(10000000);
  }

  /**
   *  Return the minimum value of this algorithm.
   *  @return The minimum is 0.
   **/
  protected Float getMinimumInternal() 
  {
    return new Float(0);
  }

  /**
   *  Return the average value of function over the whole strand.
   *  @return null is returned if this algorithm doesn't have an average or if
   *    the average can't be calculated.
   **/
  public Float getAverage() 
  {
    return null;
  }  
  
  private static String makeName(final Strand strand) 
  {
    if(strand.isForwardStrand()) 
      return "Effective Codon Number(Nc)";
    else 
      return "Reverse Effective Codon Number(Nc)";
  }
}  
