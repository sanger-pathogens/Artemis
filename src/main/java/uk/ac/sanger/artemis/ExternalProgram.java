/* ExternalProgram.java
 *
 * created: Tue Jan 26 1999
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998-2005  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/ExternalProgram.java,v 1.21 2009-08-17 15:49:25 tjc Exp $
 **/

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.components.ProgressBarFrame;
import uk.ac.sanger.artemis.components.filetree.RemoteFileNode;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.j2ssh.FileTransferProgressMonitor;
import uk.ac.sanger.artemis.j2ssh.FTProgress;

import java.io.*;
import java.text.*;
import java.util.Hashtable;
import java.util.Enumeration;

/**
 *  Each object of this class represents one external executable or script,
 *  and contains methods for invoking it.
 *
 *  @author Kim Rutherford
 *  @version $Id: ExternalProgram.java,v 1.21 2009-08-17 15:49:25 tjc Exp $
 **/

public class ExternalProgram 
{
  public static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(ExternalProgram.class);
  
  final public static int AA_PROGRAM = 0;
  final public static int DNA_PROGRAM = 1;
  final public static int APPLICATION = 2;

  /**
   *  The string to put at the end of the file that stores the counter used
   *  when creating the protein files to run the search on.  The file contains
   *  a comment line and then a line with a single integer on it.  This file
   *  is updated by setFileNumber() and read with getFileNumber().
   **/
  final private static String file_counter_filename = "file_number_counter";

  /**
   *  The name of the external program that was passed to the constructor.
   **/
  private String name;

  /** The default options to pass to this program. */
  private String program_options;

  /** One of AA_PROGRAM, DNA_PROGRAM, APPLICATION  */
  private int program_type;

  /**
   *  Create a new ExternalProgram object for the program with given name.
   *  @param name The name of the program.
   *  @param search_options The default options to use when running the
   *    program.
   *  @param program_type The type of external program to run(eg. AA_PROGRAM)
   **/
  public ExternalProgram(final String name,
                         final String program_options,
                         final int program_type) 
  {
    this.name = name;
    this.program_options = program_options;
    this.program_type = program_type;
  }

  /**
   *  Return the program_type that was passed to the constructor.
   **/
  public int getType() 
  {
    return program_type;
  }

  /**
   *  Run this program with the options that were set with setOptions().
   *  @param features The program will be run for each of these features.
   *  @param logger The log for errors, STDOUT and STDERR of the external
   *    program.
   *  @return An ExternalProgramThread object that can be used to monitor the
   *    external program after it starts.
   *  @exception IOException Thrown if an IO error occur while trying to run
   *    the program(eg the program could not be found).
   *  @exception ExternalProgramException Thrown if there is a non-IO error
   *    while attempting to run the program.
   *  @exception ReadOnlyException Thrown if the one features is in a
   *    read-only entry.
   *  @exception EntryInformationException Thrown if there is no qualifier
   *    that can be used to store the output filename
   **/
  public ExternalProgramMonitor run(final FeatureVector features,
                                    final Logger logger)
      throws IOException, ExternalProgramException, EntryInformationException,
             ReadOnlyException 
  {
    final StringVector sequence_file_names = new StringVector();

    // sequence_file_names will be set by prepareRun()
    final File file_of_filenames = prepareRun(features, sequence_file_names);

    try
    {
      final String[] arguments;
      switch (program_type)
      {
      case DNA_PROGRAM:
        // fall through
      case AA_PROGRAM:
        arguments = new String[]
        { file_of_filenames.getPath(), getProgramOptions() };
        break;
      case APPLICATION:
        arguments = new String[]
        { file_of_filenames.getPath(), };
        break;
      default:
        throw new Error("internal error - unknown program type");
      }

      final Process process = startProgram("run_" + getRealName(), arguments);

      //
      new ProgressBarFrame(1, getName());
      return new ProcessMonitor(process, getName(), logger);
    }
    catch (SecurityException e)
    {
      // re-throw as an ExternalProgramException
      throw new ExternalProgramException("SecurityException while running "
          + getName() + ": " + e.getMessage());
    }
  }
  

  /**
   *  Write sequence files for each of the given features and add a
   *  /something_file qualifier
   *  @param features Files will be written for each of these features.
   *  @param sequence_file_names The names of each sequence files will be
   *    returned in this Vector.
   *  @return a File representing a file that contains the names of each of
   *    the newly created files.
   **/
  private File prepareRun(final FeatureVector features,
                          final StringVector sequence_file_names)
      throws IOException, ExternalProgramException, EntryInformationException,
             ReadOnlyException 
  {
    final String new_qualifier_name = getName() + "_file";

    if(getType() != APPLICATION) 
    {
      // do this first so that we don't get an exception half way through the
      // list
      for(int i = 0 ; i < features.size() ; ++i) 
      {
        final Feature this_feature = features.elementAt(i);

        if(this_feature.getEntry().isReadOnly()) 
          throw new ReadOnlyException();

        final EntryInformation entry_information =
          this_feature.getEntry().getEntryInformation();

        if(!entry_information.isValidQualifier(this_feature.getKey(),
                                                 new_qualifier_name)) 
        {
          final String message = this_feature.getKey() + " cannot have " +
                                 new_qualifier_name + " as a qualifier";
          throw new EntryInformationException(message);
        }
      }
    }

    if(features.size() == 0) 
      return null;

    final NumberFormat number_format = NumberFormat.getNumberInstance();

    number_format.setMaximumIntegerDigits(5);
    number_format.setMinimumIntegerDigits(5);
    number_format.setGroupingUsed(false);

    // save the directory of the first feature we see - used for writing the
    // file of filenames
    File first_directory = null;

    // stores the number of features in each directory
    final Hashtable<File, Long> feature_count_hash = new Hashtable<File, Long>();

    for(int i = 0 ; i < features.size() ; ++i) 
    {
      final Feature this_feature = features.elementAt(i);
      final Entry this_feature_entry = this_feature.getEntry();

      File base_directory;

      base_directory = getBaseDirectoryFromEntry(this_feature_entry);

      if(base_directory == null) 
        base_directory = new File(".");

      if(first_directory == null) 
        first_directory = base_directory;

      if(feature_count_hash.containsKey(base_directory)) 
      {
        final long new_number =
         ((Long) feature_count_hash.get(base_directory)).longValue() + 1;

        feature_count_hash.put(base_directory, new Long(new_number));
      } 
      else 
        feature_count_hash.put(base_directory, new Long(1));
    }

    // store the file number to use for the next sequence file - the key is
    // a File containing the directory of there Entry that contains the
    // Feature and the value is the next file number to use

    Entry entry = features.elementAt(0).getEntry();
    RemoteFileNode node = null;

    if(((DocumentEntry)entry.getEMBLEntry()).getDocument()
                            instanceof RemoteFileDocument)
    {
      RemoteFileDocument nodeDoc =
               (RemoteFileDocument)(((DocumentEntry)entry.getEMBLEntry()).getDocument());
      node = nodeDoc.getRemoteFileNode();
    }  
    
    final Hashtable<File, Long> file_number_hash = new Hashtable<File, Long>();

    for(final Enumeration<File> e = feature_count_hash.keys() ;
         e.hasMoreElements() ;) 
    {
      final File directory = e.nextElement();

      final long old_file_number = getFileNumber(directory, node);

      final long feature_count =
       ((Long) feature_count_hash.get(directory)).longValue();

      file_number_hash.put(directory, new Long(old_file_number));

      setFileNumber(directory, old_file_number + feature_count, node);
    }

    // write the sequences out
    for(int i = 0 ; i < features.size() ; ++i) 
    {
      final Feature this_feature = features.elementAt(i);
      final Entry this_feature_entry = this_feature.getEntry();

      // the directory where we will write the program results.  defaults to
      // the current directory
      File base_directory;

      base_directory = getBaseDirectoryFromEntry(this_feature_entry);

      if(base_directory == null) 
        base_directory = new File(".");

      final File program_directory = new File(base_directory, getName());

      makeDirectory(program_directory);

      final String entry_name;

      {
        final String test_name = this_feature_entry.getName();
        if(test_name == null) 
          entry_name = "no_name";
        else 
          entry_name = test_name;
      }


      final long old_file_number =
       ((Long) file_number_hash.get(base_directory)).longValue();

      file_number_hash.put(base_directory, new Long(old_file_number + 1));

      final String new_file_name =
        getName() + File.separatorChar + entry_name + ".seq." +
        number_format.format(old_file_number);

      final File write_file;

      final String new_file_name_full =
        new File(base_directory,
                  File.separatorChar + new_file_name).getPath();

      sequence_file_names.add(new_file_name_full);

      write_file = new File(new_file_name_full);

      final Writer writer = new FileWriter(write_file);

      switch(program_type) 
      {
        case DNA_PROGRAM:
          this_feature.writeBasesOfFeature(writer);
          break;
        case AA_PROGRAM:
          this_feature.writeAminoAcidsOfFeature(writer);
          break;
        case APPLICATION:
          this_feature.writeNative(writer);
          break;
        default:
          throw new Error("internal error - unknown program type");
      }

      writer.close();

      if(program_type != APPLICATION) 
      {
        setFeatureQualifier(this_feature, new_qualifier_name, new_file_name);
        
        if( ((DocumentEntry)entry.getEMBLEntry()).getDocument()
             instanceof DatabaseDocument )
        {
          setProteinFeatureQualifier(this_feature, new_qualifier_name, new_file_name);
        }
      }
    }

    // write the sequence file names to a file and then run the program with
    // that file as an argument
    final File file_of_filenames =
      new File(new File(first_directory, getName()),
                getName() + "_" +
                "file_of_filenames." +
               (getFileNumber(first_directory, node) - 1));

    final Writer filenames_writer = new FileWriter(file_of_filenames);

    final PrintWriter filenames_printwriter =
      new PrintWriter(filenames_writer);

    for(int i = 0 ; i < sequence_file_names.size() ; ++i) 
      filenames_printwriter.println(sequence_file_names.elementAt(i));

    filenames_printwriter.close();
    filenames_writer.close();

    logger4j.info(name + " input file: "+file_of_filenames.getCanonicalPath());


    return file_of_filenames;
  }
  
  /**
   * Add a qualifier link to the program results
   * @param this_feature
   * @param new_qualifier_name
   * @param new_file_name
   * @throws ReadOnlyException
   * @throws EntryInformationException
   */
  private void setFeatureQualifier(final Feature this_feature, 
                                   final String new_qualifier_name,
                                   final String new_file_name) 
          throws ReadOnlyException, EntryInformationException
  {
    uk.ac.sanger.artemis.io.Qualifier new_qualifier;
    
    if(getName().startsWith("fast") || getName().startsWith("blast"))
    {
      new_qualifier = this_feature.getQualifierByName(new_qualifier_name);
      String db = getProgramOptions().trim();

      int ind;
      if((ind = db.indexOf(' ')) > -1)
        db.substring(0, ind);

      if(new_qualifier == null)
        new_qualifier = new uk.ac.sanger.artemis.io.Qualifier(
            new_qualifier_name, db + ":" + new_file_name + ".out");
      else
      {
        // search for existing 'db:' values
        final StringVector values = new_qualifier.getValues();
        for(int j = 0; j < values.size(); j++)
        {
          String value = (String) values.get(j);
          if(value.startsWith(db + ":"))
            values.remove(value);
        }

        new_qualifier.addValues(values);
        values.add(db + ":" + new_file_name + ".out");
        new_qualifier = new uk.ac.sanger.artemis.io.Qualifier(
            new_qualifier_name, values);
      }
    }
    else
      new_qualifier = new uk.ac.sanger.artemis.io.Qualifier(
          new_qualifier_name, new_file_name + ".out");
    
    this_feature.setQualifier(new_qualifier);
  }
  
  /**
   * In database mode add a link on the polypeptide to the results file.
   * @param this_feature
   * @param new_qualifier_name
   * @param new_file_name
   */
  private void setProteinFeatureQualifier(final Feature this_feature, 
                                          final String new_qualifier_name,
                                          final String new_file_name)
  {
    try
    {
      if(this_feature.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL)
          && this_feature.getEmblFeature() instanceof GFFStreamFeature)
      {
        GFFStreamFeature gffFeature = (GFFStreamFeature)(this_feature.getEmblFeature());
        ChadoCanonicalGene chado_gene = gffFeature.getChadoGene();
        String transcriptName = chado_gene.getTranscriptFromName(GeneUtils.getUniqueName(gffFeature));
        Feature protein = 
          (Feature)chado_gene.getProteinOfTranscript(transcriptName).getUserData();
        setFeatureQualifier(protein, new_qualifier_name, new_file_name);
      }
    }
    catch(Exception e)
    {
      logger4j.error(e.getMessage());
    }
  }

  /**
   *  Return the name of this ExternalProgram, as passed to the constructor.
   **/
  public String getName() 
  {
    return name;
  }

  /**
   *  Start a new external program with exec.
   *  @param name the name of the program to start.
   *  @param arguments the arguments to pass to the new program. (can be null)
   *  @return A Process object for the new program.
   **/
  public static Process startProgram(final String name,
                                     final String [] arguments)
      throws SecurityException, ExternalProgramException, IOException 
  {

    // try to get the program from the artemis directory/jar file
    InputStream code_stream = ExternalProgram.class.getResourceAsStream(name);

    // try the etc directory
    if(code_stream == null) 
      code_stream = ExternalProgram.class.getResourceAsStream("/etc/" + name);

    if(code_stream == null) 
    {
      // the code isn't in the artemis directory/jar file so just call exec()

      final String [] real_arguments;

      if(arguments == null) 
        real_arguments = new String [1];
      else
      {
        real_arguments = new String [arguments.length + 1];

        for(int i = 0 ; i < arguments.length ; ++i) 
          real_arguments[i + 1] = arguments[i];
      }

      real_arguments[0] = name;

      return Runtime.getRuntime().exec(real_arguments);
    } 
    else
    {
      final String [] sh_arguments;

      if(arguments == null) 
        sh_arguments = new String [2];
      else 
      {
        sh_arguments = new String [arguments.length + 2];

        for(int i = 0 ; i < arguments.length ; ++i) 
          sh_arguments[i + 2] = arguments[i];
      }

      sh_arguments[0] = "/bin/sh";
      sh_arguments[1] = "-s";

      final Process process = Runtime.getRuntime().exec(sh_arguments);
      final OutputStream out_stream = process.getOutputStream();

      final int BUFFER_SIZE = 10000;
      final byte [] buffer = new byte [BUFFER_SIZE];

      int read_size = 0;

      while((read_size = code_stream.read(buffer)) != -1) 
        out_stream.write(buffer, 0, read_size);

      code_stream.close();
      out_stream.close();

      return process;
    }
  }


  /**
   *  Read the number to be used when creating the next output file.  This
   *  number is appended to file names to attempt to make them unique.  The
   *  number is read from directory + "/" + prog_name + "/" +
   *  file_counter_filename.  If the file doesn't exist, it is created and
   *  initialised.
   **/
  private long getFileNumber(final File directory,
                             final RemoteFileNode node)
      throws IOException 
  {
    try 
    {
      final Reader file_reader;
      if(node == null)
        file_reader =
          new FileReader(new File(directory, File.separatorChar +
                                  getName() + File.separatorChar +
                                  file_counter_filename));
      else
      {
        String fs  = "/";    // assume ssh to unix server

        String dir = node.getRootDir()+ fs +
                     node.getFullName();
        int index  = dir.lastIndexOf(fs);
        dir = dir.substring(0,index) + fs + getName();

        FileTransferProgressMonitor monitor =
                   new FileTransferProgressMonitor(null);
        FTProgress progress = monitor.add(node.getFile());
        byte[] contents = node.getFileContents(progress, dir+fs+ 
                                               file_counter_filename);
        monitor.close();

        if(contents == null)
        {
          logger4j.debug("getFileNumber() creating "+dir+
                                fs+file_counter_filename);

          node.mkdir(dir);
          return setFileNumber(directory, 1, node);
        }

        file_reader = new StringReader(new String(contents));

        logger4j.debug("getFileNumber()\n"+new String(contents));
      }

      final BufferedReader reader = new BufferedReader(file_reader);

      // the first line should be a comment
      final String comment_line = reader.readLine();

      if(comment_line == null || comment_line.length() == 0) 
        return setFileNumber(directory, guessNumber(directory), node);

      if(!comment_line.startsWith("#")) 
        return setFileNumber(directory, guessNumber(directory), node);

      final String number_line = reader.readLine();

      try
      {
        return Integer.parseInt(number_line);
      }
      catch(NumberFormatException e) 
      {
        return setFileNumber(directory, guessNumber(directory), node);
      }
      finally
      {
        file_reader.close();
        reader.close();
      }
    }     
    catch(FileNotFoundException e)
    {
      // create a new file_number_counter file
      return setFileNumber(directory, guessNumber(directory), node);
    }
    catch(IOException e) 
    {
      return setFileNumber(directory, guessNumber(directory), node);
    }
  }

  /**
   *  Write the given number to the file_number file in the given directory
   * (eg. directory + "/blastp/" + new_file_number).
   **/
  protected long setFileNumber(final File directory,
                               final long new_file_number, final RemoteFileNode node)
      throws IOException 
  {
    makeDirectory(new File(directory, getName()));

    File local_file = new File(directory, File.separatorChar +
                               getName() + File.separatorChar +
                               file_counter_filename);

    final FileWriter file_writer = new FileWriter(local_file);
    final PrintWriter print_writer = new PrintWriter(file_writer);

    print_writer.println("# the file is machine generated - do not edit");
    print_writer.println(new_file_number);
    print_writer.close();
    file_writer.close();
    
    if(node != null)
    {
      final String fs = "/";  // assume unix ssh server
      String dir = node.getRootDir()+ fs +
                   node.getFullName();
      int index  = dir.lastIndexOf(fs);
      dir = dir.substring(0,index) + fs +
                         getName() + fs;

      logger4j.debug("setFileNumber() "+
                     local_file.getCanonicalPath()+" --> "+dir);

      node.put(dir, local_file, null, true);     
    }

    return new_file_number;
  }

  /**
   *  Create a directory path from the given directory(if it does not
   *  exist).
   **/
  protected void makeDirectory(final File directory)
      throws IOException 
  {
    if(directory.exists()) 
    {
      if(directory.isDirectory()) 
      {
        // all is ok
        if(directory.canWrite()) 
          return;
        else 
          throw new IOException("Cannot write to: " +
                                 directory.getAbsolutePath());
      }
      else 
        throw new IOException(directory.getAbsolutePath() +
                               " is not a directory");
    }
    else 
      directory.mkdirs();
  }

  /**
   *  Try to read the current file names in the results directory and guess
   *  the number we are up to.  This can be slow.
   **/
  private long guessNumber(final File directory) 
  {
    try 
    {
      long biggest_so_far = -1;

      final String[] file_names = new File(directory, getName()).list();

      if(file_names == null) 
        return 1;

      for(int i = 0 ; i < file_names.length ; ++i) 
      {
        final int last_dot = file_names[i].lastIndexOf('.');

        final String number_string;

        if(last_dot == -1) 
          number_string = file_names[i];
        else 
          number_string = file_names[i].substring(last_dot + 1);

        try 
        {
          final long number = Long.parseLong(number_string);

          if(number > biggest_so_far) 
            biggest_so_far = number;
        }
        catch(NumberFormatException e) {}
      }

      if(biggest_so_far == -1) 
        return 1;
      else 
        return biggest_so_far;
    }
    catch(SecurityException e) 
    {
      // give up
      return 0;
    }
  }

  /**
   *  Return the options that will be used when the program is run.
   **/
  public String getProgramOptions() 
  {
    if(program_options.equals("-")) 
      return "";
    else 
      return program_options;
  }

  /**
   *  Set the options that will be used when the program is run.
   **/
  public void setProgramOptions(final String program_options) 
  {
    this.program_options = program_options;
  }

  /**
   *  Return the name of this program without modifiers(the bits after the
   *  "+" in the name).  eg. if the name is "blastp+go" this method will be
   *  "blastp"
   **/
  private String getRealName() 
  {
    final int plus_index = getName().indexOf('+');

    if(plus_index > 0) 
      return getName().substring(0, plus_index);
    else 
      return getName();
  }

  /**
   *  Return the dirtectory that the given entry was read from.
   **/
  private File getBaseDirectoryFromEntry(final Entry entry) 
  {
    final uk.ac.sanger.artemis.io.Entry embl_entry = entry.getEMBLEntry();

    if(embl_entry instanceof DocumentEntry) 
    {
      final DocumentEntry document_entry =(DocumentEntry) embl_entry;

      if(document_entry.getDocument() instanceof FileDocument) 
      {
        final FileDocument file_document =
         (FileDocument) document_entry.getDocument();

        if(file_document.getFile().getParent() != null) 
          return new File(file_document.getFile().getParent());
      }
    }
    if(((DocumentEntry)entry.getEMBLEntry()).getDocument()
       instanceof RemoteFileDocument ||
       ((DocumentEntry)entry.getEMBLEntry()).getDocument()
       instanceof DatabaseDocument)
      return new File(System.getProperty("user.dir"));

    return null;
  }


}
