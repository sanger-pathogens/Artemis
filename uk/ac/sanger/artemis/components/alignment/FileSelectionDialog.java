package uk.ac.sanger.artemis.components.alignment;


import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URL;
import java.util.List;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.components.StickyFileChooser;
import uk.ac.sanger.artemis.components.variant.VCFview;

/**
 * File selection panel to allow input of DNA sequences
 */
public class FileSelectionDialog extends JDialog
{
  private static final long serialVersionUID = 1L;
  
  private JPanel dialog = new JPanel(new GridBagLayout());
  private GridBagConstraints c = new GridBagConstraints();
  private int row = 0;
  private List<JTextField> bamFields = new Vector<JTextField>(30);
  private JTextField referenceField = new JTextField(30);
  private boolean statusOK = false;
  
  public FileSelectionDialog(String fileNames[])
  {
    for(int i=0; i<fileNames.length; i++)
    {
      JTextField bamField = new JTextField(fileNames[i]);
      bamFields.add(bamField);
    }
    statusOK = true;
  }
  
  /**
   * Constructor to display any given input files and options provided and
   * @param f
   * @param showReferenceOption
   * @param program - program name
   * @param fileType - type of file (e.g. bam, vcf)
   */
  public FileSelectionDialog(Frame f, boolean showReferenceOption, 
      final String program, final String fileType)
  {
    super(f, program+" :: Select Files", true);

    addBamField(fileType, null);
    
    JButton addMoreFiles = new JButton("Add More");
    addMoreFiles.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        addBamField(fileType, null);
      }
    });
    
    row += 100;
    c.gridx = 1;
    c.gridy = row;
    dialog.add(addMoreFiles, c);

    if(showReferenceOption)
    {
      c.gridy = ++row;
      c.gridx = 0;
      dialog.add(new JLabel(" Reference sequence file (optional): "), c);
      c.gridy = ++row;
      dialog.add(referenceField, c);
      JButton selectReference = new JButton("Select...");
      addActionListener(selectReference, referenceField, fileType);
      c.gridx = 1;
      dialog.add(selectReference, c);
    }
    
    JButton okButton = new JButton("OK");
    okButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        statusOK = true;
        FileSelectionDialog.this.dispose();
      }
    });
    c.gridy = ++row;
    dialog.add(okButton, c);
    getContentPane ().add (dialog, "South");
    
    row = 1;
    
    pack();
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));
    setVisible (true);
  }
  
  /**
   * Add a text field to the dialog for adding in a path
   * to a BAM file.
   */
  private void addBamField(String fileType, String fileName)
  {
    JTextField bamField = new JTextField(30);
    bamFields.add(bamField);
    
    if(fileName != null)
      bamField.setText(fileName);
    bamField.setPreferredSize(
        new Dimension(200,bamField.getPreferredSize().height));
    c.gridy = row;
    c.gridx = 0;
    c.anchor = GridBagConstraints.WEST;
    dialog.add(new JLabel(" "+fileType+" file: "), c);
    c.gridy = ++row;
    dialog.add(bamField, c);
    c.gridx = 1;
    JButton selectBam = new JButton("Select...");
    addActionListener(selectBam, bamField, fileType);
    dialog.add(selectBam, c);
    
    pack();
  }
  
  /**
   * Add action listener to a button.
   * @param fileSelectionButton
   * @param tfield
   * @param fileType
   */
  private void addActionListener(final JButton fileSelectionButton, 
                                 final JTextField tfield,
                                 final String fileType)
  {
    final javax.swing.filechooser.FileFilter filter =
        new javax.swing.filechooser.FileFilter()
    {
      public boolean accept(final File file) 
      {
        if(file.isDirectory() ||
           file.getName().matches(VCFview.VCFFILE_SUFFIX) ||
           file.getName().matches(BamView.BAM_SUFFIX))
          return true;
        return false;
      }

      public String getDescription() 
      {
        return "BAM / VCF files";
      }
    };
    
    fileSelectionButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        StickyFileChooser fileChooser = new StickyFileChooser();
        fileChooser.setMultiSelectionEnabled(true);
        fileChooser.setFileFilter(filter);
        
        int status = fileChooser.showOpenDialog(null);
        if(status == StickyFileChooser.CANCEL_OPTION)
          return;
        
        
        File[] files = fileChooser.getSelectedFiles();
        tfield.setText(files[0].getAbsolutePath());
        
        if(files.length > 1)
        {
          for(int i=1; i<files.length; i++)
            addBamField(fileType, files[i].getAbsolutePath());
        }
      }
    });
  }
  
  /**
   * Return true only if this is a list of filenames.
   * @param filename
   * @return
   */
  public static boolean isListOfFiles(String filename) 
  {
    if(filename.startsWith("http") || filename.startsWith("ftp"))
    {
      try
      {
        URL urlBamIndexFile = new URL(filename);
        InputStreamReader reader = new InputStreamReader(urlBamIndexFile.openStream());
        
        BufferedReader is = new BufferedReader(reader);
        String s = is.readLine();
        is.close();
        reader.close();
        
        if(s != null && (s.trim().startsWith("http") || s.trim().startsWith("ftp")))
          return true;
      }
      catch (IOException e){}

      return false;
    }

    try
    {
      File f = new File(filename);
      if (f.exists())
      {
        FileInputStream fstream = new FileInputStream(f);
        DataInputStream in = new DataInputStream(fstream);
        BufferedReader br = new BufferedReader(new InputStreamReader(in));

        String line = br.readLine().trim();
        br.close();
        fstream.close();
        
        if(line.startsWith("http") || line.startsWith("ftp"))
          return true;
        f = new File(line);
        if (f.exists())
          return true;
      }
    }
    catch (IOException e)
    {
      return false;
    }
    
    return false;
  }
  
  /**
   * Read a list of file names from a file.
   * @param filename
   * @return
   */
  public static List<String> getListOfFiles(String filename)
  {
    List<String> bamFiles = new Vector<String>();
    try
    {
      if(filename.startsWith("http"))
      {
        try
        {
          URL urlBamIndexFile = new URL(filename);
          InputStreamReader reader = new InputStreamReader(urlBamIndexFile.openStream()); 
          addFilesToListFromReader(reader, bamFiles);
          reader.close();
        }
        catch (IOException e){}
      }
      else
      {
        File f = new File(filename);
        if (f.exists())
        {
          FileInputStream fstream = new FileInputStream(f);
          InputStreamReader reader = new InputStreamReader(fstream); 
          addFilesToListFromReader(reader, bamFiles);
          reader.close();
        }
      }
    }
    catch (IOException e){}
    return bamFiles;
  }
  
  /**
   * Read a list of files in and add them to the collection.
   * @param in
   * @param bamFiles
   * @throws IOException
   */
  private static void addFilesToListFromReader(Reader in, List<String> bamFiles) 
          throws IOException
  {
    BufferedReader br = new BufferedReader(in);
    String line;
    while ((line = br.readLine()) != null)
    {
      line = line.trim();
      if(!line.equals(""))
        bamFiles.add(line);
    }
    br.close();
  }
  
  /**
   * Get the BAM or VCF files as a <code>List</code> of <code>String</code>'s.
   * @return
   */
  public List<String> getFiles(String patternStr)
  {
    if(!statusOK)
      return new Vector<String>();
    Pattern p = Pattern.compile(patternStr);
    
    List<String> files = new Vector<String>();
    for(int i=0; i<bamFields.size(); i++)
    {
      String file = bamFields.get(i).getText();
      if(file != null && !file.equals(""))
      {
        if(isListOfFiles(file))
        {
          List<String> filesInList = getListOfFiles(file);
          for(int j=0; j<filesInList.size(); j++)
          {
            Matcher m = p.matcher(filesInList.get(j));
            if(m.matches())
              files.add(filesInList.get(j));
          }
        }
        else
        {
          Matcher m = p.matcher(file);
          if(m.matches())
            files.add(file);
        }
      }
    }
    return files;
  }
  
  /**
   * Get the reference
   * @return
   */
  public String getReferenceFile()
  {
    return referenceField.getText();
  }
      
}