package uk.ac.sanger.artemis.components.alignment;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.Toolkit;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.SwingWorker;

import net.sf.samtools.SAMFileReader;

public class MappedReads
{
  private JProgressBar progressBar;
  private JLabel progressTxt = new JLabel();

  private FeatureVector features;
  private String refName;
  private Hashtable<String, SAMFileReader> samFileReaderHash;
  private List<String> bamList;
  private Vector<String> seqNames;
  private Hashtable<String, Integer> offsetLengths;
  private boolean concatSequences;
  private Hashtable<String, Integer> seqLengths;
  private int sequenceLength;
  private SAMRecordPredicate samRecordFlagPredicate;
  private SAMRecordMapQPredicate samRecordMapQPredicate;
  private boolean contained;
  private boolean useIntrons;
  private JDialog dialog = new JDialog((JFrame)null, "Calculating", true);;
    
  private int mappedReads[];
    
  /**
   * Calculate the total number of mapped reads.
   * @param refName
   * @param samFileReaderHash
   * @param bamList
   * @param seqNames
   * @param offsetLengths
   * @param concatSequences
   * @param seqLengths
   * @param sequenceLength
   */
  public MappedReads(
      final FeatureVector features,
      final String refName,
      final Hashtable<String, SAMFileReader> samFileReaderHash,
      final List<String> bamList, 
      final Vector<String> seqNames,
      final Hashtable<String, Integer> offsetLengths,
      final boolean concatSequences,
      final Hashtable<String, Integer> seqLengths, 
      final int sequenceLength,
      final SAMRecordPredicate samRecordFlagPredicate,
      SAMRecordMapQPredicate samRecordMapQPredicate,
      final boolean contained, 
      final boolean useIntrons,
      final boolean allRefSeqs)
  {
    this.features = features;
    this.refName = refName;
    this.samFileReaderHash = samFileReaderHash;
    this.bamList = bamList;
    this.seqNames = seqNames;
    this.offsetLengths = offsetLengths;
    this.concatSequences = concatSequences;
    this.seqLengths = seqLengths;
    this.sequenceLength = sequenceLength;
    this.samRecordFlagPredicate = samRecordFlagPredicate;
    this.samRecordMapQPredicate = samRecordMapQPredicate;
    this.contained = contained;
    this.useIntrons = useIntrons;
    
    progressBar = new JProgressBar(0, sequenceLength);
    progressBar.setValue(0);
    progressBar.setStringPainted(true);

    JPanel panel = new JPanel(new BorderLayout());
    progressTxt.setText("Total number of mapped reads");
    panel.add(progressTxt, BorderLayout.NORTH);
    panel.add(progressBar, BorderLayout.CENTER);

    dialog.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    panel.setOpaque(true);
    dialog.setContentPane(panel);
    dialog.pack();
    centerDialog();

    CalculateTotalMappedReads cmr = new CalculateTotalMappedReads(allRefSeqs);
    cmr.start();
    dialog.setVisible(true);
  }

  /**
   * Read count for selected features.
   * @param features
   * @param refName
   * @param samFileReaderHash
   * @param bamList
   * @param seqNames
   * @param offsetLengths
   * @param concatSequences
   * @param seqLengths
   * @param samRecordFlagPredicate
   * @param samRecordMapQPredicate
   * @param contained
   * @param useIntrons
   */
  public MappedReads(
      final FeatureVector features, 
      final String refName,
      final Hashtable<String, SAMFileReader> samFileReaderHash,
      final List<String> bamList, 
      final Vector<String> seqNames,
      final Hashtable<String, Integer> offsetLengths,
      final boolean concatSequences,
      final Hashtable<String, Integer> seqLengths,
      final SAMRecordPredicate samRecordFlagPredicate,
      final SAMRecordMapQPredicate samRecordMapQPredicate,
      final boolean contained, 
      final boolean useIntrons)
  {
    this.features = features;
    this.refName = refName;
    this.samFileReaderHash = samFileReaderHash;
    this.bamList = bamList;
    this.seqNames = seqNames;
    this.offsetLengths = offsetLengths;
    this.concatSequences = concatSequences;
    this.seqLengths = seqLengths;
    this.samRecordFlagPredicate = samRecordFlagPredicate;
    this.samRecordMapQPredicate = samRecordMapQPredicate;
    this.contained = contained;
    this.useIntrons = useIntrons;
    
    progressBar = new JProgressBar(0, features.size());
    progressBar.setValue(0);
    progressBar.setStringPainted(true);

    JPanel panel = new JPanel(new BorderLayout());
    progressTxt.setText("Number of mapped reads for "+features.size()+" features");
    panel.add(progressTxt, BorderLayout.NORTH);
    panel.add(progressBar, BorderLayout.CENTER);

    panel.setOpaque(true);
    dialog.setContentPane(panel);
    dialog.pack();
    centerDialog();
    
    CalculateMappedReads cmr = new CalculateMappedReads();
    cmr.start();
    dialog.setVisible(true);
  }

  private void centerDialog()
  {
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    final int x_position =(screen.width - dialog.getSize().width) / 2;
    int y_position =(screen.height - dialog.getSize().height) / 2;

    if(y_position < 10) 
      y_position = 10;
    dialog.setLocation(new Point(x_position, y_position));  
  }
  
  class CalculateMappedReads extends SwingWorker
  {
    Hashtable<String, List<Float>> featureReadCount;
    public Object construct()
    {
      featureReadCount = new Hashtable<String, List<Float>>();
      for (int i = 0; i < features.size(); i++)
      {
        Feature f = features.elementAt(i);
        progressBar.setValue(i);

        int start = f.getRawFirstBase();
        int end = f.getRawLastBase();
        float fLen = BamUtils.getFeatureLength(f);
        List<Float> sampleCounts = new Vector<Float>();

        for (int j = 0; j < bamList.size(); j++)
        {
          String bam = bamList.get(j);
          float cnt = 0;
          if (!useIntrons && f.getSegments().size() > 1)
          {
            for (int k = 0; k < f.getSegments().size(); k++)
            {
              start = f.getSegments().elementAt(k).getRawRange().getStart();
              end = f.getSegments().elementAt(k).getRawRange().getEnd();
              cnt += BamUtils.getCount(start, end, bam, refName, samFileReaderHash,
                  seqNames, offsetLengths, concatSequences, seqLengths,
                  samRecordFlagPredicate, samRecordMapQPredicate, contained);
            }
          }
          else
            cnt = BamUtils.getCount(start, end, bam, refName, samFileReaderHash, seqNames,
                offsetLengths, concatSequences, seqLengths,
                samRecordFlagPredicate, samRecordMapQPredicate, contained);

          if (mappedReads != null)
            cnt = (cnt / (((float) mappedReads[j] / 1000000.f) * (fLen / 1000.f)));

          sampleCounts.add(cnt);
        }
        featureReadCount.put(f.getSystematicName(), sampleCounts);
      }
      return null;
    }
    
    public void finished() 
    {
      DecimalFormat df = new DecimalFormat("0.00##");
      StringBuffer buff = new StringBuffer();
      for (int j = 0; j < bamList.size(); j++)
      {
        String bam = bamList.get(j);
        buff.append("#BAM: " + bam);
        if (mappedReads != null)
          buff.append(" Mapped Reads/million: "
              + df.format(((float) mappedReads[j]) / 1000000.f));
        buff.append("\n");
      }
      buff.append("\n");

      Object[] readKey = featureReadCount.keySet().toArray();
      Arrays.sort(readKey);

      for (Object fId : readKey)
      {
        buff.append(fId + "\t");
        List<Float> cnts = featureReadCount.get(fId);
        for (int i = 0; i < cnts.size(); i++)
          buff.append(df.format(cnts.get(i)) + (i < cnts.size() - 1 ? "\t" : ""));
        buff.append("\n");
      }

      FileViewer viewer;
      if (mappedReads != null)
        viewer = new FileViewer("RPKM", true, false, true);
      else
        viewer = new FileViewer("Read Count", true, false, true);
      viewer.getTextPane().setText(buff.toString());
      
      dialog.dispose();
    }
  }
  
  class CalculateTotalMappedReads extends SwingWorker
  {
    private boolean useAllRefSeqs;
    CalculateTotalMappedReads(boolean useAllRefSeqs)
    {
      this.useAllRefSeqs = useAllRefSeqs;
    }
    
    public Object construct()
    {
      mappedReads = new int[bamList.size()];
      if(concatSequences || !useAllRefSeqs)
        calc(refName, sequenceLength);
      else
      {
        for (String name : seqNames)
        {
          progressTxt.setText(name);
          int thisLength = seqLengths.get(name);
          progressBar.setValue(0);
          progressBar.setMaximum(thisLength);
          calc(name, thisLength);
        }
      }
      
      progressBar.setValue(0);
      progressBar.setMaximum(features.size());
      progressTxt.setText("RPKM values for "+features.size()+" features");
      CalculateMappedReads cmr = new CalculateMappedReads();
      cmr.start();
      return null;
    }
    
    private void calc(String refName, int sequenceLength)
    {
      int MAX_BASE_CHUNK = 2000 * 60;
      boolean contained = false;
      for (int i = 0; i < sequenceLength; i += MAX_BASE_CHUNK)
      {
        progressBar.setValue(i);
        int sbegc = i;
        int sendc = i + MAX_BASE_CHUNK - 1;

        for (int j = 0; j < bamList.size(); j++)
        {
          String bam = bamList.get(j);
          if (concatSequences)
          {
            int len = 0;
            int lastLen = 1;
            for (String name : seqNames)
            {
              int thisLength = seqLengths.get(name);
              len += thisLength;

              if ((lastLen >= sbegc && lastLen < sendc)
                  || (len >= sbegc && len < sendc)
                  || (sbegc >= lastLen && sbegc < len)
                  || (sendc >= lastLen && sendc < len))
              {
                int offset = offsetLengths.get(name);
                int thisStart = sbegc - offset;
                if (thisStart < 1)
                  thisStart = 1;
                int thisEnd = sendc - offset;
                if (thisEnd > thisLength)
                  thisEnd = thisLength;

                mappedReads[j] += BamUtils.count(bam, samFileReaderHash, name,
                    thisStart, thisEnd, samRecordFlagPredicate,
                    samRecordMapQPredicate, contained);

              }
              lastLen = len;
            }
          }
          else
          {
            mappedReads[j] += BamUtils.count(bam, samFileReaderHash, refName, sbegc,
                sendc, samRecordFlagPredicate, samRecordMapQPredicate,
                contained);
          }
        }
      }
    }
  }
}
