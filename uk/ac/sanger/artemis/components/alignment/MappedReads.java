package uk.ac.sanger.artemis.components.alignment;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.Toolkit;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;

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
  private HashMap<String, Integer> offsetLengths;
  private boolean concatSequences;
  private HashMap<String, Integer> seqLengths;
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
      final HashMap<String, Integer> offsetLengths,
      final boolean concatSequences,
      final HashMap<String, Integer> seqLengths, 
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
      final HashMap<String, Integer> offsetLengths,
      final boolean concatSequences,
      final HashMap<String, Integer> seqLengths,
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
  
  

  /**
   * Search for new features based on a threshold of read counts in the intergenic 
   * and anti-sense regions of existing annotations. This should exclude rRNA and 
   * tRNA regions. If two or more BAM files are loaded it should create features 
   * based on the combined read peak span.
   * @param feature_display
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
   */
  public MappedReads(
      final FeatureDisplay feature_display,
      final String refName,
      final Hashtable<String, SAMFileReader> samFileReaderHash,
      final List<String> bamList, 
      final Vector<String> seqNames,
      final HashMap<String, Integer> offsetLengths,
      final boolean concatSequences,
      final HashMap<String, Integer> seqLengths,
      final SAMRecordPredicate samRecordFlagPredicate,
      final SAMRecordMapQPredicate samRecordMapQPredicate,
      final int threshold,
      final int minSize,
      final boolean contained)
  {
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
    
    progressBar = new JProgressBar(0, feature_display.getSequenceLength());
    progressBar.setValue(0);
    progressBar.setStringPainted(true);

    JPanel panel = new JPanel(new BorderLayout());
    progressTxt.setText("");
    panel.add(progressTxt, BorderLayout.NORTH);
    panel.add(progressBar, BorderLayout.CENTER);

    panel.setOpaque(true);
    dialog.setTitle("Search");
    dialog.setContentPane(panel);
    dialog.pack();
    centerDialog();
    
    CalculateNewFeatures cmr = new CalculateNewFeatures(feature_display, refName, threshold, minSize);
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
        featureReadCount.put(ReadCountDialog.getFeatureName(f), sampleCounts);
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
    
    private void calc(final String refName, final int sequenceLength)
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
  
  /**
   * Find new features from read count peaks.
   */
  class CalculateNewFeatures extends SwingWorker
  {
    private EntryGroup entryGroup;
    private Bases bases;
    private String refSeq;
    private int threshold;
    private int minSize;
    
    CalculateNewFeatures(final FeatureDisplay feature_display, 
        final String refSeq, final int threshold,
        final int minSize)
    {
      entryGroup = feature_display.getEntryGroup();
      bases = feature_display.getBases();
      this.refSeq = refSeq;
      this.threshold = threshold;
      this.minSize = minSize;
    }
    
    public Object construct()
    {
      int MAX_BASE_CHUNK = 2000 * 60;
      Key excluseKeys[] = { new Key("rRNA"), new Key("tRNA") };      
      int seqlen = bases.getLength();

      final int beg = 1;
      final int end = seqlen;

      int fwdStart = -1;
      int revStart = -1;
      final List<MarkerRange> fwdMarkers = new Vector<MarkerRange>();
      final List<MarkerRange> revMarkers = new Vector<MarkerRange>();
      for (int i = 0; i < bamList.size(); i++)
      {
        for(int j=beg; j<end; j+=MAX_BASE_CHUNK)
        {
          progressBar.setValue((j + (i*end)) / bamList.size());
          if(j > end)
            continue;
          int start = j;
          int stop  = j+MAX_BASE_CHUNK;

          try
          {
            int cnt[][] = new int[stop-start+1][2];
            for (int row = 0; row < cnt.length; row++)
              for (int col = 0; col < 2; col++)
                cnt[row][col] = 0;
            
            if (concatSequences)
            {
              for (String name : seqNames)
              {
                int len = seqLengths.get(name);
                int offset = offsetLengths.get(name);
                
                if( (start >= offset && start <= offset+len) ||
                    (stop >= offset  && start <= offset+len) )
                {
                  cnt =
                    BamUtils.countOverRange(bamList.get(i), samFileReaderHash, 
                        name, start-offset, stop-offset, cnt,
                        samRecordFlagPredicate, samRecordMapQPredicate);
                }
              }
            }
            else
            {
              cnt =
                BamUtils.countOverRange(bamList.get(i), samFileReaderHash, 
                    refSeq, start, stop, cnt,
                    samRecordFlagPredicate, samRecordMapQPredicate);
            }
            
            for(int k=0; k<cnt.length; k++)
            {
              final Range r = new Range(start+k, start+k+1);
              
              // find forward strand potential features
              fwdStart = findFeatures(cnt[k][0], true, fwdStart, j+k,
                r, excluseKeys, fwdMarkers, entryGroup);
            
              // find reverse strand potential features
              revStart = findFeatures(cnt[k][1], false, revStart, j+k,
                r, excluseKeys, revMarkers, entryGroup);
            }
            
          }
          catch (OutOfRangeException e1)
          {
            e1.printStackTrace();
          }
        }
      }

      final Entry newEntry = entryGroup.createEntry ("align_"+threshold+"_"+minSize);
      createFeatures(fwdMarkers, true, newEntry);
      createFeatures(revMarkers, false, newEntry);
      return null;
    }
    
    
    private void createFeatures(final List<MarkerRange> markers, final boolean isFwd, final Entry newEntry)
    {
      final Key key = Key.CDS;
/*      if(entryGroup.getDefaultEntry() != null &&
         entryGroup.getDefaultEntry().getEMBLEntry() instanceof GFFDocumentEntry)
        key = new Key("region");
      else
        key = new Key("misc_feature");*/
      
      // merge overlapping markers
      List<MarkerRange> featureMarkers = new Vector<MarkerRange>();
      List<Integer> ignoreMarkers = new Vector<Integer>();
      for(int i=0; i<markers.size(); i++)
      {
        if(ignoreMarkers.contains(i))
          continue;
        MarkerRange mi = markers.get(i);
        for(int j=i+1; j<markers.size(); j++)
        {
          if(ignoreMarkers.contains(j))
            continue;
          MarkerRange mj = markers.get(j);
          if(mi.overlaps(mj))
          {
            mi = mi.combineRanges(mj, false);
            ignoreMarkers.add(j);
          }
        }
        featureMarkers.add(mi);
      }
      
      for(MarkerRange r: featureMarkers)
      {
        try
        {
          newEntry.createFeature(key, 
              ( isFwd ? r.createLocation() : r.createLocation().getComplement()), null);
        }
        catch (ReadOnlyException e1)
        {
          e1.printStackTrace();
        }
        catch (EntryInformationException e1)
        {
          e1.printStackTrace();
        }
        catch (OutOfRangeException e1)
        {
          e1.printStackTrace();
        }
      }
    }
    
    /**
     * Search for new features based on read counts.
     * @param cnt             read count
     * @param isFwd           strand to find features on is the forward strand if true
     * @param fStart          start of a new feature or -1 if not identified yet
     * @param pos             current base position
     * @param range           current base range
     * @param excluseKeys     feature keys of regions to avoid looking in
     * @param markers         list of new features
     * @param entryGroup      entry group
     * @return
     */
    private int findFeatures(final int cnt,
                             final boolean isFwd,
                             int fStart, 
                             final int pos,
                             final Range range,  
                             final Key excluseKeys[],
                             final List<MarkerRange> markers,
                             final EntryGroup entryGroup)
    {
      if(cnt > threshold && fStart == -1)  // START FEATURE
      {
        boolean exclude = false;
        try
        {
          final FeatureVector features =
              entryGroup.getFeaturesInRange(range);
          for(int k=0; k<features.size(); k++)
          {
            Feature f = features.elementAt(k);
            if( f.isProteinFeature() && ! (isFwd ^ f.isForwardFeature()) ) // on same strand
              return fStart;
            
            for(int l=0; l<excluseKeys.length; l++)
              if(f.getKey().equals(excluseKeys[l]))
                exclude = true;
          }
        }
        catch (OutOfRangeException e1){ }

        if(!exclude)
          fStart = range.getStart();
      }
      else if(cnt < threshold && fStart != -1)
      {
        try
        {
          boolean exclude = false;
          final FeatureVector features = entryGroup.getFeaturesInRange(range);
          for(int k=0; k<features.size(); k++)
          {
            Feature f = features.elementAt(k);
            for(int l=0; l<excluseKeys.length; l++)
              if(f.getKey().equals(excluseKeys[l]))
                exclude = true;
          }
          
          if(!exclude &&
             (range.getStart()-fStart) >= minSize)
          {
            final MarkerRange marker = new MarkerRange(
              bases.getForwardStrand(), fStart, range.getStart());
            markers.add( marker );
          }
        }
        catch (OutOfRangeException e1){}
        
        return -1;
      }
      else if (fStart != -1)
      {
        try
        {
          FeatureVector features = entryGroup.getFeaturesInRange(range);
          for(int k=0; k<features.size(); k++)
          {
            boolean exclude = false;
            Feature f = features.elementAt(k);
            for(int l=0; l<excluseKeys.length; l++)
              if(f.getKey().equals(excluseKeys[l]))
                exclude = true;
            
            if( exclude ||
                (f.isProteinFeature() && ! (isFwd ^ f.isForwardFeature())) ) // on same strand
            {
              if((f.getRawFirstBase()-fStart) >= minSize)
              {
                final MarkerRange marker = new MarkerRange(
                  bases.getForwardStrand(), fStart, f.getRawFirstBase() );
                markers.add( marker );
              }
              return -1;
            }
          }
        }
        catch (OutOfRangeException e)
        {
          e.printStackTrace();
        }
        
      }
      return fStart;
    }
    
    
    public void finished() 
    {
      dialog.dispose();
    }
  }
}
