package uk.ac.sanger.artemis.components.alignment;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.text.DecimalFormat;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Set;
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

public class MappedReads
{
  private JProgressBar progressBar;
  private JLabel progressTxt = new JLabel();

  private FeatureVector features;
  private BamView bamView;

  private boolean contained;
  private boolean useIntrons;
  private boolean useStrandTag = false;
  private JDialog dialog = new JDialog((JFrame)null, "Calculating", true);
    
  /**
   * Calculate the total number of mapped reads.
   * @param bamView
   * @param features
   * @param sequenceLength
   * @param contained
   * @param useIntrons
   * @param allRefSeqs
   * @param useStrandTag
   */
  public MappedReads(
      final BamView bamView,
      final FeatureVector features,
      final int sequenceLength,
      final boolean contained, 
      final boolean useIntrons,
      final boolean allRefSeqs,
      final boolean useStrandTag)
  {
    this.features = features;
    this.bamView = bamView;
    this.contained = contained;
    this.useIntrons = useIntrons;
    this.useStrandTag = useStrandTag;

    showProgress(sequenceLength, "Total number of mapped reads");
    CalculateTotalMappedReads cmr = new CalculateTotalMappedReads(allRefSeqs, sequenceLength);
    cmr.start();
    dialog.setVisible(true);
  }

  /**
   * Read count for selected features.
   * @param bamView
   * @param features
   * @param contained
   * @param useIntrons
   * @param useStrandTag
   */
  public MappedReads(
      final BamView bamView,
      final FeatureVector features, 
      final boolean contained, 
      final boolean useIntrons,
      final boolean useStrandTag)
  {
    this.bamView = bamView;
    this.features = features;
    this.bamView = bamView;
    this.contained = contained;
    this.useIntrons = useIntrons;
    this.useStrandTag = useStrandTag;

    showProgress(features.size(), "Number of mapped reads for "+features.size()+" features");
    CalculateMappedReads cmr = new CalculateMappedReads(null);
    cmr.start();
    dialog.setVisible(true);
  }

  /**
   * Search for new features based on a threshold of read counts in the intergenic 
   * and anti-sense regions of existing annotations. This should exclude rRNA and 
   * tRNA regions. If two or more BAM files are loaded it should create features 
   * based on the combined read peak span.
   * @param bamView
   * @param groupsFrame
   * @param threshold
   * @param minSize
   * @param minBams
   * @param readsOnOppositeStrand
   * @param contained
   */
  public MappedReads(  
      final BamView bamView,
      final GroupBamFrame groupsFrame,
      final int threshold,
      final int minSize,
      final int minBams,
      final boolean readsOnOppositeStrand,
      final boolean contained)
  {
    this.bamView = bamView;
    this.contained = contained;
    
    final FeatureDisplay feature_display = bamView.getFeatureDisplay();
    showProgress(feature_display.getSequenceLength(), "");
    final CalculateNewFeatures cmr = new CalculateNewFeatures(
        feature_display, (String) bamView.getCombo().getSelectedItem(), groupsFrame, 
        threshold, minSize, minBams, readsOnOppositeStrand);
    cmr.start();
    dialog.setVisible(true);
  }
  
  private void showProgress(int progressBarMax, String txt)
  {
    progressBar = new JProgressBar(0, progressBarMax);
    progressBar.setValue(0);
    progressBar.setStringPainted(true);

    JPanel panel = new JPanel(new BorderLayout());
    progressTxt.setText(txt);
    panel.add(progressTxt, BorderLayout.NORTH);
    panel.add(progressBar, BorderLayout.CENTER);

    panel.setOpaque(true);
    dialog.setTitle("Search");
    dialog.setContentPane(panel);
    dialog.pack();
    centerDialog();
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
    private Hashtable<String, List<ReadCount>> featureReadCount;
    private int mappedReads[];
    
    public CalculateMappedReads(final int mappedReads[])
    {
      this.mappedReads = mappedReads;
    }
    
    public Object construct()
    {
      featureReadCount = BamUtils.calculateMappedReads(bamView, features, 
          contained, useIntrons, useStrandTag, mappedReads, progressBar);
      return null;
    }
       
    public void finished() 
    {
      final DecimalFormat df;
      if(mappedReads != null)
        df = new DecimalFormat("0.000");
      else
        df = new DecimalFormat("0");

      // use two possible titles used for the truncated names and
      // the full names of BAM's when saved to file
      final StringBuilder hdr = new StringBuilder();
      final StringBuilder title = new StringBuilder();
      final StringBuilder titleToSave = new StringBuilder();
      final StringBuilder body = new StringBuilder();
      
      for (int j = 0; j < bamView.bamList.size(); j++)
      {
        String bam = bamView.bamList.get(j);
        hdr.append("#BAM: " + bam);
        if (mappedReads != null)
        {
          final DecimalFormat df2 = new DecimalFormat("0.000000");
          hdr.append(" Mapped Reads/million: "
              + df2.format(((float) mappedReads[j]) / 1000000.f));
        }
        hdr.append("\n");
      }
      hdr.append("\n");

      final Enumeration<String> enumFeat = featureReadCount.keys();
      int n = 0;
      while (enumFeat.hasMoreElements())
      {
        final String fId = enumFeat.nextElement();
        if(n == 0)
        {
          for (int i = 0; i < fId.length(); i++)
          {
            titleToSave.append(" ");
            title.append(" ");
          }
          titleToSave.append("\t");
          title.append("\t");
          for (String bam: bamView.bamList)
          {
            String name = new File(bam).getName();
            titleToSave.append(name+"\t");
            if(mappedReads == null)
            {
              if(name.length() > 21)
                 name = "  "+name.substring(0, 19) + "...";
              else
              {
                int len = 21-name.length();
                name = "  "+name;
                for (int i = 0; i < len; i++)
                  name = name+" ";
              }
            }
            else if(mappedReads != null)
            {
              if(name.length() > 28)
                 name = "      "+name.substring(0, 26) + "...";
              else
              {
                int len = 28-name.length(); 
                name = "      "+name;
                for (int i = 0; i < len; i++)
                  name = name+" ";
              }
            }
            title.append(name+"\t");
          }
          
          title.append("\n");
          titleToSave.append("\n");
          for (int i = 0; i < fId.length(); i++)
          {
            title.append(" ");
            titleToSave.append(" ");
          }
          title.append("\t");
          titleToSave.append("\t");
          for (int j = 0; j < bamView.bamList.size(); j++)
          {
            if(mappedReads != null)
            {
              title.append("      Sense   Antisense       Total\t");
              titleToSave.append("      Sense   Antisense       Total\t");
            }
            else
            {
              title.append("  Sense Antisense  Total\t");
              titleToSave.append("  Sense Antisense  Total\t");
            }
          }
          title.append("\n");
          titleToSave.append("\n");
        }
        
        body.append(fId + "\t");
        final List<ReadCount> cnts = featureReadCount.get(fId);
        
        for (ReadCount c: cnts)
        {
          pad(body, c.senseCnt, df);
          body.append(" ");
          pad(body, c.antiCnt, df);
          body.append(" ");
          pad(body, c.senseCnt+c.antiCnt, df);
          body.append("\t");
        }
        body.append("\n");
        n++;
      }

      final FileViewer viewer;
      final ActionListener saveListener = new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          FileViewer.writeToFile(hdr.toString()+
              titleToSave.toString()+body.toString());
        }
      };
      if (mappedReads != null)
        viewer = new FileViewer("RPKM", true, false, true, saveListener);
      else
        viewer = new FileViewer("Read Count", true, false, true, saveListener);
      viewer.getTextPane().setText(hdr.toString()+title.toString()+body.toString());
      
      dialog.dispose();
    }
  }
  
  private static void pad(StringBuilder buff, float f, final DecimalFormat df)
  {
    if(f < 10)
      buff.append("      ");
    else if(f < 100)
      buff.append("     ");
    else if(f < 1000)
      buff.append("    ");
    else if(f < 10000)
      buff.append("   ");
    else if(f < 100000)
      buff.append("  ");
    else if(f < 1000000)
      buff.append(" ");
    buff.append(df.format(f));
  }
  
  class CalculateTotalMappedReads extends SwingWorker
  {
    private boolean useAllRefSeqs;
    private int sequenceLength;
    
    CalculateTotalMappedReads(boolean useAllRefSeqs, final int sequenceLength)
    {
      this.useAllRefSeqs = useAllRefSeqs;
      this.sequenceLength = sequenceLength;
    }
    
    public Object construct()
    {
      int mappedReads[] = null;
      if(bamView.isConcatSequences() || !useAllRefSeqs)
      {
        String refName = (String) bamView.getCombo().getSelectedItem();
        mappedReads = BamUtils.calc(bamView, refName, sequenceLength, 
            useStrandTag, progressBar);
      }
      else
      {
        for (String name : bamView.getSeqNames())
        {
          progressTxt.setText(name);
          int thisLength = bamView.getSeqLengths().get(name);
          progressBar.setValue(0);
          progressBar.setMaximum(thisLength);
          mappedReads = BamUtils.calc(bamView, name, thisLength, 
              useStrandTag, progressBar);
        }
      }
      
      progressBar.setValue(0);
      progressBar.setMaximum(features.size());
      progressTxt.setText("RPKM values for "+features.size()+" features");
      CalculateMappedReads cmr = new CalculateMappedReads(mappedReads);
      cmr.start();
      return null;
    }
  }
  
  /**
   * Find new features from read count peaks.
   */
  class CalculateNewFeatures extends SwingWorker
  {
    private EntryGroup entryGroup;
    private Bases bases;
    private int threshold;
    private int minSize;
    private int minBams;
    private boolean readsOnOppositeStrand;
    private GroupBamFrame groupsFrame;
    
    CalculateNewFeatures(final FeatureDisplay feature_display,
        final String refSeq, 
        final GroupBamFrame groupsFrame, 
        final int threshold,
        final int minSize,
        final int minBams,
        final boolean readsOnOppositeStrand)
    {
      entryGroup = feature_display.getEntryGroup();
      bases = feature_display.getBases();

      this.groupsFrame = groupsFrame;
      this.threshold = threshold;
      this.minSize = minSize;
      this.minBams = minBams;
      this.readsOnOppositeStrand = readsOnOppositeStrand;
    }
    
    class MarkerObj
    {
      MarkerRange r;
      int bamIdx;
      MarkerObj(MarkerRange r, int bamIdx)
      {
        this.r = r;
        this.bamIdx = bamIdx;
      }
    }
    
    public Object construct()
    {
      final int MAX_BASE_CHUNK = 2000 * 60;
      Key excluseKeys[] = { new Key("rRNA"), new Key("tRNA") };      

      final int beg = 1;
      final int end = bases.getLength();

      int fwdStart = -1;
      int revStart = -1;
      final List<MarkerObj> fwdMarkers = new Vector<MarkerObj>();
      final List<MarkerObj> revMarkers = new Vector<MarkerObj>();
      for (short i = 0; i < bamView.bamList.size(); i++)
      {
        if(bamView.hideBamList.contains(i))
          continue;
        
        for(int j=beg; j<end; j+=MAX_BASE_CHUNK)
        {
          progressBar.setValue((j + (i*end)) / bamView.bamList.size());
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
            
            if (bamView.isConcatSequences())
            {
              for (String name : bamView.getSeqNames())
              {
                int len = bamView.getSeqLengths().get(name);
                int offset = bamView.getOffsetLengths().get(name);
                
                if( (start >= offset && start <= offset+len) ||
                    (stop >= offset  && start <= offset+len) )
                {
                  int thisStart = start-offset;
                  if(thisStart < 1)
                    thisStart = 1;
                  int thisEnd = stop-offset;
                  if(thisEnd > len)
                    thisEnd = len;
                  
                  int concatShift = 0;
                  if(offset > start)
                    concatShift = offset-start;

                  cnt =
                    BamUtils.countOverRange(bamView, bamView.bamList.get(i), 
                        thisStart, thisEnd, concatShift, cnt);
                }
              }
            }
            else
              cnt =
                BamUtils.countOverRange(bamView, bamView.bamList.get(i),  
                    start, stop, 0, cnt);

            for(int k=0; k<cnt.length; k++)
            {
              final Range r = new Range(start+k, start+k+1);
              // find forward strand potential features
              fwdStart = findFeatures( 
                  (!readsOnOppositeStrand ? cnt[k][0] : cnt[k][1]), true, 
                  fwdStart, j+k, r, excluseKeys, fwdMarkers, entryGroup, i);

              // find reverse strand potential features
              revStart = findFeatures( 
                  (!readsOnOppositeStrand ? cnt[k][1] : cnt[k][0]), false, 
                  revStart, j+k, r, excluseKeys, revMarkers, entryGroup, i);
            }
          }
          catch (OutOfRangeException e1)
          {
            e1.printStackTrace();
          }
        }
      }

      final Entry newEntry = entryGroup.createEntry (
          "align_"+threshold+"_"+minBams+"_"+minSize+(!readsOnOppositeStrand ? "":"_opp"));
      createFeatures(fwdMarkers, true, newEntry);
      createFeatures(revMarkers, false, newEntry);
      return null;
    }
    
    
    private void createFeatures(final List<MarkerObj> markers, final boolean isFwd, final Entry newEntry)
    {
      final Key key = Key.CDS;
      
      // merge overlapping markers
      final Set<Integer> ignoreMarkers = new HashSet<Integer>();
      for(int i=0; i<markers.size(); i++)
      {
        if(ignoreMarkers.contains(i))
          continue;
        
        Set<Integer> bamIdxList = new HashSet<Integer>();
        MarkerObj mi = markers.get(i);
        MarkerRange mri = mi.r;
        bamIdxList.add(mi.bamIdx);
        for(int j=i+1; j<markers.size(); j++)
        {
          if(ignoreMarkers.contains(j))
            continue;
          MarkerObj mj = markers.get(j);
          MarkerRange mrj = mj.r;

          if(mri.overlaps(mrj))
          {
            bamIdxList.add(mj.bamIdx);
            mri = mri.combineRanges(mrj, false);
            ignoreMarkers.add(j);
          }
        }
        
        if(bamIdxList.size() >= minBams)
        {
          if(groupsFrame != null && minBams > 1)
          {
            boolean foundMinBams = false;
            Hashtable<String, Integer> groupCluster = new Hashtable<String, Integer>();
            for(int j=0; j<bamIdxList.size(); j++)
            {
              File f = new File(bamView.bamList.get(j));
              String grp = groupsFrame.getGroupName( f.getName() );
              if(groupCluster.containsKey(grp))
              {
                int val = groupCluster.get(grp)+1;
                if(val >= minBams)
                {
                  foundMinBams = true;
                  break;
                }
                groupCluster.put(grp, val);
              }
              else
                groupCluster.put(grp, 1);
            }
            if(!foundMinBams)
              continue;
          }
          
          // create a new feature
          try
          {
            newEntry.createFeature(key, 
                ( isFwd ? mri.createLocation() : mri.createLocation().getComplement()), null);
          }
          catch (ReadOnlyException e)
          {
            e.printStackTrace();
          }
          catch (EntryInformationException e)
          {
            e.printStackTrace();
          }
          catch (OutOfRangeException e)
          {
            e.printStackTrace();
          }
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
                             final List<MarkerObj> markers,
                             final EntryGroup entryGroup,
                             final int bamIdx)
    {
      if(cnt >= threshold && fStart == -1)  // START FEATURE
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
            markers.add( new MarkerObj(marker, bamIdx) );
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
                markers.add( new MarkerObj(marker, bamIdx) );
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
