package uk.ac.sanger.artemis.io;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Hides the implementation of the line groups list.
 * If we implement multi-sequence Genbank loading then 
 * this may need to change a bit later on.
 * 
 * @author kp11
 *
 */
public class LineGroupCache
{
	@SuppressWarnings("unused")
	private static final long serialVersionUID = 1L;
	
	public static final String FORWARD_REF_DELIM = "###";
	
	private LinkedHashMap<Long, LineGroup> headers = new LinkedHashMap<Long, LineGroup>();
	
	private LinkedHashMap<Long, LineGroup> misc = new LinkedHashMap<Long, LineGroup>();
	
	/**
	 * The Sequence line group.
	 */
	private LineGroup sequence;
	
	/**
	 * The FeatureTable line group.
	 */
	private FeatureTable featureTable;
	
	/**
	 * Determine which internal data structure is associated with the given group.
	 * @param group LineGroup
	 * @return LinkedHashMap
	 */
	protected LinkedHashMap<Long, LineGroup> findRelevantList(LineGroup group)
	{
		if (group instanceof FeatureHeader)
		{
			return headers;
		}
		
		return misc;
	}
	
	/**
	 * Return the feature table line group.
	 * 
	 * @return FeatureTable
	 */
	public FeatureTable getFeatureTable()
	{
		return featureTable;
	}
	
	/**
	 * Return the sequence line group.
	 * 
	 * @return LineGroup
	 */
	public LineGroup getSequence()
	{
		return sequence;
	}
	
	/**
	 * Add a line group to the end of the line group list.
	 * 
	 * @param lineGroup LineGroup
	 */
	public synchronized void add(LineGroup lineGroup)
	{
		if (lineGroup instanceof FeatureTable)
		{
			featureTable = (FeatureTable)lineGroup;
		}
		else if (lineGroup instanceof Sequence)
		{
			sequence = lineGroup;
		}
		else if (lineGroup instanceof GFFMisc && FORWARD_REF_DELIM.equals(((GFFMisc) lineGroup).getLine()))
		{
			// Remove these forward ref delimiters as we don't need them.
			return;
		}
		else
		{
			findRelevantList(lineGroup).put(lineGroup.getUniqueId(), lineGroup);
		}
	}

	/**
	 * Remove the given line group from this cache and return it.
	 * @param lineGroup LineGroup
	 * @return LineGroup
	 */
	public synchronized LineGroup remove(LineGroup lineGroup)
	{	
		LineGroup group = null;
		
		if (lineGroup instanceof FeatureTable)
		{
			group = featureTable;
			featureTable = null;
			
		}
		else if (lineGroup instanceof Sequence)
		{
			group = sequence;
			sequence = null;
		}
		else 
		{
			group = findRelevantList(lineGroup).remove(lineGroup.getUniqueId());
		}
		
		return group;
	}
    
    /**
     *  Create and return a new FeatureTable().  The new FeatureTable will be
     *  added to line_groups in the appropriate place.
     *  @return FeatureTable
     **/
    public synchronized FeatureTable createFeatureTable() 
    {
    	featureTable = new StreamFeatureTable();

    	return featureTable;
    }
    
    /**
     *  Return the text of the "header" of this Entry or null if there is no
     *  header.
     *  @return String
     **/
    public synchronized String getHeadersAsText() 
    {
      int numLineGroups = headers.size() + misc.size();
      
      List<LineGroup> groups = new ArrayList<LineGroup>(numLineGroups);
      
      groups.addAll(misc.values());
      groups.addAll(headers.values());
      
      final StringBuilder buffer = 
    		  new StringBuilder( numLineGroups*80 );

      for (LineGroup group : groups)
      {
    	  buffer.append(group.toString());
      }
      
     if(buffer.length() > 0) 
        return buffer.toString();
      else 
        return null;
    }
    
    /**
     * Return a list of all line groups in the cache.
     * @return List
     */
    public synchronized List<LineGroup> getAllLineGroups()
    {
    	List<LineGroup> allGroups = new ArrayList<LineGroup>();
    	allGroups.addAll(misc.values());
    	allGroups.addAll(headers.values());
    	//allGroups.addAll(features.values());
    	
    	if (featureTable != null)
    		allGroups.add(featureTable);
    	
    	// Add sequence at the end
    	if (sequence != null)
    		allGroups.add(sequence);
    	
    	return allGroups;
    }
    
    /**
     * Clear all miscellaneous line groups.
     */
    public synchronized void clearMiscLineGroups()
    {
    	misc.clear();
    }
    
    /**
     * Clear the cache.
     */
    public synchronized void clear()
    {
    	headers.clear();
    	misc.clear();
    	sequence = null;
    	featureTable = null;
    }
    
    /**
     * Write all line groups out to a given stream in the correct order.
     * @param entry Entry
     * @param writer Writer
     * @throws IOException
     */
    public void writeToStream(final Entry entry, final Writer writer) throws IOException
    {
    	long num = 0;
    	boolean hasSequence = false;
    	boolean hasFeatureTable = false;
    	
    	for (java.util.Map.Entry<Long, LineGroup> miscGroup : misc.entrySet())
    	{
    		miscGroup.getValue().writeToStream(writer);
    		++num;
    	}
    	
    	for (java.util.Map.Entry<Long, LineGroup> headerGroup : headers.entrySet())
        {
        	headerGroup.getValue().writeToStream(writer);
        	++num;
        }
    	
    	if (featureTable != null)
    	{
    		featureTable.writeToStream(writer);
    		hasFeatureTable = true;
    	}
    	
    	if (sequence != null)
    	{
    		if(entry instanceof GFFDocumentEntry && sequence instanceof FastaStreamSequence)
    		{
    			LineGroup.writeStartOfGFFEntry(writer);
    		}
    		
    		sequence.writeToStream(writer);
    		hasSequence = true;
    	}
    	
    	if(num == 0L && (hasFeatureTable ^ hasSequence)) 
        {
    		// Don't write out the "//" end of entry marker if we have only one
    		// LineGroup (feature table or sequence)
    		// - this makes life easier for external programs that need
    		// to read the feature table or sequence
    		return;
        }
    	
    	if(num == 0L && (hasFeatureTable && hasSequence)) 
        {
    		// Don't write out the "//" end of entry marker if this is raw or
    		// FASTA sequence
    		if(sequence instanceof RawStreamSequence ||
    			sequence instanceof FastaStreamSequence) 
    			return;
        }

    	// Else write out the end delimiter
        if(entry instanceof PublicDBDocumentEntry) 
          LineGroup.writeEndOfEMBLEntry(writer);
    	
    }

}
