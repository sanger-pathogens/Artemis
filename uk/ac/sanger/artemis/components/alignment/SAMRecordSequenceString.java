package uk.ac.sanger.artemis.components.alignment;

/**
 * Utility class that handles operations on an alignment block
 * string taking into account the possibility that it could be
 * a secondary alignment.
 * 
 * @author kp11
 *
 */
public class SAMRecordSequenceString 
{

	public final static char SECONDARY_ALIGNMENT_BASE_CHAR = '=';
	public final static char SECONDARY_ALIGNMENT_MARKER = '*';
	
	private String sequence;
	private boolean isSecondaryAlignment;
	
	public SAMRecordSequenceString(String sequence) 
	{
		this.sequence = sequence;
		this.isSecondaryAlignment = (String.valueOf(SECONDARY_ALIGNMENT_MARKER).equals(sequence));
	}

	public String substring(int beginIdx, int endIdx) 
	{
		
		String result = null;
		
		if (!isSecondaryAlignment()) 
		{
			result = sequence.substring(beginIdx, endIdx);
		} 
		else
		{
			result =  new String(new char[endIdx-beginIdx]).replace('\0', SECONDARY_ALIGNMENT_BASE_CHAR);
		}
		
		return result;

	}
	
	public char charAt(int idx) 
	{
		
		char result;
		
		if (!isSecondaryAlignment()) 
		{
			result = sequence.charAt(idx);
		}
		else
		{
			result = SECONDARY_ALIGNMENT_BASE_CHAR;
		}
		
		return result;
	}

	public String getSequence() 
	{
		return sequence;
	}
	
	public int length() 
	{
		return sequence.length();
	}
	
	public boolean isSecondaryAlignment() {
		return isSecondaryAlignment;
	}
}
