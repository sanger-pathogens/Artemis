
package uk.ac.sanger.artemis.io;

import java.lang.IndexOutOfBoundsException;
import java.lang.String;
import java.util.Iterator;

import org.biojava.bio.symbol.Edit;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.bio.seq.io.SymbolListCharSequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.seq.DNATools;
import org.biojava.utils.ChangeVetoException;

import uk.ac.sanger.artemis.util.ReadOnlyException;

public class BioJavaSequence implements Sequence
{
    private org.biojava.bio.symbol.SymbolList symbols;

    private int aCount;
    private int cCount;
    private int gCount;
    private int tCount;

    private boolean sequenceIsStale = true;

    public BioJavaSequence(final org.biojava.bio.symbol.SymbolList symbols)
    {
	this.symbols = symbols;
    }
    
    public int length()
    {
	return symbols.length();
    }

    public int getACount()
    {
	if (sequenceIsStale)
	    countSymbols();

	return aCount;
    }

    public int getCCount()
    {
	if (sequenceIsStale)
	    countSymbols();

	return cCount;
    }

    public int getGCount()
    {
	if (sequenceIsStale)
	    countSymbols();

	return gCount;
    }

    public int getTCount()
    {
	if (sequenceIsStale)
	    countSymbols();

	return tCount;
    }

    public int getOtherCount()
    {
	if (sequenceIsStale)
	    countSymbols();

	return (symbols.length() - aCount - cCount - gCount -tCount);
    }

  org.biojava.bio.symbol.SymbolList getSymbolList () {
    return symbols;
  }

    public String getSubSequence(int index1, int index2)
    {
	String subSeq = "";

	try
	{
	    subSeq = symbols.subStr(index1, index2);
	}
	catch (IndexOutOfBoundsException ioe)
	{
	    System.err.println("An error occurred while extracting subsequence "
			       + ioe.getMessage());
	    ioe.printStackTrace();
	}

	return subSeq;
    }

    public void setFromString(String seqString)
        throws ReadOnlyException, IllegalSymbolException
    {
      Edit ed = new Edit(1, length (), DNATools.createDNA(seqString));
      try {
        symbols.edit(ed);
      } catch (ChangeVetoException e) {
        throw new ReadOnlyException ("cannot set sequence - readonly");
      } catch (IllegalAlphabetException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      }
    }

    private void countSymbols()
    {
      final SymbolListCharSequence slcs = new SymbolListCharSequence (symbols);
      
      int a = 0;
      int c = 0;
      int g = 0;
      int t = 0;

      for (int i = 0 ; i < slcs.length () ; ++i) {
        char token = slcs.charAt(i);

	    switch (token)
	    {
		case 'a':
		    a++;
		    break;

		case 'c':
		    c++;
		    break;

		case 'g':
		    g++;
		    break;

		case 't':
		    t++;
		    break;

		default:
		    break;
	    }
	}

	aCount = a;
	cCount = c;
	gCount = g;
	tCount = t;

        sequenceIsStale = false;
    }
}
