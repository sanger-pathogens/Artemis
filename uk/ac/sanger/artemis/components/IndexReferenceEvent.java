package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.ChangeEvent;

public class IndexReferenceEvent extends ChangeEvent 
{
  private static final long serialVersionUID = 1L;

  public IndexReferenceEvent(Object source)
  {
    super(source);
  }

}
