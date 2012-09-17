package uk.ac.sanger.artemis.components;

public interface IndexReferenceListener extends uk.ac.sanger.artemis.ChangeListener
{
  /**
   *  Invoked when a component scrolls or changes the scale.
   **/
  void indexReferenceChanged (IndexReferenceEvent event);
}
