package uk.ac.sanger.artemis.components.alignment;

import net.sf.samtools.SAMRecord;


class BamViewRecord
{
  protected SAMRecord sam;
  protected short bamIndex = -1;
  BamViewRecord(final SAMRecord sam, final short bamIndex)
  {
    this.sam = sam;
    this.bamIndex = bamIndex;
  }
}