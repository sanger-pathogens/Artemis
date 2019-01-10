#! /bin/bash

rm -f RPKM_sorted_test_1.bam RPKM_sorted_test_1.bam.bai
rm -f RPKM_sorted_test_2.bam RPKM_sorted_test_2.bam.bai

echo "Creating test bam file 1"
samtools view -bS RPKM_test_1.sam > RPKM_test_1.bam
samtools sort -T /tmp/tmp1.sorted -o RPKM_sorted_test_1.bam RPKM_test_1.bam
samtools index RPKM_sorted_test_1.bam RPKM_sorted_test_1.bam.bai

echo "Creating test bam file 2"
samtools view -bS RPKM_test_2.sam > RPKM_test_2.bam
samtools sort -T /tmp/tmp2.sorted -o RPKM_sorted_test_2.bam RPKM_test_2.bam
samtools index RPKM_sorted_test_2.bam RPKM_sorted_test_2.bam.bai

echo "Finished"

