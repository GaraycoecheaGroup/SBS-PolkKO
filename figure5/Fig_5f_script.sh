#!/bin/bash
#Yang Jiang, 2024

#Usage: ./MutSeq_analysis_YJ.sh input_R1.fq input_R2.fq ref.fa barcode.txt num_threads
    #order of input can not be changed
    #input fq: fastq files before multiplexing
    #barcode.txt: tsv with sample_id (column 1) and barcode (column 2), the header is needed, the file need to end with an empty line
    #parameter for fqtk demux may need to change depende on experimental design


#loop for processing each sample
tail -n +2 $4 | while IFS=$'\t' read -r row _; do
  
  mkdir "$row"
  mkdir "$row"/demux_output

  fqtk demux -i $1 $2 -r 150T 16M8B126T -o "$row"/demux_output/ -s $4 

  #1st alingment
  bwa mem -t "$5" $3 "$row"/demux_output/"$row".R1.fq.gz "$row"/demux_output/"$row".R2.fq.gz | samtools sort -@ "$5" -o "$row"/"$row".bam
  
  samtools index "$row"/"$row".bam

  #dedup using UMI
  umi_tools dedup -I "$row"/"$row".bam --paired -S "$row"/"$row".dedup.bam --umi-separator=":" --method=unique
  
  # name sort is needed for the next step
  samtools sort -@ "$5"  -n -o "$row"/"$row".nsorted.bam "$row"/"$row".dedup.bam
  samtools index "$row"/"$row".dedup.bam

  #remove the soft-clipped part of reads
  fgbio ClipBam --upgrade-clipping -i "$row"/"$row".nsorted.bam -o "$row"/"$row".clipped.bam -m "$row"/"$row"_clipping_mx.txt -r $3
  samtools sort -@ "$5" -o "$row"/"$row".clipped_sorted.bam "$row"/"$row".clipped.bam
  samtools index "$row"/"$row".clipped_sorted.bam

  #conver bam to fastq
  bedtools bamtofastq -i "$row"/"$row".clipped.bam -fq "$row"/"$row"_clean.R1.fq -fq2 "$row"/"$row"_clean.R2.fq

  #merge read pair
  bbmerge.sh in1="$row"/"$row"_clean.R1.fq in2="$row"/"$row"_clean.R2.fq out="$row"/"$row"_merged.fq outu="$row"/"$row"_unmerged.fq pfilter=1 ihist=ihist.txt
  gzip "$row"/"$row"_merged.fq

  rm "$row"/"$row"_unmerged.fq "$row"/"$row".bam "$row"/"$row".bam.bai "$row"/"$row".dedup.bam "$row"/"$row".dedup.bam.bai "$row"/"$row".clipped.bam "$row"/"$row"_clean.R2.fq "$row"/"$row"_clean.R1.fq "$row"/"$row"_clipping_mx.txt ihist.txt "$row"/"$row".clipped_sorted.bam "$row"/"$row".clipped_sorted.bam.bai "$row"/"$row".nsorted.bam

done