#!/bin/bash
#Yang Jiang, 25.07.2024
#put this file aunder the same folder with R1.fq.gz R2.fq.gz
#Usage: ./MutSeq_analysis_YJ.sh ref.fa barcode regions_to_plot.bed num_threads


if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
    echo "Error: All arguments must be provided."
    echo "Usage: ./MutSeq_analysis_YJ.sh ref.fa barcode regions_to_plot.bed num_threads"
    exit 1
fi

for i in $(find ./ -maxdepth 1 -type f -name "*.fq.gz" | while read F; do basename -s .fq.gz $F | rev | cut -c 4- | rev; done | sort | uniq)

    do

        mkdir $i
        echo -e "sample_id\tbarcode\n"$i"_demux\t"$2"" > "$i"/"$i"_metadata.tsv

        #demultiplexing
        fqtk demux -i "$i"_R1.fq.gz "$i"_R2.fq.gz -r 150T 16M8B126T -o "$i"/ -s "$i"/"$i"_metadata.tsv
        
        #1st alingment
        bwa mem -v 2 -t "$4" $1 "$i"/"$i"_demux.R1.fq.gz "$i"/"$i"_demux.R2.fq.gz | samtools sort -@ "$4" -o "$i"/"$i".bam
        
        samtools index "$i"/"$i".bam
        
        #dedup using UMI
        umi_tools dedup -I "$i"/"$i".bam --paired -S "$i"/"$i".dedup.bam --umi-separator=":" --method=unique
        
        # name sort is needed for the next step
        samtools sort -@ "$4"  -n -o "$i"/"$i".nsorted.bam "$i"/"$i".dedup.bam
        samtools index "$i"/"$i".dedup.bam
        
        #remove the soft-clipped part of reads
        fgbio ClipBam --upgrade-clipping -i "$i"/"$i".nsorted.bam -o "$i"/"$i".clipped.bam -m "$i"/"$i"_clipping_mx.txt -r $1
        samtools sort -@ "$4" -o "$i"/"$i".clipped_sorted.bam "$i"/"$i".clipped.bam
        samtools index "$i"/"$i".clipped_sorted.bam
        
        #conver bam to fastq
        bedtools bamtofastq -i "$i"/"$i".clipped.bam -fq "$i"/"$i"_clean.R1.fq -fq2 "$i"/"$i"_clean.R2.fq
        
        #merge read pair
        bbmerge.sh in1="$i"/"$i"_clean.R1.fq in2="$i"/"$i"_clean.R2.fq out="$i"/"$i"_merged.fq outu="$i"/"$i"_unmerged.fq pfilter=1 ihist=ihist.txt
        
        #2nd alignment
        bwa mem -t "$4" $1 "$i"/"$i"_merged.fq | samtools sort -@ "$4" -o "$i"/"$i".re_aligned.bam
        samtools index "$i"/"$i".re_aligned.bam
        
        #get the per-base piled up information 
        samtools mpileup -d 0 -q 20 -Q 30 -f $1 "$i"/"$i".re_aligned.bam -l $3 --no-output-ends --no-output-ins --no-output-del > "$i"/"$i".mpileup.txt
        
        #optional lines for compress/ remove redundent files
        gzip "$i"/"$i"_unmerged.fq
        gzip "$i"/"$i"_merged.fq
        gzip "$i"/"$i".mpileup.txt
        
    


done
