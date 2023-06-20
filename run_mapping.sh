##@param1: reference fasta sequence filename
##@param2: fastq1 file
##@param3: fastq2 file
##@param4: name of aligners: biscuit/bismark/bwa-meth/bsmap
##@param5: output folder name


##
ref=$1
fastq1=$2
fastq2=$3
aligner=$4
output_folder=$5
if [ "$aligner" == "bismark" ]
then
    echo "Running bismark ..."
        bismark_folder="/common/bermanblab/bin/wgbs/bismark_v0.14.5"
        # mapping bismark 10 mismatches per 100bp --score_min L,0,-0.8 : bowtie2
        $bismark_folder/bismark $ref --score_min L,0,-0.8 -q -p 12 -1 $fastq1  -2 $fastq2 -X 1000 --multicore 8 --unmapped --ambiguous -o $sample"_BISMARK"
    /common/bermanblab/bin/samtools flagstat $sample"_BISMARK"/$sample"_R1_001_val_1.fq.gz_bismark_bt2_pe.bam" > $sample"_R1_001_val_1.fq.gz_bismark_bt2_pe.bam.flagstat.metric.txt"
        $bismark_folder/deduplicate_bismark --bam $fastq1"_bismark_bt2_pe.bam"
        $bismark_folder/bismark_methylation_extractor --cytosine_report --paired-end --comprehensive  $fastq1"_bismark_bt2_pe.deduplicated.bam" --genome_folder $ref
elif [ "$aligner" == "bwa-meth" ]
then
    bwa_meth_folder="/common/bermanblab/bin/wgbs"
    #python $bwa_meth_folder/bwameth.py --threads 16 --reference $ref --prefix $output_folder/$output_folder  $fastq1 $fastq2
    /common/bermanblab/bin/samtools flagstat $output_folder/$output_folder".bam" > $output_folder/$output_folder".bam.flagstat.metric.txt"
elif [ "$aligner" == "biscuit" ]
then
    biscuit_path=/common/bermanblab/bin/wgbs/biscuit-master/bin/biscuit
    $biscuit_path align -t 8 $ref $fastq1 $fastq2 > $output_folder/$output_folder".bam"
    /common/bermanblab/bin/samtools flagstat $output_folder/$output_folder".bam" > $output_folder/$output_folder".bam.flagstat.metric.txt"
else
    echo "Running bsmap"
    bsmap_folder="/common/bermanblab/bin/wgbs/bsmap-2.90/"
    $bsmap_folder/bsmap -a $fastq1 -b $fastq2 -o $fastq1"_bsmap.bam" -d $ref -s 16 -V 2 -v 10
fi
