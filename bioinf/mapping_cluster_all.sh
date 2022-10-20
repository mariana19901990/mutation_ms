#!/bin/bash

#Mapping all reads to the reference


#Path to the reference genome
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta
savedir=/scratch/project_2000350/genomics/mutacc/alignments

#List of sample names
#samples=( ANC1A ANC2a L1G40 L2G40 L3G40 L4G40 L5G40 L10G40 L13G40 L14G40 L15G40 L16G40 L17G40 L18G40 L19G40 L26G40 L27G40 L28G40 L29G40 L30G40 L31G40 L32G40 L33G40 L34G40 L35G40 L36G40 L37G40 L38G40 L39G40 L40G40 )
samples=( L6G40 L7G40 L8G40 L9G40 L11G40 L12G40 L20G40 L21G40 L22G40 L23G40 L24G40 L25G40 )

for s in ${samples[@]} 
do
fwd_fastq=/scratch/project_2000350/genomics/mutacc/raw/${s}_1.fq.gz
rev_fastq=/scratch/project_2000350/genomics/mutacc/raw/${s}_2.fq.gz
sample_name=${s}
#check that file names are corrent
#echo $fwd_fastq
#echo $rev_fastq
#echo $sample_name

### Begin mapping reads ###
echo Mapping reads of $sample_name againts reference genome
#Note that since this would generate a lot of intermediate files I'm overwriting files up until final step.

### Using BWA
bwa mem $reference $fwd_fastq $rev_fastq > $savedir/tempfile.PE.sam

#Convert to a .bam file
samtools view -S -b $savedir/tempfile.PE.sam > $savedir/tempfile.PE.bam
#Sorting reads
samtools sort $savedir/tempfile.PE.bam -o $savedir/tempfile.PE.sorted.bam
#index bam
samtools index $savedir/tempfile.PE.sorted.bam

### Fix mate information ###
picard FixMateInformation VALIDATION_STRINGENCY=LENIENT INPUT=$savedir/tempfile.PE.sorted.bam OUTPUT=$savedir/tempfile.mate_fixed.bam

#Sort the resulting bam
samtools sort $savedir/tempfile.mate_fixed.bam -o $savedir/tempfile.mate_fixed.sorted.bam

### Add or replace read-groups ###
#Note that sample names is used here in output and RGID!
picard AddOrReplaceReadGroups I=$savedir/tempfile.mate_fixed.sorted.bam O=$savedir/$sample_name.RG.bam SORT_ORDER=coordinate RGID=$sample_name RGLB=my_library RGPL=illumina RGSM=$sample_name RGPU=who_cares CREATE_INDEX=True VALIDATION_STRINGENCY=LENIENT

done

exit 0

