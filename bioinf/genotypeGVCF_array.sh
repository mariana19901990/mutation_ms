#!/bin/bash

#Using haplotype caller to call genotypes one by one

#Path to the reference genome
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta
outputdir=/scratch/project_2000350/genomics/mutacc/genotyping
variantdb=/scratch/project_2000350/genomics/mutacc/genotyping/genoDB
intervals=/scratch/project_2000350/genomics/mutacc/genotyping/interval_list #List of genomic intervals, used in parallelization

#Paths to programs
path_GATK=/projappl/project_2000350/Genomics/gatk-4.2.0.0


#Note annoying syntax '-V gendb:///home/path/myDB' is correct three (3) forward slashes!
$path_GATK/gatk --java-options "-Xmx4g" GenotypeGVCFs  \
   -R $reference -V gendb://$variantdb --sample-ploidy 2 -all-sites -L $intervals/${1}.interval_list -O $outputdir/${1}.vcf.gz

#Note that resulting vcf.files need to be combined and indexed

