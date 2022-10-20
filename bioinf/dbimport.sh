#!/bin/bash

outputdir=/scratch/project_2000350/genomics/mutacc/genotyping/genoDB
inputdir=/scratch/project_2000350/genomics/mutacc/genotyping
path_intervals=/scratch/project_2000350/genomics/Neurospora_reference/all.list

#Paths to programs
path_GATK=/projappl/project_2000350/Genomics/gatk-4.2.0.0


$path_GATK/gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
      -V $inputdir/ANC1A.g.vcf.gz \
      -V $inputdir/ANC2a.g.vcf.gz \
      -V $inputdir/L1G40.g.vcf.gz \
      -V $inputdir/L2G40.g.vcf.gz \
      -V $inputdir/L3G40.g.vcf.gz \
      -V $inputdir/L4G40.g.vcf.gz \
      -V $inputdir/L5G40.g.vcf.gz \
      -V $inputdir/L10G40.g.vcf.gz \
      -V $inputdir/L13G40.g.vcf.gz \
      -V $inputdir/L14G40.g.vcf.gz \
      -V $inputdir/L15G40.g.vcf.gz \
      -V $inputdir/L16G40.g.vcf.gz \
      -V $inputdir/L17G40.g.vcf.gz \
      -V $inputdir/L18G40.g.vcf.gz \
      -V $inputdir/L19G40.g.vcf.gz \
      -V $inputdir/L26G40.g.vcf.gz \
      -V $inputdir/L27G40.g.vcf.gz \
      -V $inputdir/L28G40.g.vcf.gz \
      -V $inputdir/L29G40.g.vcf.gz \
      -V $inputdir/L30G40.g.vcf.gz \
      -V $inputdir/L31G40.g.vcf.gz \
      -V $inputdir/L32G40.g.vcf.gz \
      -V $inputdir/L33G40.g.vcf.gz \
      -V $inputdir/L34G40.g.vcf.gz \
      -V $inputdir/L35G40.g.vcf.gz \
      -V $inputdir/L36G40.g.vcf.gz \
      -V $inputdir/L37G40.g.vcf.gz \
      -V $inputdir/L38G40.g.vcf.gz \
      -V $inputdir/L39G40.g.vcf.gz \
      -V $inputdir/L40G40.g.vcf.gz \
      -L $path_intervals --genomicsdb-workspace-path $outputdir
