# mutation_ms
Data and scripts for manuscript: "Chromatin structure influences rate and spectrum of spontaneous mutations in Neurospora crassa" by Mariana Villalba de la Pena, Pauliina Summanen, Martta Liukkonen, and Ilkka Kronholm

Preprint is available in [bioRxiv](https://doi.org/10.1101/2022.03.13.484164)

This was a fairly large and complex project. I've separated the different pipelines, input files, and analyses into different files.

## Bioinformatics pipeline to map reads and call mutations

We used a national supercluster (CSC) to process the short read sequencing data. The script files have been mainly written so that they work on the cluster. If you want to use them, you have to modify them so that they work in your environment. Nevertheless you can extract the GATK, BWA etc. commands and run them on your system.

The script files are in folder /bioinf
1. Map short reads to the reference genome: file /bioinf/mapping_cluster_all.sh contains the commands used for BWA and post processing of the files
2. Use haplotypecaller in the GATK pipeline to call genotypes for each sample (GVCF file): file /bioinf/call_genotype.sh contains the haplotypecaller commands
3. Consolidate the samples into a databse: file /bioinf/dbimport.sh contains the commands to make the database
4. Jointly call genotypes from the sample database: file /bioinf/genotypeGVCF_array.sh contains commands to call genotypes jointly, chromsome by chromosome for parallelization purposes
5. Postprocess the resulting VCF files: file /bioinf/postprocess_vcf.sh

Further VCF post processing
 - Extract .gz file

 - Remove INFO.LEAC, INFO:MLEAF columns from vcf file
 
 	~/bcftools-1.11/bcftools annotate -x INFO/MLEAC,INFO/MLEAF all.samples.vcf -o all.samples.rm.vcf

 - For nat pop also run (PID column too long, and Neurospora is haploid anyway)
 
	~/bcftools-1.11/bcftools annotate -x FORMAT/PID all.samples.rm.vcf -o all.samples.rm2.vcf
	
	~/bcftools-1.11/bcftools annotate -x INFO/MLEAC,INFO/MLEAF,FORMAT/PID all.samples.vcf -o all.samples.natpop_comp.vcf

6. Then make an indexed database from the resulting VCF file using wormtable
 - Convert to wormtable database

	vcf2wt --truncate all.samples.rm.vcf allsamples.wt
	
	wtadmin add allsamples.wt CHROM+POS

  - Can look what are the contents of wormtable columns with
        wtadmin show allsamples.wt

7. Filter mutations called by the GATK pipeline: file /bioinf/mutations_neuro_MA.py contains a python scripts that uses wormtable to filter for candidate mutations
8. Genotyping of natural population samples uses the same pipeline but for step 7. another script is used: /bioinf/natpop_allpos.py to call all sites in the sample of natural populations

File /bioinf/retrieve.py contains a script that is called from the R-script MAanalysis.R that retrieves the trinucleotide context the mutation. Used the generated wormtable database 

## Estimating the number of mitoses

The analysis scripts and data files for estimating the number of mitoses that the MA lines went through are in folder /mitoses
- File /mitoses/mitoses.R contains R scripts used to estimate the number of mitoses
- File /mitoses/colonies_on_plates2.csv is a data file for estimating the number of conidia produced by cultures
- File /mitoses/sorbose_colonies.csv is a data file for sizes of colonies on sorbose plates
- File /mitoses/nuclei_counts.csv is a data file for numbers of nuclei within the mycelium
- File /mitoses/tt_diameters.csv is a data file for test tube diameters

## Statistical analysis of mutation data

- File curated_mutations_final.csv contains the mutations that occurred in the MA lines. This is the final set of mutations after manual inspection in IGV and Sanger confirmation.
- File mutationscripts.R contains R scripts that are needed to process some of the mutation data

Folder /data contains .RData files that contain intermediate steps and results of statistical analyses.

## Chromatin modification data
Data files for chromatin modifications are in folder /chromatin
- File /chromatin/duplicated.csv contains data about the duplicated regions defined by Wang et al. 2020
- File /chromatin/centromeres.csv contains data about the centromere locations
- File /chromatin/2489.H3K27.domains.duprm.bed is a bed file that contains regions where H3K27me3 occurs
- File /chromatin/2489.H3K27_exK9.bed is a bed file that contains regions of H3K27me3, where regions overlapping H3K9me3 have been excluded
- File /chromatin/2489.H3K9.domains.duprm.bed is a bed file that contains regions where H3K9me3 occurs
- File /chromatin/2489.H3K36.domains.duprm.bed is a bed file that contains regions where H3K36me2 occurs
- File /chromatin/2489.euchromatin.bed is a bed file that contains regions of the genome where H3K9me3, H3K27me3 have been excluded

## Trinucleotide count data
Data files for trinucleotide counts for different regions of the genome
- File Nc_trinuc_freq_whole_genome.txt contains trinucleotide counts over the whole genome
- File Nc_trinuc_freq_euchromatin.txt contains trinucleotide counts for euchromatic regions
- File Nc_trinuc_freq_centromere.txt contains trinucleotide counts for centromeric regions
- File Nc_trinuc_freq_H3K27_exK9.txt contains trinucleotide counts for regions where only H3K27me3 occurs
- File Nc_trinuc_freq_H3K9.txt contains trinucleotide counts for regions where H3K9me3 occurs

## Data for called sites
The original data files for the numbers of called sites are very large. The file /data/calledsites.RData contains the number of called sites that has been processed and can be directly loaded into R
