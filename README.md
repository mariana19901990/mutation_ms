# mutation_ms
Data and scripts for manuscript: "Chromatin structure influences rate and spectrum of spontaneous mutations in Neurospora crassa" by Mariana Villalba de la Pena, Pauliina Summanen, Martta Liukkonen, and Ilkka Kronholm

This was a fairly large and complex project. I've separated the different pipelines, input files, and analyses into different files.

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
