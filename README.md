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


