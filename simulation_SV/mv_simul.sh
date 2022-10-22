SPTH='/home/marvilla/programs/SURVIVOR/Debug'
SPTH_1='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/survivor'
SPTH_r='/home/marvilla/Dropbox/Neurospora/mut_ac/raw_data'
SPTH_b='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/bam_files'

#Stop at any error
set -ue
 
cd ${SPTH}

#Create simulated data (germline)
#create a parameter file
#./SURVIVOR simSV parameter_file

#18SV
./SURVIVOR simSV ~/Dropbox/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta parameter_file 0.001 0 ${SPTH_1}/simulated_neurospora

#read simulator Using dwgsim
dwgsim -N 10000 -C 30 -c 0 -S 2 -r 0 -F 0 -R 0 -1 150 -2 150 -e 0.003 -E 0.003 -Q 2 ${SPTH_1}/simulated_neurospora.fasta ${SPTH_r}/dwg

##########################################
#30SV
./SURVIVOR simSV ~/Dropbox/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta parameter_file_1 0.001 0 ${SPTH_1}/simulated_neurospora_1

#read simulator Using dwgsim
dwgsim -N 10000 -C 30 -c 0 -S 2 -r 0 -F 0 -R 0 -1 150 -2 150 -e 0.003 -E 0.003 -Q 2 ${SPTH_1}/simulated_neurospora_1.fasta ${SPTH_r}/dwg_1

############################################
#40SV
./SURVIVOR simSV ~/Dropbox/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta parameter_file_2 0.001 0 ${SPTH_1}/simulated_neurospora_2

#read simulator Using dwgsim
dwgsim -N 10000 -C 30 -c 0 -S 2 -r 0 -F 0 -R 0 -1 150 -2 150 -e 0.003 -E 0.003 -Q 2 ${SPTH_1}/simulated_neurospora_2.fasta ${SPTH_r}/dwg_2

##########################################
#50SV
./SURVIVOR simSV ~/Dropbox/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta parameter_file_3 0.001 0 ${SPTH_1}/simulated_neurospora_3

#read simulator Using dwgsim
dwgsim -N 10000 -C 30 -c 0 -S 2 -r 0 -F 0 -R 0 -1 150 -2 150 -e 0.003 -E 0.003 -Q 2 ${SPTH_1}/simulated_neurospora_3.fasta ${SPTH_r}/dwg_3

##########################################
#100SV
./SURVIVOR simSV ~/Dropbox/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta parameter_file_4 0.001 0 ${SPTH_1}/simulated_neurospora_4

#read simulator Using dwgsim
dwgsim -N 10000 -C 30 -c 0 -S 2 -r 0 -F 0 -R 0 -1 150 -2 150 -e 0.003 -E 0.003 -Q 2 ${SPTH_1}/simulated_neurospora_4.fasta ${SPTH_r}/dwg_4
