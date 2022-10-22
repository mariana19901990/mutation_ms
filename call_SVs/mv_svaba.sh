SAMPLE=$1
ANC=$2
S_PATH='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/bam_files'
REF=/home/marvilla/Dropbox/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta

#Stop at any error
set -ue

#Go the Svaba directory
cd /home/marvilla/Dropbox/Neurospora/mut_ac/SNV/SvABA

#Run SvABA
#Usage: svaba run -t <BAM/SAM/CRAM> -G <reference> -a myid [OPTIONS]
#svaba run -t ${S_PATH}/${SAMPLE}_al.sort.bam -n ${S_PATH}/${ANC}_al.sort.bam -a ${SAMPLE} -G ${REF}

#run without the ancester
svaba run -t ${S_PATH}/${SAMPLE}_al.sort.bam -a ${SAMPLE} -G ${REF}
echo "SvABA done"
