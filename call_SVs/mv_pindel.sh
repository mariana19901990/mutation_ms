SAMPLE=$1
ANC=$2
SPATH='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/bam_files'
PATH_p='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/pindel'
#stop at any error
set -ue

cd ${PATH_p}
#create a config file
INS_SAMPLE=$(samtools stats ${SPATH}/${SAMPLE}_al.bam | grep 'insert size average' |grep -oP '(?<=insert size average:\t)[0-9]+')

INS_ANC=346

#echo -e "/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/bam_files/${SAMPLE}_al.sort.bam \t${INS_SAMPLE} \t${SAMPLE}\n/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/bam_files/${ANC}_al.sort.bam \t${INS_ANC} \t${ANC}" > ${SAMPLE}.config.txt
#Just the sample
echo -e "/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/bam_files/${SAMPLE}_al.sort.bam \t${INS_SAMPLE} \t${SAMPLE}" > ${SAMPLE}.config.txt

pindel -f ~/Dropbox/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta -i ${SAMPLE}.config.txt -c ALL -o ${SAMPLE} -g 

pindel2vcf -P ${SAMPLE} -r ~/Dropbox/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta -R neurospora_crassa_or74a -d 20101123 -v ${SAMPLE}.vcf

#346
#319
