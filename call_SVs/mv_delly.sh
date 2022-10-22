FILE=$1
#ANC=$2
#ANCESTER OPTIONS ANC1A OR ANC2a
S2PTH='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/delly'
S3PTH='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/bam_files'


#Stop at any error
set -ue

#Delly call use:
#delly call  -o t1.bcf -g hg19.fa tumor1.bam control1.bam

delly call -o ${S2PTH}/${FILE}.bcf -g ~/Dropbox/Neurospora/reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta ${S3PTH}/${FILE}_al.sort.bam 

#echo -e "${FILE}\ttumor \n${ANC}\tcontrol" > ${S2PTH}/${FILE}_samples.tsv

#delly filter -f somatic -p -o ${S2PTH}/${FILE}.somatic.bcf -s ${S2PTH}/${FILE}_samples.tsv ${S2PTH}/${FILE}.bcf

#Transform the bcf to vcf
bcftools view ${S2PTH}/${FILE}.bcf > ${S2PTH}/${FILE}.vcf

rm ${S2PTH}/${FILE}_samples.tsv

echo "delly All done"


