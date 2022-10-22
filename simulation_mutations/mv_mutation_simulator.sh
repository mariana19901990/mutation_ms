sample=$1

#Stop at any error
set -ue
#simulate, based in the RMT file, higher mutation rate in H3K9 domains
bash create_RMT2.sh ${sample} > rmtfile
mutation-simulator neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta rmt rmtfile
mv neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_ms.fasta neurospora_crassa_ms_${sample}.fasta
mv neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_ms.vcf neurospora_crassa_ms_${sample}.vcf

echo "RMT simulation done"

# Simulate SNPs uniformly across all sequences
mutation-simulator neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta args -sn 0.000002 -titv 1.08 -s neurospora_crassa -n ${sample} 
mv neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_ms.fasta neurospora_crassa_ms_${sample}.fasta
mv neurospora_crassa_or74a_12_supercontigs_mtDNA_mata_ms.vcf neurospora_crassa_ms_${sample}.vcf

echo "uniformly simulation done"

#Simulate reads
dwgsim -N 10000 -C 30 -c 0 -S 2 -r 0 -F 0 -R 0 -1 150 -2 150 -e 0.003 -E 0.003 -Q 2 neurospora_crassa_ms_${sample}.fasta ${sample}

echo "read simulation done"
