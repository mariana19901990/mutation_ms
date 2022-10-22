SAMPLE=$1
SPTH='/home/marvilla/programs/SURVIVOR/Debug'
SPTH_1='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/survivor'
SPTH_d='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/delly'
SPTH_l='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/lumpy'
SPTH_p='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/pindel' 
SPTH_s='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/SvABA'

#Stop at any error
set -ue
 
cd ${SPTH}
#USE ./SURVIVOR eval caller.vcf simulated.bed 10 eval_res

./SURVIVOR eval ${SPTH_d}/${SAMPLE}.vcf ${SPTH_1}/simulated_neurospora.bed 50 ${SPTH_1}/eval_delly_${SAMPLE}
./SURVIVOR eval ${SPTH_l}/${SAMPLE}_lumpy_gt.vcf ${SPTH_1}/simulated_neurospora.bed 50 ${SPTH_1}/eval_lumpy_${SAMPLE}
./SURVIVOR eval ${SPTH_s}/${SAMPLE}.svaba.sv.vcf ${SPTH_1}/simulated_neurospora.bed 50 ${SPTH_1}/eval_svaba_${SAMPLE}
./SURVIVOR eval ${SPTH_p}/${SAMPLE}.vcf ${SPTH_1}/simulated_neurospora.bed 50 ${SPTH_1}/eval_pindel_${SAMPLE}
#./SURVIVOR eval ${SPTH_1}/${SAMPLE}_merged.vcf ${SPTH_1}/simulated_neurospora_${NUM}.bed 50 ${SPTH_1}/eval_merged_${SAMPLE}

echo 'evaluation done'
#somatic_vcf
