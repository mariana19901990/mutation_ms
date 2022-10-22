FILE=$1
#SPATH='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/pindel'
SPATH='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/lumpy'
cd /home/marvilla/programs/snpEff

#This is to filter the lumpy results
#java -jar SnpSift.jar filter --file ${SPATH}/${FILE}_lumpy_gt.vcf "(GEN[0].GT != GEN[1].GT) && (GEN[0].DP >= 10 && GEN[*].DP >= 10 && GEN[1].AO == 0 && GEN[*].GQ >= 25)"  > ${SPATH}/${FILE}_somatic.vcf 

#This is to filter without the ancester
java -jar SnpSift.jar filter --file ${SPATH}/${FILE}_lumpy_gt.vcf "(GEN[*].DP >= 10 && GEN[*].DP >= 10 && GEN[*].AO == 0 && GEN[*].GQ >= 25)"  > ${SPATH}/${FILE}_somatic.vcf 

#- Only consider variants that have sufficient depth (DP) in tumor and
#normal. This number will depend on your sequencing depth.
#- Remove all variants that have any alt (AO) evidence in the normal.
#- Keeps variants with enough depth in the tumor. This number will depend on
#your sequencing depth and sample purity.

#This is to filter the pindel results
#java -jar SnpSift.jar filter --file ${SPATH}/${FILE}.vcf "(SVLEN >= 30 && SVLEN < 1000000) && (GEN[0].GT != GEN[1].GT) && (GEN[0].AD[*] >= 10) && (GEN[1].AD[*] >= 10) && (GEN[1].GT != '0/0')" > ${SPATH}/${FILE}_filter.vcf 


