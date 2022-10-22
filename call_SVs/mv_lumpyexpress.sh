FILE=$1
SbPATH='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/bam_files'
SPATH='/home/marvilla/Dropbox/Neurospora/mut_ac/SNV/lumpy'


#activate conda environmet
#conda init bash
#conda activate py2

#Stop at any error
set -ue

#Run tradictional lumpy
#On the sample
samtools view ${SbPATH}/${FILE}_al.bam \
   | tail -n+100000 \
   | ~/programs/lumpy-sv/scripts/pairend_distro.py \
    -r 101 \
    -X 4 \
    -N 10000 \
    -o ${SPATH}/${FILE}.lib1.histo > ${SPATH}/${FILE}_values

MEAN=$(grep -oP '(?<=mean:)[0-9\.]+' ${SPATH}/${FILE}_values)
STDEV=$(grep -oP '(?<=stdev:)[0-9\.]+' ${SPATH}/${FILE}_values)

#On the ancester
#samtools view ${SbPATH}/${ANC}_al.bam \
#    | tail -n+100000 \
#    | ~/programs/lumpy-sv/scripts/pairend_distro.py \
#    -r 101 \
#    -X 4 \
#    -N 10000 \
#    -o ${SPATH}/${ANC}.lib1.histo > ${SPATH}/${ANC}_values

#MEAN_anc=$(grep -oP '(?<=mean:)[0-9\.]+' ${SPATH}/${ANC}_values)
#STDEV_anc=$(grep -oP '(?<=stdev:)[0-9\.]+' ${SPATH}/${ANC}_values)

#Go to the binary
cd /home/marvilla/programs/lumpy-sv/bin/

./lumpy \
   -mw 4 \
   -tt 0 \
   -pe id:${FILE},bam_file:${SbPATH}/${FILE}.discordants.bam,histo_file:${SPATH}/${FILE}.lib1.histo,mean:$MEAN,stdev:$STDEV,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
   -sr id:${FILE},bam_file:${SbPATH}/${FILE}.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
   > ${SPATH}/${FILE}_lumpy.vcf

      
svtyper \
  -i ${SPATH}/${FILE}_lumpy.vcf \
  -B ${SbPATH}/${FILE}_al.sort.bam,${SbPATH}/${ANC}_al.sort.bam \
  -o ${SPATH}/${FILE}_lumpy_gt.vcf #, falta ancester file \
    

rm ${SPATH}/${FILE}.lib1.histo
#rm ${SPATH}/${ANC}.lib1.histo
rm ${SPATH}/${FILE}_values
#rm ${SPATH}/${ANC}_values

#cd /home/marvilla/programs/snpEff

#java -jar SnpSift.jar filter --file ${SPATH}/${FILE}_lumpytra_svtyper.vcf "( isVariant( GEN[${FILE}] ) && ( GEN[${ANC}].AO == 0 ) )" > ${SPATH}/${FILE}_somatic_tra.vcf

echo "Lumpy/svtyper done"

#-pe id:${ANC},bam_file:${SbPATH}/${ANC}.discordants.bam,histo_file:${SPATH}/${ANC}.lib1.histo,mean:$MEAN_anc,stdev:$STDEV_anc,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
#-sr id:${ANC},bam_file:${SbPATH}/${ANC}.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \

