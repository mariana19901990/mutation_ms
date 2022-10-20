#!/bin/bash -l
#SBATCH --job-name=haplotypecaller
#SBATCH --account=project_2000350
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH --array=0-3
#SBATCH --mail-type=END
#SBATCH --mail-user=ilkka.kronholm@jyu.fi

module load biokit

#Using haplotype caller to call genotypes one by one

#Path to the reference genome
reference=/scratch/project_2000350/genomics/Neurospora_reference/neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta
inputdir=/scratch/project_2000350/genomics/natpop/alignments
outputdir=/scratch/project_2000350/genomics/natpop/genotyping

#Paths to programs
path_GATK=/projappl/project_2000350/Genomics/gatk-4.2.0.0

#samples=(ANC1A ANC2a L1G40 L2G40 L3G40 L4G40 L5G40 L10G40 L13G40 L14G40 L15G40 L16G40 L17G40 L18G40 L19G40 L26G40 L27G40 L28G40 L29G40 L30G40 L31G40 L32G40 L33G40 L34G40 L35G40 L36G40 L37G40 L38G40 L39G40 L40G40)
#samples=(L6G40 L7G40 L8G40 L9G40 L11G40 L12G40 L20G40 L21G40 L22G40 L23G40 L24G40 L25G40)
#samples=(10948 10886 10932 1165 4498 8816 3223 8845 10908 847 10904 851 1131 8850 8819 4708 4712 6203 4824 8783 8790 3975 10928 10912 4494 3210 10923 10950 10951 10946 3211 10906 P4452 P4463 P4468 P4471 P4476 P4479 10882 10883 10884 10892 10907 10914 10915 10918 10925 10926 10927 10935 10937 10943 10983 3943 P4489)
#samples=(10948 10932 8816 8845 10904 4708 3975 3210 10946 P4452 P4476 10882 10935 3943)
#samples=(P4452)
samples=(5910 4730 1133 4716) #Extra samples

sample=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based

$path_GATK/gatk --java-options "-Xmx8g" HaplotypeCaller  \
   -R $reference -I $inputdir/$sample.RG.bam --sample-ploidy 2 -mbq 10 --output-mode EMIT_ALL_CONFIDENT_SITES -O $outputdir/$sample.g.vcf.gz -ERC GVCF 
