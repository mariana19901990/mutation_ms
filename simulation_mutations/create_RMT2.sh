sample=$1

#Stop at any error
set -ue

#create the RMT file
#get all H3K9 sites
chr1_K9=$(cat H3K9_chr1)
chr2_K9=$(cat H3K9_chr2)
chr3_K9=$(cat H3K9_chr3)
chr4_K9=$(cat H3K9_chr4)
chr5_K9=$(cat H3K9_chr5)
chr6_K9=$(cat H3K9_chr6)
chr7_K9=$(cat H3K9_chr7)

#get all H3K36 sites
#chr1_euchromatin=$(cat euchromatin_chr1)
#chr2_euchromatin=$(cat euchromatin_chr2)
#chr3_euchromatin=$(cat euchromatin_chr3)
#chr4_euchromatin=$(cat euchromatin_chr4)
#chr5_euchromatin=$(cat euchromatin_chr5)
#chr6_euchromatin=$(cat euchromatin_chr6)
#chr7_euchromatin=$(cat euchromatin_chr7)
#meta-information
echo "fasta=neurospora_crassa_or74a_12_supercontigs_mtDNA_mata.fasta"
echo "md5=c1c5127702f7c5de0a6ec16684498c29"
echo "species_name=Neurospora crassa"
echo "sample_name=${sample}"
echo "titv=1.08"
#standard
echo "std"
echo "it None"
echo "sn 0.000003"
#range definitions
echo "chr 1"
echo "it None"
echo "$chr1_K9"
#echo "$chr1_euchromatin"
echo "chr 2"
echo "it None"
echo "$chr2_K9"
#echo "$chr2_euchromatin"
echo "chr 3"
echo "it None"
echo "$chr3_K9"
#echo "$chr3_euchromatin"
echo "chr 4"
echo "it None"
echo "$chr4_K9"
#echo "$chr4_euchromatin"
echo "chr 5"
echo "it None"
echo "$chr5_K9"
#echo "$chr5_euchromatin"
echo "chr 6"
echo "it None"
echo "$chr6_K9"
#echo "$chr6_euchromatin"
echo "chr 7"
echo "it None"
echo "$chr7_K9"
#echo "$chr7_euchromatin"