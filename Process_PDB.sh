#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100g
#SBATCH --time=5-20:00:00
#SBATCH --job-name=PDB_Process
#SBATCH --output=jobname.out.%j
#SBATCH --partition=bmh
eval "$(conda shell.bash hook)"
conda activate PDB_Tools
pathtorepo="./"
gff="./TestData/Testset.gff"
CDS="./TestData/Testset_cds.fasta"
transcipt_list="./TestData/Testset_transcripts.txt"
PDB_dir="./TestData/PDB"

#### If you have a snpEFF annotated population VCF set to TRUE and provide path to VCF, if not put False
PopData="TRUE"
snpEffVCF="./TestData/Testset.vcf.gz"
threads=10


###Now that arguments are set we can run through the pipeline



##Combine freq files with snpEFF string
vcftools --gzvcf $snpEffVCF --freq 
zcat $snpEffVCF | cut -f 8 |  grep -v "#" > snpstring.txt
grep -v "#" out.frq | cut -f 1,2 > freq_subcol_left.txt
grep -v "#" out.frq | cut -f 6-  > freq_subcol_right.txt
paste freq_subcol_left.txt snpstring.txt freq_subcol_right.txt > Freq_snpEFF.txt
rm freq_subcol_left.txt freq_subcol_right.txt snpstring.txt

###Initialize directories
mkdir -p prot_files
mkdir -p Freq
mkdir -p dssp
mkdir -p CDS
mkdir -p prot_files_mutation
mkdir -p gff

#subset freq table, cds, gff for each transcripta
if [ "$PopData" = "TRUE" ]; then
	cat $transcipt_list | xargs -P $threads -I % bash $pathtorepo/subsetFreq.sh % Freq_snpEFF.txt
fi
cat $transcipt_list | xargs -P $threads -I % bash $pathtorepo/subsetCDS.sh % $CDS
cat $transcipt_list | xargs -P $threads -I % bash $pathtorepo/subsetgff.sh % $gff
#Calculate surface area of each Amino ACis
conda activate dssp
cat $transcipt_list | xargs -P $threads -I % bash $pathtorepo/mkdssp.sh % $PDB_dir
conda deactivate
#Run the linchpin R script that Matches all input statistics into a nice table
cat $transcipt_list | xargs -P $threads -I % Rscript $pathtorepo/PDB_stats.R % $PDB_dir $PopData


if [ "$PopData" = "TRUE" ]; then
#Get summary table of Population data & Protein Structure
	cat prot_files_mutation/* | head -n 1   > Protein_popstats.txt
	cat prot_files_mutation/* | grep -v RefAllele >> Protein_popstats.txt
	bgzip -f Protein_popstats.txt
fi

#Get summary table of Protein Structure Info
cat prot_files/* |head -n 1  > Protein_prot_files_mutation.txt
cat prot_files/* | grep -v RefAllele >> Protein_prot_files_mutation.txt
bgzip -f Protein_prot_files_mutation.txt

