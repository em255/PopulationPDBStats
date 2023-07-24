conda activate PDB_Tools
pathtorepo="path/to/repo"
gff="path/to/transcripts.gff"
snpEffVCF="path/to/snpEFF.vcf.gz"
transcipt_list="path/to/transcipt_list.txt"
freqsnpEff="path/to/snpEFF_Freqfile"
threads=10

##Make freq files with snpEFF string
vcftools --gzvcf $snpEffVCF --freq 
cut -f 8 $snpEffVCF |  grep -v "#" > snpstring.txt
grep -v "#" out.frq | cut -f 1,2 > freq_subcol_left.txt
grep -v "#" out.frq | cut -f 6-  > freq_subcol_right.txt
paste freq_subcol_left.txt snpstring.txt freq_subcol_right.txt > Freq_snpEFF.txt
rm freq_subcol_left.txt freq_subcol_right.txt snpstring.txt

mkdir -p prot_files
mkdir -p Freq
mkdir -p dssp
mkdir -p CDS
mkdir -p pLDDT
mkdir -p gff
cat $transcipt_list | xargs -P $threads -I % bash $pathtorepo/subsetFreq.sh % Freq_snpEFF.txt
cat $transcipt_list | xargs -P $threads -I % bash $pathtorepo/subsetCDS.sh %
cat $transcipt_list | xargs -P $threads -I % bash $pathtorepo/subsetgff.sh % $gff
cat $transcipt_list | xargs -P $threads -I % Rscript $pathtorepo/PDB_to_bp.R %

grep -v RefAllele pLDDT/* > Protein_popstats.txt
grep -v RefAllele prot_files/* > Protein_pLDDT.txt
