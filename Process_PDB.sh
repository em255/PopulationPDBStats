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


#Get summary tables
grep -v RefAllele pLDDT/* > Protein_popstats.txt
head -n 1  pLDDT/* > Protein_popstats.txt_header
cat Protein_popstats.txt_header Protein_popstats.txt > Protein_popstats.txt_temp
mv Protein_popstats.txt_temp Protein_popstats.txt
rm Protein_popstats.txt_header

#Get summary tables
grep -v RefAllele prot_files/* > Protein_pLDDT.txt
head -n 1  prot_files/* > Protein_pLDDT.txt_header
cat Protein_pLDDT.txt_header Protein_pLDDT.txt > Protein_pLDDT.txt_temp
mv Protein_pLDDT.txt_temp Protein_pLDDT.txt
rm Protein_pLDDT.txt_header
