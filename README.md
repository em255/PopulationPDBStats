# PopulationPDBStats

This repo provides a pipeline that generates a genomewide analysis of protein structure along all available proteins.  It is recommended to provide a population vcf that has been annotated by snpEFF to associate mutations to protein structure.  


## Required Inputs

1. File containing primary transcripts of your genome or transcripts of interest
2. Directory with gzipped PDB (Folded Protein) files whose names correspond to the transcripts
3. gff file corresponding to transcripts with annotated CDS regions
4. Fasta with CDS sequences with names matching transcript names
5. (Optional) snpEFF annotated VCF file

## Running the Pipeline

1. Install conda environments:
   ```bash
   conda env create -f PDB_Tools.yml
   conda env create -f dssp.yml
   ```
(If dssp is giving you issues it may be required to install with some channel hierarchy with this command: `conda create -n dssp --strict-channel-priority -c anaconda -c salilab 'libboost=1.73.0' dss`)

2. Edit paths in "Process_PDB.sh" to give paths to input gff, transcript list, cds fasta (and snpEFF annotated vcf if applicable).
3. Run Pipeline
```bash
   sbatch Process_PDB.sh
```
## Running Test data
1. Test data is set to run with the master script "Process_PDB.sh" from within the repo.  It has sbatch header for SLURM environments.
```bash
sbatch Process_PDB.sh
```
## Output Tables
Output Table(s) show positions in genes with their protein structure information.  Column names in "Protein_pLDDT.txt.gz" refer to:

```bash
Chr - Chromosome 
Pos - Position in chromosome
transcript - Transcript associated to the protein
RefAllele - The erference Allele at that position
AA_Pos - The amino acid position associated with that nucleotide
CodonPos - The codon Position associated with that nucleotide
AminoAcid - Amino Acid identity
mean_pLDDT - The average "confidence" of all atoms contained within the amino acid
Gene_pLDDT - The average "confidnece" along all amino acids over the entire protein
AminoAcid_1let - 1 Letter abbreviation of amino acid
ASA - Available surface area
maxASA - Maximum available surface area for a given amino acid              
rASA - relative available surface area (ASA/maxASA)
FourFoldDegenerate - Does the nucleotide belong to a four-fold degenerate codon position (no single nucleotide polymorphism will change amino acid)
AA_ProportionPosition - Where does the amino acid lie proportionally along protein
```

The table "Protein_popstats.txt.gz" contain columns:
```bash
alleles - Allele of the mutation, multi-allelic snps are represented as multiple rows in the table
freqs - frequency of that allele
MAF - the overall minor allele frequency of all mutations at this position
Effect - snpEFF effect class (LOW:synonymous;MODERATE:nonsynonymous,in frame indels; HIGH:Frameshift, early stop, etc)
Alt_AA - New amino Acid if it is strictly a nonsynonymous mutation
Ref_FunctionalClass - Functional class of reference amino acid
Alt_FunctionalClass - Functional class of alternate amino acid
```

## Analyze Output Tables
1. Run R Markdown analysis and set working directory to director containing output tables "Protein_pLDDT.txt.gz" and "Protein_popstats.txt.gz" (if using vcf)
2. R Markdown generates plots that associate protein structure with mutations.  

