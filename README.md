# PopulationPDBStats
Correlate population allele frequencies and functional mutation effects with protein structure metrics.


# Required Inputs
1. snpEFF annotated VCF file
2. File containing primary transcripts or transcripts of interest
3. Directory with gzipped PDB (Folded Protein) files whose names correspond to the transcripts
4. gff file corresponding to transcripts with annotated CDS regions
5. Fasta with CDS sequences with names matching transcript names

# Running pipeline
-Edit paths in `Process_PDB.sh` to give paths to input vcf,gff, transcript list, cds fasta.
-Make sure all PDB files are in a PDB folder in your current directory
