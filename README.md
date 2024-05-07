# PopulationPDBStats

Correlate population allele frequencies and functional mutation effects with protein structure metrics.

## Required Inputs

1. File containing primary transcripts or transcripts of interest
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
