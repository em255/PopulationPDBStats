#Get a CDS fasta for each transcript
samtools faidx Ref/IRGSP-1.0_cds_2023-03-15.fasta $1 > CDS/${1}.cds.fasta
