#Get a CDS fasta for each transcript
samtools faidx $2 $1 > CDS/${1}.cds.fasta
