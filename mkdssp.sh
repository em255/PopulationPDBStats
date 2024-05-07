zcat $2/${1}.pdb.gz | grep ATOM > $2/${1}.pdb
mkdssp -i $2/${1}.pdb | grep -A 100000 "RESIDUE AA" | cut -c 1-5,13-14,35-39 | grep -v "\\!" | awk '{print $1,$2,$3}' > dssp/${1}.dssp
