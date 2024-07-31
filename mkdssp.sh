zcat $2/${1}.pdb.gz | grep ATOM   > $2/${1}.pdb
mkdssp -i $2/${1}.pdb | grep -A 100000 "RESIDUE AA" | cut -c 1-5,13-14,35-39 | grep -v "\\!" | awk '{print $1,$2,$3}' > dssp/${1}.dssp
sed -i -E 's/ ([A])([0-9])/ \1 \2/g' $2/${1}.pdb
sed -i -E 's/([0-9]+)-([0-9]+)/\1 -\2/g' $2/${1}.pdb
