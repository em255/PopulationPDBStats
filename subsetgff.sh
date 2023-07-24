grep $1 $2 | awk -v OFS="\t" '{if($3=="CDS"){print}}'  > gff/${1}.gff
