head -n 1 $2 | awk -v OFS="\t" '{print $1,$2,$3,NA,NA,NA,NA,NA,NA}' > Freq/${1}.frq
grep $1 $2 > Freq/${1}.frq2
cat Freq/${1}.frq2 >> Freq/${1}.frq
rm Freq/${1}.frq2
