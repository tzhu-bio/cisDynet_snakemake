awk -vFS="\t" '{if(FNR == 1){NN++;} a[FNR]=(a[FNR]?a[FNR]FS:"")(NN<2?$1 FS $6:$6)} END {for(i=1;i<=FNR;i++) print a[i]}' $1 | sed '1d' > $2
