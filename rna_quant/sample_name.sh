ls $1 |sed 's|quant/||g' |sed 's/.genes.results//g'|paste -s -d "\t"|awk '{print "Gene""\t"$0}' > $2
