#!/bin/bash
###Beatriz Navarro
###Script to convert a concatenated, sorted bed file to the table format required for TE_insertion_merger.py
###

if [ "$#" -ne 1 ]; then
    echo "$0 all_cat.bed"
else
        echo "Processing your results";

input=$1

cat $input | sed '/_non-reference_/ s/$/\tnon-reference/' | sed '/_reference_/ s/$/\treference/' > tmp.1
cat $input | awk {'print $4'} | sed 's/_non-reference/\t/' | sed 's/_reference/\t/' | awk {'print $1'} | paste tmp.1 - > tmp.2 
cat $input | awk {'print $4'} | sed 's/ngs_te_mapper/ngs-te-mapper/' | sed 's/_non-reference/\t/' | sed 's/_reference/\t/' | awk {'print $2'} |sed 's/_/\t/g' | paste tmp.2 - > tmp.3
sed -i 's/|/\t/g' tmp.3
mv tmp.3 $(basename $input .bed).tab
rm tmp.*

tab=$(basename $input .bed).tab

###Discard calls non-evidence of absence
grep -v 'nonab' $tab > $(basename $tab .tab).tab2

###Reorder the columns for TE_insertion_merger.py
cat $(basename $tab .tab).tab2 | awk {'print $1,$2,$3,$6,$9"_"$10,$11,$8,$7,$12,$13'} > $(basename $tab .tab).tab3


fi
