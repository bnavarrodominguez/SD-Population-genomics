#!/bin/bash 
#####Beatriz Navarro Dominguez
#####Generate TSV file for McClintock from RepeatMasker GFF

gff=$1
te_fasta=$2

# Change the name column so it's the same than in the example data of McClintock
sed "s/Target=/ID=/g" $gff > tmp

# Change the second column so it's the same than in the example data of McClintock

sed -i "s/dispersed_repeat/transposable_element/g" tmp

#Drop the last columns and change the 6 for a dot as in the example data
awk {'print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t\.\t"$7"\t"$8"\t"$9"__"NR'} tmp > $(basename $gff .gff).mcclintock.gff 

# Remove commented lines

sed -i -e '/^[ \t]*#/d' $(basename $gff .gff).mcclintock.gff 

#Remove the "-int" labels added by repeatmasker? while generating the GFF

sed -i 's/\-int//' $(basename $gff .gff).mcclintock.gff

# Generate tsv file from the Rmasker database

#grep ">" specieslib_mod2_centromere_v2.fasta | sed "s/>//" | sed "s/#/\t/" | sed "s/\//_/" > dmel.chromosomes.fa.TE.mcClintock.tsv

awk {'print $9'} $(basename $gff .gff).mcclintock.gff |sed "s/ID=//" > tmp

awk {'print $9'} $(basename $gff .gff).mcclintock.gff |sed "s/ID=//" | sed "s/__.*$//"  >  tmp2

paste tmp tmp2 > $(basename $gff .gff).mcclintock.tsv

#rm tmp*

# Clean up fasta headers (from # until the end of the line)
sed 's/\#.*//' $te_fasta > $(basename $te_fasta .fa).mcclintock.fa

