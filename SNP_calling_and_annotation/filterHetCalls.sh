#!/bin/bash 
#####Beatriz Navarro Dominguez
#####Filter heterozygous calls in a VCF file

#####Usage

display_usage() {
	echo "Usage: $0 file.vcf.gz \n Requires:GATK\n"
	}

# if wrong number of arguments supplied, display usage 
	if [  "$#" -ne 1 ] 
	then 
		display_usage
		exit 1
	fi 
 



######My variables

vcf=$1
vcfVF=$(basename $vcf).VF.vcf.gz
vcfSV=$(basename $vcf .vcf.gz).SV.vcf.gz


#####Label heterozygous calls 

gatk VariantFiltration \
-V $vcf \
-O $vcfVF \
--genotype-filter-expression "isHet == 1" \
--genotype-filter-name "isHetFilter"

#### Convert heterozygous calls to N
gatk SelectVariants \
-V $vcfVF \
--set-filtered-gt-to-nocall \
-O $vcfSV



