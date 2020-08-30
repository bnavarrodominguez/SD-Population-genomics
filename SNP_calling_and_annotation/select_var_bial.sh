#!/bin/bash 
#####Beatriz Navarro Dominguez
#####Keep only biallelic & polimorphic SNPs

######My variables

input=$1
output=$2


#####Usage

display_usage() {
	echo "Usage: $0 input.vcf output.vcf \n Requires: bcftools, samtools\n"
	}

# if wrong number of arguments supplied, display usage 
	if [  "$#" -ne 2 ] 
	then 
		display_usage
		exit 1
	fi 
 

####Keep only biallelic & polimorphic SNPs

bcftools view -v snps -m2 -M2 --min-ac 1:minor $input | awk '/^#/||$7=="PASS"' | awk '$1 ~ /^#/ {print $0;next} {if ($4 ~ /A|C|T|G/ && $5 ~ /A|C|T|G/) print $0}' | bgzip  >   $output


