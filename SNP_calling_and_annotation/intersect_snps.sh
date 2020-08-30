#!/bin/bash 
#####Beatriz Navarro Dominguez
#####Intersect SNPs: Find shared & private SNPs within 3 VCF files

#####Variables
vcf1=$1
vcf2=$2
vcf3=$3


#### Usage
display_usage() {
	echo "Usage: $0 file1.vcf.gz file2.vcf.gz file3.vcf.gz \n Requires: bcftools. VCF files must be tabix indexed (tabix -f file.vcf.gz). \n"
	}

# if wrong number of arguments supplied, display usage 
	if [  "$#" -ne 3 ] 
	then 
		display_usage
		exit 1
	fi 
 


### Run

lab1=$(basename $vcf1 .vcf.gz)
lab2=$(basename $vcf2 .vcf.gz)
lab2=$(basename $vcf3 .vcf.gz)


bcftools isec -p isec_privates_"$lab1" -f "PASS" -C $vcf1 $vcf2 $vcf3

bcftools isec -p isec_privates_"$lab2" -f "PASS" -C $vcf2 $vcf3 $vcf1

bcftools isec -p isec_privates_"$lab3" -f "PASS" -C $vcf3 $vcf2 $vcf1

bcftools isec -p isec_sharedby2 -f "PASS" -n+2 $vcf3 $vcf2 $vcf1


