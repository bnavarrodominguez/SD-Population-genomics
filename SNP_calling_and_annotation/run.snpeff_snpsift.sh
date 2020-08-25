#!/bin/bash 
#####Beatriz Navarro Dominguez
#####Run SNPeff/SNPsift in VCF files 

#######Load modules
module load java
module load samtools

######My variables
path=$1
vcf=$2
outdir=$3


######Usage
display_usage() {
	echo -e "Usage: $0 snpEff_path vcf out_dir \n"
	}

# if wrong number of arguments supplied, display usage 
	if [  "$#" -ne 3 ] 
	then 
		display_usage
		exit 1
	fi 
 
# check whether user had supplied -h or --help . If yes display usage 
	if [[ ( "$#" == "--help") ||  $# == "-h" ]] 
	then 
		display_usage
		exit 0
	fi 


###### Create out dir

mkdir -p $outdir;

##### Annotate VCF with snpEff & stat summary in csv file. Consider only canonical isoforms. 

java -Xmx4g -jar $path/snpEff.jar -canon -c $path/snpEff.config dmel_r6.12 -csvStats $outdir/$(basename $vcf .vcf).stats.canon.csv -v $vcf > $outdir/$(basename $vcf .vcf).canon.snpeff.vcf

#### Move summary to outdir
mv snpEff_summary.html $outdir

##### Parse anno
cat $outdir/$(basename $vcf .vcf).canon.snpeff.vcf | perl $path/scripts/vcfEffOnePerLine.pl | java -jar $path/SnpSift.jar extractFields -e "\""."\"" - CHROM POS ID REF ALT AF DP MQ "EFF[*].EFFECT" "EFF[*].IMPACT" "EFF[*].FUNCLASS" "EFF[*].GENE" "EFF[*].BIOTYPE" "EFF[*].CODING" "ANN[*].GENEID" > $outdir/$(basename $vcf .vcf).eff_opl.tab

bgzip -c $outdir/$(basename $vcf .vcf).canon.snpeff.vcf > $outdir/$(basename $vcf .vcf).canon.snpeff.vcf.gz && rm $outdir/$(basename $vcf .vcf).canon.snpeff.vcf
