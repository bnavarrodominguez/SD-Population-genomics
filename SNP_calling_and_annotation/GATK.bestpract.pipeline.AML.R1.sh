# Amanda Larracuente
#!/bin/bash
#
# To run in SLURM script, load samtools, picard, bwa mem and GATK modules
# Create a list of readfile prefixes with your sample IDs to give it on the command line: 
# e.g. List.txt
# sample1
# sample2
#
# where in the same directory you have sample1_1.fastq and sample1_2.fastq
#
#

#####	Software 
# Replace with paths on your system and appropriate memory requirements

function GATK {
    module load gatk/3.5
    java -Xmx60g -Djava.io.tmpdir=`pwd`/tmp -jar /home/alarracu/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar $@
    module unload gatk/3.5
}

function GATK_Filter {
    module load gatk/3.5
    java -Xmx60g -Djava.io.tmpdir=`pwd`/tmp -jar /home/alarracu/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
    --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' $@
    module unload gatk/3.5
}

function PICARD {
    module load java/1.9.126
    java -Xmx60g -Djava.io.tmpdir=`pwd`/tmp -jar /software/picard/2.0.1/dist/picard.jar $@
    module unload java/1.9.126
}

function SAMTOOLS {
    module load samtools/1.3.1
    samtools $@
    module unload samtools/1.3.1
}
function BWA {
    module load bwa/0.7.9a
    bwa $@
    module unload bwa/0.7.9a
}
#######	Files
# Reads in list of samples you want to process on the command line
# Set reference path

INFILE=$1
#REFERENCE="/home/alarracu/dmel-all-chromosome-r6.03.numerical.fasta"
REFERENCE="/home/alarracu/Drosophila_melanogaster.BDGP6.v88.dna.toplevel.NUM.fa"

echo "input $INFILE"
PREFIX=${INFILE//.list.txt/}
echo "prefix is $PREFIX"
########  Prep reference file

if [ -f "$REFERENCE.bwt" ]; then
        echo "$REFERENCE.bwt exists, moving on..."
else
	echo "Prepping bwa index..."
        BWA index $REFERENCE  
fi

if [ -f "$REFERENCE.dict" ]; then
	echo "$REFERENCE.dict exists, moving on..."
else 
	echo "Prepping reference genome..."
	
	# Create sequence dictionary using Picard Tools.
	# the following command produces a SAM-style header file describing the contents of our fasta file.

PICARD \
CreateSequenceDictionary \
REFERENCE=$REFERENCE \
OUTPUT="$REFERENCE.dict"

	echo "created sequence dictionary ${REFERENCE%.fa}.dict for the reference genome."
fi

########	Create the fasta index file.
if [ -f "$REFERENCE.fai" ]; then
	echo "$REFERENCE.fai exists, moving on..."
else
	echo "indexing the reference genome...."
	SAMTOOLS faidx $REFERENCE
	echo "Reference genome is now ready for GATK."
fi

########	Map and process reads
#Loop through all lines of your INFILE. This contains the list of sample names

while IFS='' read -r linename || [[ -n "$linename" ]]; do
    echo "Working on sample: $linename"
	RG="@RG\tID:$PREFIX\tSM:$linename\tPL:illumina\tLB:lib1\tPU:unit1"
	#map reads

	BWA mem -t 12 -M -R $RG $REFERENCE $linename"_1.fastq" $linename"_2.fastq" > $linename".sam"

	if [ -f "$linename.sam" ]; then
	#sort sam and convert to bam
PICARD SortSam \
INPUT=$linename.sam \
OUTPUT=$linename.sorted.bam \
SORT_ORDER=coordinate
echo "skip test"
	else 
		echo "$linename.sam not created"
		exit 1
	fi


	if [ -f "$linename.sorted.bam" ]; then
	#Mark duplicates in Picard
PICARD MarkDuplicates \
INPUT=$linename.sorted.bam \
OUTPUT=$linename.nodup.bam \
METRICS_FILE=$linename.dup.metrics.txt
echo "skip test"
	else 
		echo "$linename.sorted.bam not created"
		exit 1
	fi

	#index bam file in Picard
PICARD \
BuildBamIndex \
INPUT=$linename.nodup.bam 


if [ -f "$linename.nodup.bam" ]; then
	#analyze patterns of covariation in dataset
	#use SNP data from DPGP1 ensembl release 88
	#/home/alarracu/Drosophila_melanogaster.num.vcf has new contig names (numbered)  
	#get just 2L and 2R (1 and 2)
	
SAMTOOLS view -b -h -o $linename.2L2Rnodup.bam $linename.nodup.bam -L /scratch/alarracu_lab/incoming/reads/PS/targets.sam.list
	echo "begin recalibration..."  
#index bam file

PICARD \
BuildBamIndex \
INPUT=$linename.2L2Rnodup.bam

GATK \
-T BaseRecalibrator \
-R $REFERENCE \
-I $linename.2L2Rnodup.bam \
-nct 12 \
-L /scratch/alarracu_lab/incoming/reads/targets.list \
-knownSites /home/alarracu/Drosophila_melanogaster.num.sorted.vcf \
-o $linename.recal_data.table
 
else 
	echo "$linename.nodup.bam not created"
	exit 1
fi

	if [ -f "$linename.recal_data.table" ]; then
	#recalibrate based on first recalibration table
	echo "continue recalibration..."
GATK \
-T BaseRecalibrator \
-R $REFERENCE \
-I $linename.2L2Rnodup.bam \
-L /scratch/alarracu_lab/incoming/reads/targets.list \
-nct 12 \
-knownSites /home/alarracu/Drosophila_melanogaster.num.sorted.vcf \
-BQSR $linename.recal_data.table \
-o $linename.post.recal_data.table 

else 
	echo "$linename.recal_data.table not created"
	exit 1
fi

	echo "creating plots..."
	#make plots
	mkdir $PREFIX.recalibration.plots

GATK \
-T AnalyzeCovariates \
-R $REFERENCE \
-L /scratch/alarracu_lab/incoming/reads/targets.list \
-before $linename.recal_data.table \
-after $linename.post.recal_data.table \
-plots $PREFIX.recalibration.plots/$linename.recalibration_plots.pdf

if [ -f "$linename.post.recal_data.table" ]; then  
echo "Applying recalibration..."

GATK \
-T PrintReads \
-R $REFERENCE \
-I $linename.2L2Rnodup.bam \
-L /scratch/alarracu_lab/incoming/reads/targets.list \
-BQSR $linename.recal_data.table \
-o $linename.recal.bam

else 
	echo "$linename.post.recal_data.table not created"
	exit 1
fi    

if [ -f "$linename.recal.bam" ]; then  
echo "begin calling SNPs..."
	#Call variants, just on targets.list (e.g. 2R, 2L)

PICARD \
BuildBamIndex \
INPUT=$linename.recal.bam

GATK \
-T HaplotypeCaller \
-R $REFERENCE \
-I $linename.recal.bam \
-L /scratch/alarracu_lab/incoming/reads/targets.list \
--genotyping_mode DISCOVERY \
--emitRefConfidence GVCF \
--output_mode EMIT_ALL_CONFIDENT_SITES \
-ploidy 1 \
-stand_call_conf 30 \
-variant_index_type LINEAR \
-variant_index_parameter 128000 \
-o $linename.all_variants.vcf

else 
	echo "$linename.recal.bam not created"
	exit 1
fi    

done < "$INFILE"


############### 	Variant analysis ##############
echo "begin variant analysis..."

while IFS='' read -r line || [[ -n "$line" ]]; do
	echo "creating all variants list, added: $line"
	echo "$line.all_variants.vcf ">>$PREFIX.vcf.list

done < "$INFILE"

##Merge all
echo "merge vcf files..."

GATK \
-T GenotypeGVCFs \
-R $REFERENCE \
--variant $PREFIX.vcf.list \
-o $PREFIX.All_variants.vcf \
-allSites    

rm $PREFIX.vcf.list

if [ -f "$PREFIX.All_variants.vcf" ]; then  

echo "exclude indels from vcf files..."
GATK \
-T SelectVariants \
-nt 24 -V $PREFIX.All_variants.vcf \
--selectTypeToExclude INDEL \
-R $REFERENCE \
-o $PREFIX.All_SNPs.vcf

GATK \
-T SelectVariants \
-nt 24 -V $PREFIX.All_variants.vcf \
-selectType INDEL \
-R $REFERENCE \
-o $PREFIX.All_INDELSs.vcf
 
 else 
	echo "$PREFIX.All_variants.vcf not created"
	exit 1
fi	

echo "final SNP filtering"

GATK_Filter \
-T VariantFiltration \
-R $REFERENCE \
-V $PREFIX.All_SNPs.vcf \
--filterName "filter" \
-o $PREFIX.All_filtered_snps.vcf
 	
echo "Done"
