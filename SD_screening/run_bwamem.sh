sample=$1
reference=$2
prefix=$3
dir="/scratch/bnavarr2/analysis/sd_in2Rmal_dpgp3/mapping_dpgp3_to_chr2R/bams"


if [ "$#" -ne 2 ]; then
    echo "Usage: $0 sample_name reference out_prefix"
	exit 1
fi


####Build bwa index

if [ -f "$reference.bwt" ]; then
        echo "$reference.bwt exists, moving on..."
else
        echo "Prepping bwa index..."
        bwa index $reference
fi

if [ -s "$sample.sorted.bam.bai" ]; then
        echo "$sample".sorted.bam" already exists: exiting"
	exit 1
fi


echo "Working on sample: $sample"
	
#echo "
###define header
#RG="@RG\tSM:$sample\tPL:illumina\tLB:lib1\tPU:unit1"
RG="@RG\tID:$prefix\tSM:$sample\tPL:illumina\tLB:lib1\tPU:unit1"

###map reads

if [[ ! -e "$sample.sam" ]]; then
	echo "mapping $sample to $reference" ; 
	bwa mem -t 12 -M -R $RG $reference $sample"_pass_1.fastq.gz" $sample"_pass_2.fastq.gz" > $sample".sam"
	else
	echo "mapping already done, moving on..."
fi &&

if [ -s "$sample.sam" ]; then
	echo "$sample.sam done... sorting and converting to bam... "
        
	#sort sam and convert to bam
	samtools view -bS $sample".sam" | samtools sort -@12 - -o $sample".sorted.bam" && 
	samtools index $sample".sorted.bam" 
	else
		echo "$sample".sam" not created"
		exit 1
	fi
if [ -s "$sample.sorted.bam.bai" ]; then 
	echo "$sample".sorted.bam" generated correctly, All done!" ;
	rm $sample".sam";
	else
                echo "$sample".sorted.bam" not created"
                exit 1
fi



#"
