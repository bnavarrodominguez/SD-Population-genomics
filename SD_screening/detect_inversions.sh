bam=$1
sample=$(basename $bam .sorted.bam)

#samtools view -H $bam > $sample".head"
samtools view -f 1 -F 14 -q 30 $bam | awk '{if($9>3000000 && $3==1 || $9>3000000 && $3==2) print $0 }' | awk '{print $3"\t"$4"\t"$4+1"\t"$7"\t"$8"\t"$8+1}' > $sample".sel.bedpe"

if [ -s "$sample.sel.bedpe" ]; then
        echo "$sample".sel.bedpe" generated correctly, All done!" ;
        else
                echo "$sample".sel.bedpe" not created"
                exit 1
fi

