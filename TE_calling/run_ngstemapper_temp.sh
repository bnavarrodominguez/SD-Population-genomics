#!/bin/bash 
#####Beatriz Navarro Dominguez
#####Run ngs_te_mapper and TEMP through McClintock


#!/bin/bash -l

run_dir=/scratch/bnavarr2/analysis/sd_te/mcclintock_running
bin_dir=/scratch/bnavarr2/scripts/mcclintock
reference=dmel.chromosomes.fa
te_fasta=specieslib_mcClintock.fa
gff_file=dmel.chromosomes.fa.TE.mcClintock.gff
tsv_file=dmel.chromosomes.fa.TE.mcClintock.tsv
reads_1=$1
reads_2=$2
threads=8
out_dir=$run_dir/$(basename $reads_1 _1.fastq)_temp

# Run the pipeline 
cd $run_dir

bash $bin_dir/mcclintock.sh -i -m "TEMP" -r $run_dir/$reference -c $run_dir/$te_fasta -g $run_dir/$gff_file -t $run_dir/$tsv_file -1 $run_dir/$reads_1 -2 $run_dir/$reads_2 -p $threads -o $out_dir

