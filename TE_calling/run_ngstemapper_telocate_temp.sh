#!/bin/bash 
#####Beatriz Navarro Dominguez
#####Run ngs_te_mapper, TEMP and TE-locate through McClintock

#!/bin/bash -l

run_dir=/home/user/mcclintock_running
bin_dir=/bin/mcclintock/
reference=reference.fa
te_fasta=te.mcclintock.fa
gff_file=te.mcclintock.gff
tsv_file=te.mcclintock.tsv
reads_1=$1
reads_2=$2
out_dir=$run_dir/$(basename $reads_1 _1.fastq)_ngs_temp_telocate
threads=8

# Run the pipeline 
cd $run_dir

bash $bin_dir/mcclintock.sh -i -m "ngs_te_mapper TEMP TE-locate" -r $run_dir/$reference -c $run_dir/$te_fasta -g $run_dir/$gff_file -t $run_dir/$tsv_file -1 $run_dir/$reads_1 -2 $run_dir/$reads_2 -p $threads -o $out_dir

