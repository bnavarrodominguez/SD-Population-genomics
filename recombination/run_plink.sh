#!bin/bash
#run.plink.sh

if [ "$#" -ne 6 ]; then
    echo "$0 vcf_file outdir chr coord_start coord_end thin-count(=0 for no thinning)"
        exit 2
fi


vcf_file=$1
outdir=$2
chr=$3
coord_start=$4
coord_end=$5
thin=$6


mkdir $outdir ;
cd $outdir;
ln -s $vcf_file . ;


if (( "$thin" == 0 ))
then
        echo "No thinning..."
        plink --threads 16 --vcf $(basename $vcf_file) --set-hh-missing --mind 0.1 --geno 0.1 --biallelic-only --maf 0.2 --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 3 --chr $chr --from-bp $coord_start --to-bp $coord_end

else

      echo "Thinning to $thin SNPs"
        plink --threads 16 --vcf $(basename $vcf_file) --set-hh-missing --mind 0.1 --geno 0.1 --biallelic-only --maf 0.2 --thin-count $thin --r2 --ld-window-r2 0 --ld-window 999999 --ld-window-kb 3 --chr $chr --from-bp $coord_start --to-bp $coord_end


fi

