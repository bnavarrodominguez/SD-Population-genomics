# SNP calling and annotation

## GATK
Depends: samtools, picard, bwa mem and GATK
Create a list of readfile prefixes with your sample IDs to give it on the command line: 

e.g. reads.txt
sample1
sample2

where in the same directory you have sample1_1.fastq and sample1_2.fastq

```
$ sh GATK.bestpract.pipeline.AML.R1.sh reads.txt reference.fasta
```

## SNPeff

Dependencies: SNPeff
```
$ snpeff_snpsift.sh path_to_SNPeff input.vcf out_dir
```

## Keep only biallelic & polymorphic SNPs

 Requires: bcftools, samtools


```
$  select_var_bial.sh input.vcf output.vcf
```

## Filter heterozygous calls

Requires: GATK

```
$ filterHetCalls.sh input.vcf.gz
```

## Intersect SNPs 

 Requires: bcftools. 
 VCF files must be tabix indexed (`tabix -f file.vcf.gz`). 

```
$ intersect_snps.sh file1.vcf.gz file2.vcf.gz file3.vcf.gz 
```

# Population genomics

Simplify VCF (for each sample)

```
$ SimplifyVCF_Basic.readFiltfromVCF.pl population_1.simple.vcf
$ SimplifyVCF_Basic.readFiltfromVCF.pl population_2.simple.vcf
$ SimplifyVCF_Basic.readFiltfromVCF.pl population_3.simple.vcf
```

Run population genomics analyses

```
$ Windows_Basic.pl sample_vcf_list sample_name_list gff_file min_TajD_depth min_Fst_depth max_depth window_size outfile
```
1. sample_vcf_list: comma-separated list of (simplified) vcfs for each sample
2. sample_name_list: comma-separated list of names for each sample
3. gff_file: gff-formatted list of masked sites; excluded *completely* from analysis (if nothing needs to be masked, give an empty text file instead)
4. min_tajD_depth:  minimum sample depth required to include a variant in Tajima's D calculations (required for focal sample)
5. min_Fst_depth minimum sample depth required to include a site in Fst calculations (required for all samples)
6. max_depth: maximum possible sample depth across all samples
7. window_size: size (in bp) of non-overlapping windows to consider
8. outfile: name of the output file


Example:

```
$ Windows_Basic.v2.2.pl population_1.simple.vcf,population_2.simple.vcf,population_3.simple.vcf pop1,pop2,pop3 mask_sites.gff 8 8 1000 10000 out.txt
```

# Simulations

## ABC (infer time of the sweep)

Depends: ms (http://home.uchicago.edu/rhudson1/source/mksamples.html)

Script to generate *m* posterior samples under a sweep model (absolute bottleneck of N=1 at a time t 4N<sub>e</sub> generations). 

We simulated with values of S <sub>sim</sub> drawn from a uniform distribution ±5% of S <sub>obs</sub> .We considered a prior uniform distribution of time of the sweep (t) ranging from 0 to 4N e generations. 

The rejection sampling algorithm is as follows: 
1.  Draw S <sub>sim</sub and t from prior distributions; 
2. Simulate 1000 samples using the coalescent under a selective sweep model
3. Calculate average summary statistics for drawn S <sub>sim</sub> and t
4. Accept or reject chosen parameter values conditional on |π <sub>obs</sub> − π <sub>sim</sub> | ≤ ε, |D <sub>obs</sub> − D <sub>sim</sub> | ≤ ε;
5. Return to step 1 and continue simulations until m desired samples from the joint posterior probability distribution are collected.

```
$ run_ms.py m out_prefix sumstats

```
1. m: Number of posterior samples
2. Out prefix
3. sumstats: file that contains π <sub>obs</sub>, S <sub>obs</sub> and D <sub>obs</sub> (separated by spaces)

## Model fitting

Generate m samples under three models:
- a sweep (absolute bottleneck of N = 1) at a time t
- constant population size
- exponential growth

```
$ run_ms.model_fitting.py Nsims OutPrefix ObsS tsweep alpha
```

1. Nsims: number of simulations
2. Out prefix
3. ObsS: observed S
4. tsweep: time of the sweep, estimated by abc
5. Alpha: exponential growth factor

# TE calling

## Generate TSV file for McClinctock

```
$ convert_gff_and_tsv_uniqIDs.sh gff fasta
```

1. gff: from RepeatMasker output
2. TE.fasta: fasta file with TE sequences

## Run McClinctock

```
$ run_ngstemapper_telocate_temp.sh
$ run_ngstemapper_temp.sh 
$ run_popoolation.sh
$ run_retroseq.sh
```
## Merge close TE insertions

Depends: [TE_insertion_merger.py](https://github.com/KamilSJaron/reproductive_mode_TE_dynamics/blob/master/empirics/TE_insertion_merger.py) from KamilSJaron

Concatenate bed files from McClinctock and merge overlapping evidence within a distance smaller than *d* bp

```
$ cat *.bed > all_te_calls.bed

## Format as a table for TE_insertion_merger.py
$ get_tab_from_bed.sh all_te_calls.bed

$ TE_insertion_merger.py all_te_calls.tab3  d
```

# Detect inversions from Illumina data

Read mapping (repeats in the reference must be masked)

```
$ run_bwa_mem.sh sample_name reference out_prefix
```

Detect In(2L)t, In(2R)NS, In(2R)Mal and Sd-RanGAP in the bams. 

```
$ detect_inversions.sh bamfile
$ detect_Sd.sh bamfile
$ Rscript SD_screening.R
```


# Recombination

## Linkage disequilibrium decay with distance

```
$ run_plink.sh vcf_file outdir chr coord_start coord_end thin-count(=0 for no thinning)
```
1. vcf_file with only only non-singletons (mac>1) and biallelic SNPs
2. outdir: output directory
3. chr: chromosome
4. coord_start: start coordinate
5. coord_end: end coordinate
6. thin-count: number of SNPs to reduce data, use 0 for no thinning (use all SNPs)

## Generate input for RecMin from VCF

```
$ SimplifyVCF_Basic.readFiltfromVCF.pl sample.vcf
$ Rscript recmin.R sample.simple.vcf
```

## Runs of shared and private SNPs

```
$ Rscript shared-priv_distribution.R private_sites.txt shared_sites.txt
```
