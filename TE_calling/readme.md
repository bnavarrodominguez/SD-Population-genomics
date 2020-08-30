# Transposable element calling 

## Dependencies
- McClintock & its dependencies: https://github.com/bergmanlab/mcclintock
- TE\_insertion\_merger.py: from https://github.com/KamilSJaron/reproductive_mode_TE_dynamics (last accessed Jan 31 2019)

## Use

### 0. Data needed:
- Reference genome (fasta): `reference.fa`
- Transposable element (TE) sequences (fasta): `te.fa`
- Annotation of TEs in the reference genome (GFF): `te.gff`
- Paired read libraries: `reads_1.fastq` `reads_2.fastq`

### 1. Edit GFF and TE fasta headers & generate TSV file for McClintock

```
sh convert_gff_and_tsv_uniqIDs.sh te.gff te.fa
```

### 2. Run McClintock
Edit the scripts to define run directory, directory  where McClontock is installed, reference.fa, te.gff, te.tsv. 

```
sh run_retroseq.sh reads_1.fastq reads_2.fastq &> retroseq.log
sh run_popoolation.sh reads_1.fastq reads_2.fastq &> popoolation.log
sh run_ngstemapper_telocate_temp.sh reads_1.fastq reads_2.fastq &> ngstemapper_telocate_temp.log
```

### 3.  Concatenate output bed files and convert to table format
Assumes that the name column in the bed files is formatted like this:
`TEname_reference_pop_genotype_samplename_TEcaller_sr_##`
```
#All beds in a folder named "beds_to_merge" - bed from McClinctock output folders
cat beds_to_merge/*.bed > merged.bed
#Convert to tab format for TE_insertion_merger.py
sh get_tab_from_bed.sh merged.bed
```

### 4. Merge overlapping evidence within N bp
Script from https://github.com/KamilSJaron/reproductive_mode_TE_dynamics (last accessed Jan 31 2019)

```
TE_insertion_merger.py merged.tab 500
```

### 5. Filter calls that have been predicted only by one of the TE callers


