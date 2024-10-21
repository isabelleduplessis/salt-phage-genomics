# Introduction

The goal of this research project is to characterize the viral populations found in S. alterniflora plants and surrounding soil in a salt marsh wetland environment. 

# Files

host_taxonomy/genera: Contains lists of microbial host genera known to possess biogeochemical properties in iron oxidaztion, nitrification, sulfate reduction, and sulfur oxidation.

host_taxonomy/getorders.sh: Takes taxonomic genera and uses fulllineage.dmp from NCBI taxonomy to find corresponding taxonomic orders.

saturation_nmds_plot.R: Used to produce saturation plot and NMDS plot.

abundance.R: Used to produce viral and prokaryotic abundance plots and vOTU count and rank order plots.

host_plot.R: Used to produce host phyla plot and biogeochemical properties.

## Intermediate Files for Plotting

Intermediate files required for R scripts found in intermediate_files:

votus_cov75thres.txt: Contains coverage, depth, and rpkm data for each viral Operational Taxonomic Unit (vOTU).

metadata.csv: Contains metadata for the 24 samples collected including compartment and phenotype. This data was previously provided to us.

qPCR_data.csv: Contains qPCR data used to measure absolute bacterial abundance. This data was previously provided to us.

HostPrediction.tsv: Contains results from iPHoP host prediction and corresponding taxonomy

hostproperties.tsv: Indicates biogeochemical properties matched to each host order. With column headers being "virus", "hostorder", "sulfuroxidizing", "sulfatereducing", "ironoxidizing", "nitrifying."


## Packages Used:

The following packages are required to run R scripts:

ggplot2 3.4.2

reshape2 1.4.4

vegan 2.5-6

grid 4.0.2

gridExtra 2.3

readr 2.1.4

ggpubr 0.6.0

ggvenn 0.1.10


# Methods

The following workflow was used to produce the intermediate files used for plotting.

## Packages Used:

BWA v0.7.17

Blast-Plus v2.13.0

CheckV v1.0.1

Virsorter2 v2.2.4

Samtools v1.14

CD-HIT v4.8.1

iPHoP 

Prodigal

HMMER

## Mapping raw reads to PIGEON and GOV 2.0 Database

The PIGEON Database was retrieved from:
https://datadryad.org/stash/dataset/doi:10.25338 and indexed using bwa v0.7.17:

    bwa index -p PIGEONv1.0 PIGEONv1.0.fa.gz

For each pair of raw read fastq files:

the following command was used to map the raw read to the PIGEON
database:


    bwa mem -t 4 PIGEONv1.0 R1_${samplenumber}.fastq R2_${samplenumber}.fastq > sam/${samplenumber}.sam

Next, to prepare the files for samtools coverage, the output sam files
were sorted and output as bam files using the following command for each
raw read sam file:

    samtools sort -@4 sam/${samplenumber}.sam -o ${samplenumber}.bam

The coverage for each read was then obtained using the following command
on each bam file:

    samtools coverage ${samplenumber}.bam -o rawreads_coverage/${samplenumber}_coverage

The output was a file containing columns for rname, startpos, endpost,
numreads, covbases, coverage, meandepth, meanbaseq, meanmapq. We were
interested in the coverage value, and used the following commands to
obtain sequences, that had \> 50, \> 60, \> 70, and \> 75 percent
coverage:

    for FILENAME in rawreads_coverage/*; do awk '{ if($6>50){print $0 "\t" FILENAME}}' $FILENAME \
    >> coverage50.txt; done
    for FILENAME in rawreads_coverage/*; do awk '{ if($6>60){print $0 "\t" FILENAME}}' $FILENAME \
    >> coverage60.txt; done
    for FILENAME in rawreads_coverage/*; do awk '{ if($6>70){print $0 "\t" FILENAME}}' $FILENAME \
    >> coverage70.txt; done
    for FILENAME in rawreads_coverage/*; do awk '{ if($6>75){print $0 "\t" FILENAME}}' $FILENAME \
    >> coverage75.txt; done


This was repeated for the GOV 2.0 database.

7 reads mapped to GOV with 75% coverage. Of those 7 reads, 4 of them
were unique sequences. There were some GOV sequences that multiple raw
reads mapped to.


2 reads mapped to PIGEON database with 75% coverage. 1 of the 2 PIGEON sequences matched to 1 of the 4 GOV 5KB sequences.

Therefore, we have 5 total unique sequences from raw reads that mapped to PIGEON and GOV 5KB.

## Gene Prediction using Prodigal

The file viral_virsorter_cdhit.fna contains sequences predicted to be viral using VirSorter2.

    prodigal -i viral_virsorter_cdhit.fna -o viral_virsorter_cdhit_genes -a viral_virsorter_cdhit_proteins

## Gene Annotation with Blast

Blast-Plus version 2.13.0 was used for gene annotation of the output
\*.gff file from prodigal. This was done using UniProt available dataset swissprot.fa from [swissprot](# https://www.uniprot.org/help/download). This was used as the database for the following Blast commands.

The Uniprot database was obtained using:

    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

The Uniprot database was indexed using:

    gunzip uniprot_sprot.fasta.gz
    makeblastdb -in uniprot_sprot.fasta -dbtype prot -parse_seqids -out uniprot_sprot -title "uniprot_sprot"

Next, the output from gene prediction was annotated using blast against
the Uniprot database:

    blastp -query viral_virsorter_cdhit_proteins -db uniprot_sprot -num_threads 4 -out viral_virsorter_uniprot_blast


## Increasing Sequences

The VirSorter prediction sequences can be found at:

    viral_virsorter_cdhit.fna

First, the VirSorter prediction file was indexed using bwa to prepare the database for mapping:

    bwa index -p viral_virsorter viral_virsorter_cdhit.fna

### Mapping

Next, the trimmed reads were each mapped to the database:

    bwa mem -t 4 viral_virsorter trimmed/$var1\_L002_R1_001.trimmed.fastq trimmed/$var1\_L002_R2_001.trimmed.fastq > $var1\_L002_aln-pe_001.trimmed.sam

where "$var1" is the prefix of the file.

Next, samtools was used to convert sam files to sorted bam files:

    for FILENAME in *.sam; do
        samtools sort -@4 $FILENAME -o bam/$FILENAME\.bam

To filter bam files to only include the reads that mapped:

    for FILENAME in *sam.bam; do
        samtools view -b -F 4 $FILENAME > mapped_$FILENAME; done

To convert the mapped reads to fastq files, the following commands were used:

    for FILENAME in mapped*; do
        samtools sort -n -o sort_$FILENAME $FILENAME
        bedtools bamtofastq -i $sort_FILENAME -fq r1_$FILENAME -fq2 r2_$FILENAME; done

### MEGAHIT

MEGAHIT was run on each pair of fastq files to reassemble them:

    megahit -1 r1_mapped_${samplenumber}.fq -2 r2_mapped_${samplenumber}.fq -o megahitoutput${samplenumber} -t 4

### VirSorter

Viral sequences from each assembly were extracted used VirSorter2:

    virsorter run -w ${samplenumber}.out -i megahitoutput${samplenumber}/final.contigs.fa --min-length 1500 -j 4 all

The following command was used to keep the sequences above 5KB in length:

    awk -v RS='>[^\n]+\n' 'length() >= 5000 {printf "%s", prt $0} {prt = RT}' ${samplenumber}.out/final-viral-combined.fa > ${samplenumber}_viral.fa


### CheckV

Set up CheckV on personal computer:

    conda create -n checkv -y
    conda activate checkv
    conda conda install -c conda-forge -c bioconda checkv
    checkv download_database ./
    cd checkv-db-v1.5
    export CHECKVDB=checkv-db-v1.5

Run CheckV:

    for filename in viral/*.fa; do checkv end_to_end $filename $filename\results -t 2; done

Get non-\"Not determined\" sequences, where denovo_mapped_checkv.fna contains the concatenated VirSorter outputs from the previous step:

    for filename in *.fa; do 
        grep -v "Not-determined" quality_summary.tsv > quality_summary_filtered.tsv;
        grep -A 1 -f <(awk '(NR>1) {print $1}' quality_summary_filtered.tsv) viral/${filename} >> \
    viral/$filename\.toappend.fna
    done
    cat <files> > append.fna
    sed -i '/--/d' append.fna
    cat denovo_mapped_checkv.fna toappend/append.fna > denovo_mapped_checkv_new.fna

### CDHIT

    module load cdhit/4.8.1-n4smhr
    cd-hit-est -i denovo_mapped_checkv_new.fna -o saltviralseqs.fna -c 0.95 -s 0.8

This resulted in 769 final viral sequences.

## Abundance & Diversity Analysis

### BWA index

First, the viral sequences were indexed using BWA:

    bwa index -p saltviralseqs saltviralseqs.fna

### BWA mem

    bwa mem -t 4 saltviralseqs rawreads/JK03-${samplenumber}_*R1*.fastq rawreads/JK03-${samplenumber}_*R2*.fastq > sample_${samplenumber}.sam

### Samtools Coverage

    samtools view -bS sample_${samplenumber}.sam > sample_${samplenumber}.bam
    samtools sort sample_${samplenumber}.bam > sample_${samplenumber}_sorted.bam
    samtools depth sample_${samplenumber}_sorted.bam > sample_${samplenumber}_depth.txt
    samtools coverage sample_${samplenumber}_sorted.bam > sample_${samplenumber}_coverage.txt

### RPKM and Mean Depth

    for i in {1..24}; do totreads=$(samtools view -c sample_${i}_sorted.bam); echo "sample_${i},$totreads"; done > totalreads.csv

    for ARG in {1..24}; do totreads=$(grep "sample_$ARG," totalreads.csv | cut -f2 -d ","); awk -F'\t' 'NR==1 || $6>=75' "sample_${ARG}_coverage.txt" | grep -v "^#" | sed -e "s/^/sample_$ARG\t/" | awk -v t=$totreads '{print $1"\t"$2"\t"$7"\t"$8"\t"($5*10**9)/($4*t)}'; done > votus_cov75thres.txt;
    sed -i '1 i\#sample\tvotus\tcoverage\tmeandepth\trpkm' votus_cov75thres.txt 

    for ARG in {1..24}; do samtools view -b -F 4 sample_${ARG}_sorted.bam > sample_${ARG}_mapped.bam; done
    for i in {1..24}; do totreads=$(samtools view -c sample_${i}_mapped.bam); echo -e "sample_${i}\t$totreads"; done > totalreads_mapped.csv

## Host Prediction with iPhop

Separate the file of 769 sequences into several smaller ones to run
iphop:

    awk 'BEGIN {n_seq=0;} 
    /^>/ {if(n_seq%100==0){file=sprintf("seq%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < seq.fna

The output is 8 files with up to 100 sequences. For each of the output
files, the following command was run:

    iphop predict --fa_file seq0.fa --db_dir iphop_db/Sept_2021_pub/ --out_dir iphop_output0

The results of each file were concatenated. The following command output the number of sequences that belonged to each unique host taxonomy:

    awk -F ',' '{print $3}' host_prediction_to_genus_m90_all.csv | sort | uniq -c | sort -nr > hostfreq_genus.txt

55 unique hosts were identified at the genus level, and 64 out of 769 viral sequences were matched to a host genus.


## Relevant Gene Annotation - searching for auxiliary metabolic genes

Gene Prediction with Prodigal

    prodigal -i saltviralseqs.fna -o saltviralgenes -a saltviralproteins

Annotation with HMMER

    while read -r line; do sbatch --export=id=$line hmmr2.sbatch; done < relevant_gene_annotation/relevant_metabolic_genes.txt

hmmr2.sbatch:

    hmmsearch --domtblout $id\.domtblout --tblout $id\.tblout kofam_profiles/$id\.hmm saltviralproteins.faa > $id\1.out

Format output

    ls *domtblout | while read file; do relevant_gene_annotation/parse_hmmr_domtblout.sh $file; done 

Get significant hits

    cat *domtblout.tab | awk -F'\t' '$7<0.001' |  awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"($17-$16+1)/($6)}' | awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5]+=$6}END{for(i in a) print i"\t"a[i]}' | awk -F'\t' '$6>=0.5' > significant_hits

This resulted in 3 significant hits.
