# salt-phage-genomics

# Introduction

# Files

host_taxonomy
uses fulllineage.dmp

figures

saturation_nmds_plot.R:\
Contains R script used to produce saturation plot and NMDS plot. Requires intermediate files votus_cov75thres.txt and metadata.csv.\
This script requires the following packages:\
        ggplot2 3.4.2\
        reshape2 1.4.4\
        vegan 2.5-6\
        grid 4.0.2\
        gridExtra 2.3

abundance.R:\
Contains R script used to produce viral and prokaryotic abundance plots and vOTU count and rank order plots. Requires intermediate files votus_cov75thres.txt, metadata.csv, and qPCR_data.csv.\
This script requires the following packages:\
        ggplot2 3.4.2\
        reshape2 1.4.4\
        readr 2.1.4\
        grid 4.0.2\
        gridExtra 2.3\
        ggpubr 0.6.0

host_plot.R:\
Contains R script used to produce host phyla plot and biogeochemical properties. Requires intermediate files votus_cov75thres.txt, metadata.csv, HostPrediction.tsv, and hostproperties.tsv.\
This script requires the following packages:\
        ggplot2 3.4.2\
        readr 2.1.4\
        gridExtra 2.3\
        ggpubr 0.6.0\
        ggvenn 0.1.10



# Methods

## Mapping raw reads to PIGEON and GOV 2.0 Database

The PIGEON Database was retrieved from:
https://datadryad.org/stash/dataset/doi:10.25338 and indexed using bwa
mem:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaPIGEONindex
    module load bwa/0.7.17-r23kej
    bwa index -p PIGEONv1.0 PIGEONv1.0.fa.gz

For each file in

    /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/rawreads

the following command was used to map the raw read to the PIGEON
database:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaPIGEONindex/sam
    module load bwa/0.7.17-r23kej
    bwa mem -t 4 PIGEONv1.0 \
    /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/rawreads/JK03-#_S#_L002_R1_001.fastq \
    /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/rawreads/JK03-#_S#_L002_R2_001.fastq \
    > JK03-#_S#_L002_aln-pe_001.sam

Next, to prepare the files for samtools coverage, the output sam files
were sorted and output as bam files using the following command for each
raw read sam file:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaPIGEONindex/sam/bam
    samtools sort -@4 \
    /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaPIGEONindex/sam/JK03-#_S#_L002_aln-pe_001.sam \
    -o JK03-#_S#_L002_aln-pe_001.bam

The coverage for each read was then obtained using the following command
on each bam file:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaPIGEONindex/sam/bam/rawreads_coverage
    samtools coverage /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaPIGEONindex/sam/bam/JK03-#_S#_L002_aln-pe_001.bam \
    -o JK03-#_S#_L002_aln-pe_001_coverage

The output was a file containing columns for rname, startpos, endpost,
numreads, covbases, coverage, meandepth, meanbaseq, meanmapq. We were
interested in the coverage value, and used the following commands to
obtain sequences, that had \> 50, \> 60, \> 70, and \> 75 percent
coverage:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaPIGEONindex/sam/bam/rawreads_coverage/coverage_cutoffs

    for FILENAME in rawreads_coverage/*; do awk '{ if($6>50){print $0 "\t" FILENAME}}' $FILENAME \
    >> coverage50.txt; done
    for FILENAME in rawreads_coverage/*; do awk '{ if($6>60){print $0 "\t" FILENAME}}' $FILENAME \
    >> coverage60.txt; done
    for FILENAME in rawreads_coverage/*; do awk '{ if($6>70){print $0 "\t" FILENAME}}' $FILENAME \
    >> coverage70.txt; done
    for FILENAME in rawreads_coverage/*; do awk '{ if($6>75){print $0 "\t" FILENAME}}' $FILENAME \
    >> coverage75.txt; done


First, the GOV 2.0 larger than 5kb or circular database was indexed
using bwa:

        cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaGOVindex
        module load bwa/0.7.17-r23kej
        bwa index -p GOV2_5KB GOV2_viral_populations_larger_than_5KB_or_circular.fasta

### Mapping

Next, the raw reads were each mapped to the database:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaGOVindex
    bwa mem -t 4 \
    GOV2_5KB /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/rawreads/$var1\_L002_R1_001.fastq \
    /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/rawreads/$var1\_L002_R2_001.fastq > \
    $var1\_L002_aln-pe_001.sam

where $var1$ is the prefix of the file.

Reads with \> 75 percent coverage were found using:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaGOVindex/sam

    module load samtools/1.14-3alo66

    # find coverage using samtools
    for FILENAME in *.sam; do
            samtools sort -@4 $FILENAME -o bam/$FILENAME\.bam
            samtools coverage bam/$FILENAME\.bam -o bam/coverage/$FILENAME\_coverage.txt; done

    # get coverage above 75
    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/PIG_GOV_db/bwaGOVindex/sam/bam/coverage
    for FILENAME in *.txt; do
            awk '{ if($6>75){print $0 "\t" FILENAME}}' $FILENAME >> coverage75.txt; done

7 reads mapped to GOV with 75% coverage. Of those 7 reads, 4 of them
were unique sequences. There were some GOV sequences that multiple raw
reads mapped to.

Blastn GUI with megablast option to align two or more sequences.

1 of the 2 PIGEON sequences matched to 1 of the 4 GOV 5KB sequences.

Therefore, we have 5 total unique sequences from raw reads that mapped
to PIGEON and GOV 5KB.

## Gene Prediction using Prodigal

File: viral_virsorter_cdhit.fna

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/geneprediction
    prodigal -i /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/viral_virsorter_cdhit.fna \
    -o viral_virsorter_cdhit_genes -a viral_virsorter_cdhit_proteins



## Gene Annotation with Blast

Blast-Plus version 2.13.0 was used for gene annotation of the output
\*.gff file from prodigal. This was done using UniProt available dataset
swissprot.fa from [swissprot
link](# https://www.uniprot.org/help/download). This was used as the
database for the following Blast commands.

The Uniprot database was obtained using:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/geneprediction/annotation
    wget \
    https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

The Uniprot database was indexed using:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/geneprediction/annotation
    gunzip uniprot_sprot.fasta.gz
    makeblastdb -in uniprot_sprot.fasta \
    -dbtype prot -parse_seqids -out uniprot_sprot -title "uniprot_sprot"

Next, the output from gene prediction was annotated using blast against
the Uniprot database:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/geneprediction/annotation
    blastp -query \
    /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/geneprediction/viral_virsorter_cdhit_proteins  \
    -db uniprot_sprot -num_threads 4 -out viral_virsorter_uniprot_blast



## Increasing Sequences

The VirSorter prediction sequences can be found at:

    /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/viral_virsorter_cdhit.fna

The trimmed reads can be found at:

    /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/trimmed

First, the VirSorter prediction file was indexed using bwa to prepare
the database for mapping:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db
    module load bwa/0.7.17-r23kej
    bwa index -p viral_virsorter /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/viral_virsorter_cdhit.fna

### Mapping

Next, the trimmed reads were each mapped to the database:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db
    module load bwa/0.7.17-r23kej
    bwa mem -t 4 viral_virsorter /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/trimmed/$var1\_L002_R1_001.trimmed.fastq /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/trimmed/$var1\_L002_R2_001.trimmed.fastq > $var1\_L002_aln-pe_001.trimmed.sam

where $var1$ is the prefix of the file.

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db/sam

    module load samtools/1.14-3alo66

    for FILENAME in *.sam; do
        samtools sort -@4 $FILENAME -o bam/$FILENAME\.bam

To filter bam files to only include the reads that mapped, the following
commands were used:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db/sam/bam
    module load samtools/1.14-3alo66
    for FILENAME in *sam.bam; do
        samtools view -b -F 4 $FILENAME > mapped_$FILENAME; done

To convert the mapped reads to fastq files, the following commands were
used:

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db/sam/bam
    module load bedtools2/2.30.0-mlxpvg
    module load samtools/1.14-3alo66

    for FILENAME in mapped*; do
        samtools sort -n -o sort_$FILENAME $FILENAME
        bedtools bamtofastq -i $sort_FILENAME -fq r1_$FILENAME -fq2 r2_$FILENAME; done

### Megahit

Megahit was run on each pair of fastq files to reassemble them using the
following code (JK03-24_S115 as example).

    cd /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db/sam/bam/fastq
    module load anaconda3
    conda activate megahit

    megahit -1 r1_mapped_JK03-24_S115_L002.fq -2 r2_mapped_JK03-24_S115_L002.fq -o megahitoutput24 -t 4

### VirSorter

Viral sequences from each assembly were extracted used VirSorter2 with
the following commands for each file (JK03-24_S115 as example)

    virsorter run -w \
    /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db/virsortoutputs/24.out \
    -i /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db/sam/bam/fastq/megahitoutput24/final.contigs.fa \
    --min-length 1500 -j 4 all

The following command was used to keep the sequences above 5KB
(JK03-24_S115 as example):

    awk -v RS='>[^\n]+\n' 'length() >= 5000 {printf "%s", prt $0} {prt = RT}' \
    /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db/virsortoutputs/24.out/final-viral-combined.fa \
    > \
    /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db/sam/bam/fastq/viral/JK03-24_S115_L002_viral.fa

Resulting viral fa files from the reassemblies can be found at:

    /storage/home/hcoda1/7/iplessis3/p-jweitz3-0/viral_virsorter_db/sam/bam/fastq/viral

### Checkv

Set up checkv

    conda create -n checkv -y
    conda activate checkv
    conda conda install -c conda-forge -c bioconda checkv
    checkv download_database ./
    cd checkv-db-v1.5
    export CHECKVDB=/Users/isabelleduplessis/Downloads/checkv-db-v1.5

Run checkv

    for filename in viral/*.fa; do checkv end_to_end $filename $filename\results -t 2; done

Get non-\"Not determined\" sequences

    for filename in *.fa; do 
        cd /Users/isabelleduplessis/Downloads/checkv-db-v1.5/viral/${filename}\results 
        grep -v "Not-determined" quality_summary.tsv > quality_summary_filtered.tsv;
        grep -A 1 -f <(awk '(NR>1) {print $1}' quality_summary_filtered.tsv) viral/${filename} >> \
    viral/$filename\.toappend.fna
    done
    cat <files> > append.fna
    sed -i '/--/d' append.fna
    cat /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/denovo_mapped_checkv.fna /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/toappend/append.fna > denovo_mapped_checkv_new.fna

### CDHIT

    module load cdhit/4.8.1-n4smhr
    cd-hit-est -i denovo_mapped_checkv_new.fna -o denovo_mapped_checkv_toannotate.fna -c 0.95 -s 0.8

Final output file can be found at:

    /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/saltviralseqs.fna

## Abundance & Diversity Analysis

### BWA index

First, the viral sequences were indexed using BWA:

    module load bwa/0.7.17-r23kej
    cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/abundance
    bwa index -p denovo_mapped_checkv denovo_mapped_checkv.fna

### BWA mem

    module load bwa/0.7.17-r23kej
    cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/abundance
    bwa mem -t 4 denovo_mapped_checkv /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/rawreads/JK03-${samplenumber}_*R1*.fastq /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/rawreads/JK03-${samplenumber}_*R2*.fastq > sample_${samplenumber}.sam

### Samtools Coverage

    module load samtools

    cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/abundance
    samtools view -bS sample_${samplenumber}.sam > sample_${samplenumber}.bam
    samtools sort sample_${samplenumber}.bam > sample_${samplenumber}_sorted.bam
    samtools depth sample_${samplenumber}_sorted.bam > sample_${samplenumber}_depth.txt
    samtools coverage sample_${samplenumber}_sorted.bam > sample_${samplenumber}_coverage.txt

### RPKM and Mean Depth

    cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/abundance/bam
    for i in {1..24}; do totreads=$(samtools view -c sample_${i}_sorted.bam); echo "sample_${i},$totreads"; done > totalreads.csv
    cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/abundance/coverage
    for ARG in {1..24}; do totreads=$(grep "sample_$ARG," /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/abundance/bam/totalreads.csv | cut -f2 -d ","); awk -F'\t' 'NR==1 || $6>=75' "sample_${ARG}_coverage.txt" | grep -v "^#" | sed -e "s/^/sample_$ARG\t/" | awk -v t=$totreads '{print $1"\t"$2"\t"$7"\t"$8"\t"($5*10**9)/($4*t)}'; done > votus_cov75thres.txt;
    sed -i '1 i\#sample\tvotus\tcoverage\tmeandepth\trpkm' votus_cov75thres.txt 
    cd /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/abundance/bam
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

The results of each file were concatenated. The following command output
the number of sequences that belonged to each unique host taxonomy:

    awk -F ',' '{print $3}' host_prediction_to_genus_m90_all.csv | sort | uniq -c | sort -nr > hostfreq_genus.txt

55 unique hosts were identified at the genus level, and 64 out of 769
viral sequences were matched to a host genus. Results can be found at

    /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/iphop_results


## Relevant Gene Annotation - searching for auxiliary metabolic genes

Gene Prediction with Prodigal

    prodigal -i /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/saltviralseqs.fna -o saltviralgenes -a saltviralproteins

Annotation with HMMER

    while read -r line; do sbatch --export=id=$line hmmr2.sbatch; done < /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/relevant_gene_annotation/relevant_metabolic_genes.txt

hmmr2.sbatch:

    hmmsearch --domtblout $id\.domtblout --tblout $id\.tblout /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/kofam_profiles/$id\.hmm /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/iplessis3/saltviralproteins.faa > $id\1.out

Format output

    ls *domtblout | while read file; do /storage/coda1/p-jweitz3/0/shared/PhageGenomicsHI/salt/relevant_gene_annotation/parse_hmmr_domtblout.sh $file; done 

Get significant hits

    cat *domtblout.tab | awk -F'\t' '$7<0.001' |  awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"($17-$16+1)/($6)}' | awk '{a[$1"\t"$2"\t"$3"\t"$4"\t"$5]+=$6}END{for(i in a) print i"\t"a[i]}' | awk -F'\t' '$6>=0.5' > significant_hits

Found 3 significant hits
