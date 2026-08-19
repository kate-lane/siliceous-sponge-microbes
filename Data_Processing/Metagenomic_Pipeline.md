Summary of metagenomic processing pipeline. All commands were run on Woods Hole Oceanographic Institution Poseidon SLURM cluster.

## Trim/QC raw reads
fastqc *gz

bbduk.sh -Xmx1g in1=/raw.d/${i}.1.fastq.gz in2=/raw.d/${i}.2.fastq.gz out1=/raw.d/${i}_bbduk.1.fastq.gz out2=/raw.d/${i}_bbduk.2.fastq.gz minlen=50 qtrim=rl trimq=10 ktrim=r k=25 mink=8 ref=/home/katherine.lane/references/adapters.fa hdist=1 stats=/raw.d/${i}_bbduk.log

fastqc *bbduk*gz

## Assembly
megahit -1 /raw.d/${i}_bbduk.1.fastq.gz -2 /raw.d/${i}_bbduk.2.fastq.gz -o /assembly/${i}_megahit

#retain <1000bp contigs only, bbtools reformat.sh
reformat.sh in=/assembly/${i}_megahit/final.contigs.fa out=/assembly/${i}_min1000.fa minlength=1000


## Binning MAGs
mkdir /assembly/${i}_megahit_min1000
cp /assembly/${i}_min1000.fa /assembly/${i}_megahit_min1000/${i}_min1000.fa

bwa index /assembly/${i}_megahit_min1000/${i}_min1000.fa

bwa mem /assembly/${i}_megahit_min1000/${i}_min1000.fa /raw.d/${i}_bbduk.1.fastq.gz /raw.d/${i}_bbduk.2.fastq.gz > /assembly/${i}_megahit_min1000/aln-pe.sam

samtools view -bu /assembly/${i}_megahit_min1000/aln-pe.sam | samtools sort -@48 -l 6 -O bam -o /assembly/${i}_megahit_min1000/aln-pe.sorted.bam -

jgi_summarize_bam_contig_depths --outputDepth /assembly/${i}_megahit_min1000/depth.txt /assembly/${i}_megahit_min1000/aln-pe.sorted.bam

metabat2 -i /assembly/${i}_megahit_min1000/${i}_min1000.fa -a /assembly/${i}_megahit_min1000/depth.txt -o /assembly/${i}_megahit_min1000/bins_dir/bin

## MAGs QC
checkm2 predict --threads 30 --input /assembly/${i}_megahit_min1000/bins_dir/  --output-directory /assembly/${i}_megahit_min1000/bins_dir/checkm2_results -x fa 

#Retain only bins with Completeness > 70 and Contamination < 10

## MAG abundances
coverm genome -b /assembly/${i}_megahit_min1000/aln-pe.sorted.bam --genome-fasta-directory /assembly/${i}_megahit_min1000/bins_dir/filtered_bins/renamed_bins -x fa --methods mean relative_abundance

## Taxonomic Assignment GTDBTK
gtdbtk classify_wf --genome_dir /proj/hansel-lab/klane/baltic_2024_10_05/assembly/${i}_megahit_min1000/bins_dir/ --out_dir /proj/hansel-lab/klane/baltic_2024_10_05/gtdbtk/${i}_gtdbtk -x fa

## METABOLIC-G
perl METABOLIC-G.pl -in-gn /METABOLIC/${i}_bins -o /METABOLIC/${i}_metabolic-G

## antiSMASH
antismash "$fasta" \
  --output-dir /antismash \
  --genefinding-tool prodigal \
  --cb-general --cb-subclusters --cb-knownclusters \
  --minlength 5000 \
  --asf --pfam2go \
  --clusterhmmer \
  --cc-mibig \
  --rre \
  --tigrfam \
  --cpus 12
