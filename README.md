## Step 1
gen_reads.py: Generates the NGS reads with uniform error rate 0.01. It takes fasta file as input. Use following command to generate reads in FASTQ format.
python3 gen_reads.py ref.fa

## Step 2
Map generated reads to the Genome file. Use following commands:
bwa aln ref.fa reads.fastq 1> reads.sai
bwa samse -f reads.sam ref.fa reads.sai reads.fastq

## Step 3
get_error_rate.py: Generates error rate from mapping file (reads.sam) usf following command:
python3 get_error_rate.py reads.sam
