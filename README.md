## Oligonucleotide subsets selection by single resolution barcode identification

### NGS data processing

#### Merge
flash R1.fastq R2.fastq -o output

#### Creating a text file listing sequences and their counts based on a merged FASTQ file
awk 'NR%4==2{print $0}' input.fastq | sort | uniq -c | awk '{print $2 "\t" $1}' > output.txt

#### Aligning sequences to the reference (Aloowing mismatches of three)
/home/align/bowtie-1.3.1-linux-x86_64/bowtie-build reference.fasta reference_index
align/bowtie-1.3.1-linux-x86_64/bowtie -v 3 -x reference_index -q R1.fastq -S output.sam
samtools view -S -b output.sam > output.bam;
samtools sort output.bam -o output.bam
samtools index output.bam
samtools idxstats output_1_sorted.bam | cut -f 1,3 > output_1_result.txt

#### Measuring the Hamming distance with reference
python error_check.py

#### Downsampling the FASTQ file
seqtk sample input.fastq 0.1 > output.fastq
