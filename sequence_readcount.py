from collections import Counter

def count_reads_in_fastq(fastq_file, output_file):
    # Open the FASTQ file and read sequences
    with open(fastq_file, 'r') as f:
        sequences = []
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence lines are every 4th line starting from 1
                sequences.append(line.strip())

    # Count the occurrences of each sequence
    sequence_counts = Counter(sequences)

    # Sort sequences by their counts in descending order
    sorted_sequence_counts = sorted(sequence_counts.items(), key=lambda x: x[1], reverse=True)

    # Write the sequences and their read counts to the output file
    with open(output_file, 'w') as out_f:
        for sequence, count in sorted_sequence_counts:
            out_f.write(f"{sequence}\t{count}\n")

# Example usage
fastq_file = 'TTAC.extendedFrags.fastq'  # Replace with your merged FASTQ file path
output_file = 'TTAC_count.txt'  # Replace with desired output file path
count_reads_in_fastq(fastq_file, output_file)
