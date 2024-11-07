def trim_fastq_sequences(input_file, output_file, target_length):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        while True:
            header = infile.readline().strip()
            if not header:
                break  # End of file
            sequence = infile.readline().strip()
            plus = infile.readline().strip()
            quality = infile.readline().strip()

            # Trim the sequence and corresponding quality score
            trimmed_sequence = sequence[:target_length]
            trimmed_quality = quality[:target_length]

            # Write to output file
            outfile.write(f"{header}\n{trimmed_sequence}\n{plus}\n{trimmed_quality}\n")

# Usage
input_fastq_file = 'id20.fastq'  # Replace with the actual file path
output_fastq_file = 'id20_150.fastq'
desired_length = 150

# Run the function
trim_fastq_sequences(input_fastq_file, output_fastq_file, desired_length)

