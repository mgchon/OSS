from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import collections

def process_sequences(fastq_path, primers):
    processed_records = []
    for record in SeqIO.parse(fastq_path, "fastq"):
        seq = str(record.seq)
        for primer_pair in primers:
            start_seq, end_seq = primer_pair
            start_idx = seq.find(start_seq)
            end_idx = seq.find(end_seq, start_idx + 1)
            if start_idx != -1 and end_idx != -1:
                end_idx += len(end_seq)
                processed_seq = seq[start_idx:end_idx]
                processed_qual = record.letter_annotations["phred_quality"][start_idx:end_idx]
                new_record = SeqRecord(Seq(processed_seq), id=record.id, description="", letter_annotations={"phred_quality": processed_qual})
                processed_records.append(new_record)
                break
    return processed_records

def write_to_fastq(sequences, output_path):
    with open(output_path, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fastq")

def count_sequence_lengths(sequences):
    length_counts = collections.Counter(len(seq.seq) for seq in sequences)
    return length_counts

# Define the primers
primers = [("AGTGCAACAAGTCAATCCGT", "AATTGAATGCTTGCTTGCCG"),
           ("CGGCAAGCAAGCATTCAAT", "ACGGATTGACTTGTTGCACT")]

# Process the FASTQ file
fastq_file = "id20.fastq" # replace with your FASTQ file path
processed_records = process_sequences(fastq_file, primers)

# Write processed sequences to a new FASTQ file
output_fastq = "id20_after_trimming.fastq" # replace with your desired output file path
write_to_fastq(processed_records, output_fastq)
print(f"Processed sequences have been saved to '{output_fastq}'.")

# Filter out uncut sequences
filtered_records = [record for record in processed_records if record.seq]

# Sequence Length Analysis
seq_length_counts = count_sequence_lengths(filtered_records)

# Write lengths to a file
with open("id20_seq_cnt.txt", "w") as f:
    for length, count in seq_length_counts.items():
        f.write(f"{length} {count}\n")

print("Sequence length analysis is complete. Results are saved in 'sequence_lengths.txt'")
