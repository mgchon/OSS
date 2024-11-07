from Bio import SeqIO

def compare_sequences_to_all_references(sequence_file, reference_file, output_file):
    # Load all reference sequences from the FASTA file
    reference_seqs = []
    for record in SeqIO.parse(reference_file, "fasta"):
        reference_seqs.append(str(record.seq))
    
    # Load sequences and their counts from the sequence file
    sequences_with_counts = []
    with open(sequence_file, 'r') as f:
        for line in f:
            sequence, count = line.strip().split('\t')
            sequences_with_counts.append((sequence, int(count)))
    
    # Function to count mismatches between two sequences
    def count_mismatches(seq1, seq2):
        mismatches = []
        min_len = min(len(seq1), len(seq2))
        for i in range(min_len):
            if seq1[i] != seq2[i]:
                mismatches.append((i, seq1[i], seq2[i]))  # Record mismatch (position, base in seq1, base in seq2)
        # Handle cases where the sequences are of different lengths (insertions/deletions)
        if len(seq1) > len(seq2):
            for i in range(min_len, len(seq1)):
                mismatches.append((i, seq1[i], '-'))
        elif len(seq2) > len(seq1):
            for i in range(min_len, len(seq2)):
                mismatches.append((i, '-', seq2[i]))
        return mismatches

    # Compare each sequence to all reference sequences and find the best match
    best_matches = []
    
    for sequence, count in sequences_with_counts:
        best_match = None
        min_mismatches = float('inf')
        best_mismatch_details = []
        
        # Compare the sequence to each reference sequence
        for ref_seq in reference_seqs:
            mismatches = count_mismatches(sequence, ref_seq)
            mismatch_count = len(mismatches)
            
            # If this reference sequence has fewer mismatches, use it as the best match
            if mismatch_count < min_mismatches:
                min_mismatches = mismatch_count
                best_match = ref_seq
                best_mismatch_details = mismatches
        
        # Add the best matching result for the current sequence
        best_matches.append((sequence, count, min_mismatches, best_mismatch_details, best_match))
    
    # Sort by read count in descending order
    best_matches.sort(key=lambda x: -x[1])

    # Write the results to the output file
    with open(output_file, 'w') as out_f:
        out_f.write(f"Sequence\tRead Count\tMismatch Count\tBest Match Reference\tDifferences (Position, Sequence Base, Reference Base)\n")
        for sequence, count, mismatch_count, mismatches, best_match in best_matches:
            differences_str = "; ".join([f"{pos}:{seq_base}->{ref_base}" for pos, seq_base, ref_base in mismatches])
            out_f.write(f"{sequence}\t{count}\t{mismatch_count}\t{best_match}\t{differences_str}\n")

# Example usage
sequence_file = 'TTAC_count.txt'  # File containing sequences and counts
reference_file = 'new_reference.fasta'  # Reference FASTA file
output_file = 'TTAC_hamming_total.txt'  # Output file for best matches
compare_sequences_to_all_references(sequence_file, reference_file, output_file)
