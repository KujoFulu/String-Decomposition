from typing import List, NamedTuple
from Bio import SeqIO  # Import SeqIO for FASTA parsing
from Bio.Seq import Seq as BioSeq  # Import BioSeq for reverse complement
import pandas as pd
import multiprocessing
import os
import itertools

class ReadId:
    def __init__(self, name: str):
        self.name = name

class Seq:
    def __init__(self, read_id: ReadId, sequence: str):
        self.read_id = read_id
        self.seq = list(sequence)

class MonomerAlignment(NamedTuple):
    monomer_name: str
    read_name: str
    end_pos: int
    start_pos: int
    identity: float
    match_percentage: float
    is_reverse: bool = True

def align_part_classic_dp(read: Seq, monomers: List[Seq], 
                           ins: int = -1, 
                           delete: int = -1, 
                           match: int = 1, 
                           mismatch: int = -1) -> List[MonomerAlignment]:
    """
    Perform sequence alignment using dynamic programming
    
    Args:
        read: The input sequence to align
        monomers: List of monomer sequences to align against
        ins: Insertion penalty
        delete: Deletion penalty
        match: Matching score
        mismatch: Mismatching score
    
    Returns:
        List of MonomerAlignment objects
    """
    INF = -1000000
    monomers_num = len(monomers)
    
    # Initialize 3D DP table
    dp = [[[INF for _ in range(len(m.seq))] for m in monomers] + [[INF]] 
          for _ in range(len(read.seq))]
    
    # Initialize first row (base case for first character)
    for j, m in enumerate(monomers):
        dp[0][j][0] = match if m.seq[0] == read.seq[0] else mismatch
        
        # Handle initial deletions
        for k in range(1, len(m.seq)):
            mm_score = match if monomers[j].seq[k] == read.seq[0] else mismatch
            dp[0][j][k] = max(dp[0][j][k-1] + delete, delete * (k-1) + mm_score)
    
    # Dynamic Programming core
    for i in range(1, len(read.seq)):
        # Update termination column
        dp[i][monomers_num][0] = max(
            (dp[i][monomers_num][0], 
             max(dp[i-1][j][len(monomers[j].seq) - 1] for j in range(monomers_num)))
        )
        
        # Main DP loop
        for j in range(monomers_num):
            for k in range(len(monomers[j].seq)):
                score = INF
                mm_score = match if monomers[j].seq[k] == read.seq[i] else mismatch
                
                # Various transition possibilities
                if dp[i][monomers_num][0] > INF:
                    score = max(score, dp[i][monomers_num][0] + mm_score + k * delete)
                
                if k > 0:
                    if dp[i-1][j][k-1] > INF:
                        score = max(score, dp[i-1][j][k-1] + mm_score)
                    if dp[i-1][j][k] > INF:
                        score = max(score, dp[i-1][j][k] + ins)
                    if dp[i][j][k-1] > INF:
                        score = max(score, dp[i][j][k-1] + delete)
                
                dp[i][j][k] = score
    
    # Find best alignment
    max_score = INF
    best_m = monomers_num
    for j in range(monomers_num):
        current_score = dp[len(read.seq)-1][j][len(monomers[j].seq) - 1]
        if max_score < current_score:
            max_score = current_score
            best_m = j
    
    # Traceback to construct alignment
    ans = []
    i = len(read.seq) - 1
    j = best_m
    k = len(dp[i][j]) - 1
    monomer_changed = True
    cur_aln = None
    matches = 0
    total_length = 0
    
    while i >= 0:
        if k == len(dp[i][j]) - 1 and j != monomers_num and monomer_changed:
            cur_aln = MonomerAlignment(
                monomer_name=monomers[j].read_id.name, 
                read_name=read.read_id.name, 
                end_pos=i, 
                start_pos=i, 
                identity=dp[i][j][k], 
                match_percentage=0.0,  # Placeholder, will be updated later
                is_reverse=True
            )
            monomer_changed = False
            matches = 0
            total_length = 0
        
        if j == monomers_num:
            if i != 0:
                for p in range(len(dp[i-1])):
                    if dp[i-1][p][len(dp[i-1][p]) - 1] == dp[i][j][k]:
                        i -= 1
                        j = p
                        k = len(dp[i][j]) - 1
                        break
            else:
                i -= 1
        else:
            if k != 0 and dp[i][j][k] == dp[i][j][k-1] + delete:
                k -= 1
                total_length += 1
            elif i != 0 and dp[i][j][k] == dp[i-1][j][k] + ins:
                i -= 1
                total_length += 1
            else:
                mm_score = match if monomers[j].seq[k] == read.seq[i] else mismatch
                if mm_score == match:
                    matches += 1
                total_length += 1
                
                if i != 0 and k != 0 and dp[i][j][k] == dp[i-1][j][k-1] + mm_score:
                    i -= 1
                    k -= 1
                else:
                    monomer_changed = True
                    if (i != 0 and 
                        dp[i][monomers_num][0] + k * delete + mm_score == dp[i][j][k]):
                        cur_aln = cur_aln._replace(start_pos=i)
                        cur_aln = cur_aln._replace(identity=cur_aln.identity - dp[i][monomers_num][0])
                        cur_aln = cur_aln._replace(match_percentage=(matches / total_length) * 100)
                        ans.append(cur_aln)
                        j = monomers_num
                        k = 0
                    else:
                        cur_aln = cur_aln._replace(start_pos=i)
                        cur_aln = cur_aln._replace(match_percentage=(matches / total_length) * 100)
                        ans.append(cur_aln)
                        i -= 1
    
    return list(reversed(ans))

def process_single_read(read_record: SeqIO.SeqRecord, monomers: List[Seq], output_dir: str) -> None:
    """
    Process a single read and save alignment results to a TSV file.
    
    Args:
        read_record (SeqIO.SeqRecord): The read sequence record
        monomers (List[Seq]): List of monomer sequences
        output_dir (str): Directory to save output TSV files
    """
    # Create read sequence object
    read_id = ReadId(read_record.id)
    read = Seq(read_id, str(read_record.seq))

    # Perform alignment
    alignments = align_part_classic_dp(read, monomers)

    # Convert alignments to a pandas DataFrame
    df = pd.DataFrame([
        {
            'Monomer': aln.monomer_name,
            'Read_Name': aln.read_name,
            'Start_Position': aln.start_pos,
            'End_Position': aln.end_pos,
            'Identity': aln.identity,
            'Match_Percentage': f"{aln.match_percentage:.2f}",
            'Is_Reverse': aln.is_reverse
        } for aln in alignments
    ])

    # Create output filename (replace any characters not suitable for filenames)
    safe_filename = "".join(
        [c if c.isalnum() or c in ('-', '_') else '_' for c in read_record.id]
    )
    output_tsv = os.path.join(output_dir, f"{safe_filename}_alignments.tsv")

    # Write to TSV
    df.to_csv(output_tsv, sep='\t', index=False)
    
    print(f"Processed read {read_record.id}, saved to {output_tsv}")

def load_monomers(monomer_fasta: str) -> List[Seq]:
    """
    Load monomer sequences from a FASTA file, including both original and reverse complement.
    
    Args:
        monomer_fasta (str): Path to the monomer FASTA file
    
    Returns:
        List[Seq]: List of monomer sequences
    """
    monomer_records = SeqIO.parse(monomer_fasta, "fasta")
    monomers = []
    for record in monomer_records:
        original_seq = str(record.seq)
        reverse_complement_seq = str(BioSeq(original_seq).reverse_complement())
        
        record_id = record.id.split("_")[0]

        # Add original sequence
        monomers.append(Seq(ReadId(record_id), original_seq))
        # Add reverse complement with a modified name to distinguish it
        monomers.append(Seq(ReadId(f"{record_id}_rc"), reverse_complement_seq))
    
    return monomers

def main():
    # Configuration
    read_fasta = "../new_data/cenX_aligned_reads/cenX_aligned_reads.fasta"  # Path to the read FASTA file with multiple reads
    monomer_fasta = "test_data/DXZ1_star_monomers_new.fa"  # Path to the monomer FASTA file
    output_dir = "100_new_monomers"  # Directory to store output TSV files
    num_cores = 5  # Number of CPU cores to use
    max_reads = 100  # Maximum number of reads to process

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load monomer sequences
    monomers = load_monomers(monomer_fasta)

    # Read first 500 sequences from the read FASTA file
    read_records = list(itertools.islice(SeqIO.parse(read_fasta, "fasta"), max_reads))

    # Use multiprocessing to process reads in parallel
    with multiprocessing.Pool(processes=num_cores) as pool:
        # Prepare arguments for each read
        args = [(record, monomers, output_dir) for record in read_records]
        
        # Map the processing function to all reads
        pool.starmap(process_single_read, args)

    print(f"Processed {len(read_records)} reads across {num_cores} cores.")
    print(f"Alignment results saved in {output_dir}")

if __name__ == "__main__":
    main()
