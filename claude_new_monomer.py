from typing import List, NamedTuple
from Bio import SeqIO  # Import SeqIO for FASTA parsing
from Bio.Seq import Seq as BioSeq  # Import BioSeq for reverse complement

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
    matches: int = 0
    total_length: int = 0

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
                is_reverse=True,
                matches=matches,    
                total_length=total_length
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
                        cur_aln = cur_aln._replace(matches=matches)
                        cur_aln = cur_aln._replace(total_length=total_length)
                        ans.append(cur_aln)
                        j = monomers_num
                        k = 0
                    else:
                        cur_aln = cur_aln._replace(start_pos=i)
                        cur_aln = cur_aln._replace(match_percentage=(matches / total_length) * 100)
                        cur_aln = cur_aln._replace(matches=matches)
                        cur_aln = cur_aln._replace(total_length=total_length)
                        ans.append(cur_aln)
                        i -= 1
    
    return list(reversed(ans))

def main():
    # Parse the read sequence from a FASTA file
    read_fasta = "test_data/read.fa"  # Path to the read FASTA file
    # monomer_fasta = "test_data/DXZ1_star_monomers.fa"  # Path to the monomer FASTA file
    # output_file = "alignment_results.txt"  # Output file to save alignment results
    monomer_fasta = "test_data/DXZ1_star_monomers_new.fa"  # Path to the monomer FASTA file
    output_file = "alignment_results_new.txt"  # Output file to save alignment results

    # Read the single sequence for the read
    read_record = next(SeqIO.parse(read_fasta, "fasta"))
    read_id = ReadId(read_record.id)
    read = Seq(read_id, str(read_record.seq))

    # Read all monomer sequences
    monomer_records = SeqIO.parse(monomer_fasta, "fasta")
    monomers = []
    for record in monomer_records:
        original_seq = str(record.seq)
        reverse_complement_seq = str(BioSeq(original_seq).reverse_complement())
        
        record_id = record.id.split("_")[0] # A for A_0_DXZ1*_doubled/1978_2147/R

        # Add original sequence
        monomers.append(Seq(ReadId(record_id), original_seq))
        # Add reverse complement with a modified name to distinguish it
        monomers.append(Seq(ReadId(f"{record_id}_rc"), reverse_complement_seq))

    # Perform alignment
    alignments = align_part_classic_dp(read, monomers)

    # Save alignment results to a file
    with open(output_file, "w") as f:
        for aln in alignments:
            result = (f"Monomer: {aln.monomer_name}, "
                      f"Start: {aln.start_pos}, "
                      f"End: {aln.end_pos}, "
                      f"Identity Score: {aln.identity}, "
                      f"Match Percentage: {aln.match_percentage:.2f}%, "
                      f"Matches: {aln.matches}, "
                      f"Total Length: {aln.total_length}\n")
            f.write(result)
            print(result, end="")  # Print to console as well

    print(f"\nAlignment results have been saved to {output_file}")

if __name__ == "__main__":
    main()
