from typing import List, NamedTuple

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
    is_reverse: bool = True

def align_part_classic_dp(read: Seq, monomers: List[Seq], 
                           ins: int = -2, 
                           delete: int = -1, 
                           match: int = 2, 
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
    
    while i >= 0:
        if k == len(dp[i][j]) - 1 and j != monomers_num and monomer_changed:
            cur_aln = MonomerAlignment(
                monomer_name=monomers[j].read_id.name, 
                read_name=read.read_id.name, 
                end_pos=i, 
                start_pos=i, 
                identity=dp[i][j][k], 
                is_reverse=True
            )
            monomer_changed = False
        
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
            elif i != 0 and dp[i][j][k] == dp[i-1][j][k] + ins:
                i -= 1
            else:
                mm_score = match if monomers[j].seq[k] == read.seq[i] else mismatch
                if i != 0 and k != 0 and dp[i][j][k] == dp[i-1][j][k-1] + mm_score:
                    i -= 1
                    k -= 1
                else:
                    monomer_changed = True
                    if (i != 0 and 
                        dp[i][monomers_num][0] + k * delete + mm_score == dp[i][j][k]):
                        cur_aln = cur_aln._replace(start_pos=i)
                        cur_aln = cur_aln._replace(identity=cur_aln.identity - dp[i][monomers_num][0])
                        ans.append(cur_aln)
                        j = monomers_num
                        k = 0
                    else:
                        cur_aln = cur_aln._replace(start_pos=i)
                        ans.append(cur_aln)
                        i -= 1
    
    return list(reversed(ans))

# Example usage
def main():
    # Example setup (you would replace this with your actual data)
    read_id = ReadId("read1")
    monomer_ids = [ReadId(f"monomer{i}") for i in range(3)]
    
    read = Seq(read_id, "ATCGGATCCCGT")
    monomers = [
        Seq(monomer_ids[0], "ATCG"),
        Seq(monomer_ids[1], "GATC"),
        Seq(monomer_ids[2], "CCGT")
    ]
    
    alignments = align_part_classic_dp(read, monomers)
    
    for aln in alignments:
        print(f"Monomer: {aln.monomer_name}, "
              f"Read: {aln.read_name}, "
              f"Start: {aln.start_pos}, "
              f"End: {aln.end_pos}, "
              f"Identity: {aln.identity}")

if __name__ == "__main__":
    main()