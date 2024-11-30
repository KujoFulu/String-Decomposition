from typing import Dict

#This function is not tested
#This function should returns the dp matrix, the best final score, and the final block
def StringDecomposition(sequence: str, blocks: Dict[str, str], match: int, mismatch: int, indel: int):
    #sequence is the raw centromere seqeunce, blocks is the disctionary of monomers

    #Initialize the dp matrix
    score = {}
    for block_id, block_sequence in blocks:
        score[block_id] = [[-float('inf')] * (len(sequence) + 1) for _ in range(len(block_sequence))]
    #Initialize the common first row of the dp matrix
    score_firstrow = [0] * (len(sequence) + 1)

    #Initialize base cases
    for b in score:
        for i in range(len(score[b])):
            score[b][i] = -i*indel
    
    #for j in range(len(score_firstrow)):
        #score_firstrow[j] = 0
    
    #Recurrence
    #first row in the score matrix
    for b in score:
        for j in range(1,len(sequence)+1):
            match_score = mismatch
            if blocks[b][0] == sequence[j-1]:
                match_score = match
            score[b][0][j] = max(score_firstrow[j]-indel, score[b][0][j-1]-indel, score_firstrow[j-1]+match_score)
    
    #Other parts of the score matrix
    for b in score:
        for i in range(1,len(score[b])):
            for j in range(1,len(sequence)+1):
                match_score = mismatch
                if blocks[b][i] == sequence[j-1]:
                    match_score = match
                score[b][i][j] = max(score[i-1][j]-indel, score[b][i][j-1]-indel, score[i-1][j-1]+match_score)

    #Check the sinks, find the max score
    max_score = -float('inf')
    final_b = ""
    for b in score:
        current_score = score[b][len(score[b])-1][len(sequence)]
        if current_score > max_score:
            max_score = current_score
            final_b = b
    
    #Need further implementation to traceback
    return score, max_score, final_b