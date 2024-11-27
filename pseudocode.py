def string_decomposition(sequence, blocks, match_score=1, mismatch_penalty=-1, indel_penalty=-2):
    """
    Perform string decomposition of `sequence` using `blocks` (monomers).
    Args:
        sequence (str): The sequence to decompose.
        blocks (list of str): List of known monomers.
        match_score (int): Score for a match.
        mismatch_penalty (int): Penalty for a mismatch.
        indel_penalty (int): Penalty for insertions or deletions.

    Returns:
        tuple: Total score of the best decomposition and the sequence of monomers.
    """
    # Initialize DP table and path tracker
    sequence_len = len(sequence)
    dp = {}
    path = {}

    for block_id, block in enumerate(blocks):
        block_len = len(block)
        dp[block_id] = [[-float('inf')] * (sequence_len + 1) for _ in range(block_len + 1)]
        path[block_id] = [[None] * (sequence_len + 1) for _ in range(block_len + 1)]

        # Base case initialization
        for i in range(block_len + 1):
            dp[block_id][i][0] = -i * indel_penalty
        for j in range(sequence_len + 1):
            dp[block_id][0][j] = 0

    # Fill the DP table for all blocks
    for block_id, block in enumerate(blocks):
        block_len = len(block)
        for i in range(1, block_len + 1):
            for j in range(1, sequence_len + 1):
                match = dp[block_id][i - 1][j - 1] + (match_score if block[i - 1] == sequence[j - 1] else mismatch_penalty)
                delete = dp[block_id][i - 1][j] + indel_penalty
                insert = dp[block_id][i][j - 1] + indel_penalty

                dp[block_id][i][j] = max(match, delete, insert)

                # Track the path
                if dp[block_id][i][j] == match:
                    path[block_id][i][j] = (block_id, i - 1, j - 1)
                elif dp[block_id][i][j] == delete:
                    path[block_id][i][j] = (block_id, i - 1, j)
                elif dp[block_id][i][j] == insert:
                    path[block_id][i][j] = (block_id, i, j - 1)

    # Add transitions between blocks
    for j in range(sequence_len + 1):
        for curr_block_id, curr_block in enumerate(blocks):
            for next_block_id, next_block in enumerate(blocks):
                if dp[curr_block_id][len(curr_block)][j] != -float('inf'):
                    next_score = dp[curr_block_id][len(curr_block)][j]
                    if dp[next_block_id][0][j] < next_score:
                        dp[next_block_id][0][j] = next_score
                        path[next_block_id][0][j] = (curr_block_id, len(curr_block), j)

    # Find the best score and reconstruct the path
    best_score = -float('inf')
    best_position = None
    for block_id in range(len(blocks)):
        if dp[block_id][len(blocks[block_id])][sequence_len] > best_score:
            best_score = dp[block_id][len(blocks[block_id])][sequence_len]
            best_position = (block_id, len(blocks[block_id]), sequence_len)

    # Backtrack to reconstruct the monomer sequence
    result_sequence = []
    current = best_position
    while current:
        block_id, i, j = current
        if i == 0:  # If at the start of a block, switch
            if path[block_id][i][j]:
                prev_block_id, _, _ = path[block_id][i][j]
                if prev_block_id != block_id:
                    result_sequence.append(blocks[prev_block_id])
        current = path[block_id][i][j]

    result_sequence.reverse()

    return best_score, result_sequence


# Example usage
sequence = "ATGCTAGC"
blocks = ["ATGC", "TAGC", "GCTA", "CAGT"]
score, monomers = string_decomposition(sequence, blocks)
print(f"Score: {score}")
print(f"Monomers: {monomers}")
