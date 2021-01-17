import numpy as np
import  pandas as pd

def hamming_distance(sequence_1, sequence_2):
    if len(sequence_1) != len(sequence_2):
        print("Hamming distance is only defined for sequences of equal length.")
        return "not defined"

    mismatch_counter = 0
    for i in range(len(sequence_1)):
        if sequence_1[i] != sequence_2[i]:
            mismatch_counter += 1
    hamming_distance = mismatch_counter
    return hamming_distance

def percent_identity(align_1, align_2):
    """Calculates percent identity by dividing the numbers of matching positions with the length of the alignment."""
    matches = 0
    for i in range(len(align_1)):
        if align_1[i] == align_2[i]:
            matches+= 1
    percent_identity = matches / len(align_1)
    return percent_identity
    

def local_alignment(sequence_1, sequence_2):
    # define scores
    match_score =  2
    mismatch_score = -1
    gap_penalty = 2

    # define actions
    stop = 0
    up = 1
    left = 2
    diag = 3

    len_seq_1 = len(sequence_1)
    len_seq_2 = len(sequence_2)
    
    # initiate matrices
    F = np.zeros(shape=(len_seq_1 +1, len_seq_2 +1))
    trace = np.zeros(shape=(len_seq_1  +1 ,  len_seq_2+1))


    # fill the matrices
    for i in range(1,len_seq_1+1):
        for j in range(1, len_seq_2+1):
            if sequence_1[i-1] == sequence_2[j-1]:
                score = F[i-1][j-1] + match_score
            else:
                score = F[i-1][j-1] + mismatch_score
            trace[i][j] = diag
            
            go_up_score = F[i-1][j] - gap_penalty
            if go_up_score > score:
                score = go_up_score
                trace[i][j] = up

            go_left_score = F[i][j-1] - gap_penalty
            if go_left_score > score:
                score = go_left_score
                trace[i][j] = left

            if score > 0:
                F[i][j] = score
            else:
                F[i][j] = 0
                trace[i][j] = 0
                
    #print F matrix
    cols = [' '] + [i  for i in sequence_2]
    rows = [' '] + [i  for i in sequence_1]

    F_df = pd.DataFrame(F, index = rows, columns=cols, dtype=int)
    print("F matrix: local alignment")
    print(F_df)

    #print trace matrix
    trace_df= pd.DataFrame(trace, index = rows, columns=cols, dtype=int)
    print("trace matrix: local alignment")
    print(trace_df)

    # find starting point for trace back
    max = 0
    for i in range(F.shape[0]):
        for j in range(F.shape[1]):
            if F[i][j] >= max:
                max = F[i][j]
                max_i = i
                max_j = j

    #trace back from the position that have the highest F value
    i = max_i 
    j = max_j 

    align_1 = []
    align_2 = []

    #  take action and increment position (i,j)
    while trace[i][j] != stop or F[i][j] != 0:
        if trace[i][j] == diag:
            align_1.append(sequence_1[i-1])
            align_2.append(sequence_2[j-1])
            i -= 1
            j -= 1        
        elif trace[i][j] == left:
            align_1.append('-')
            align_2.append(sequence_2[j-1])
            j -= 1      
        elif trace[i][j] == up:
            align_1.append(sequence_1[i-1])
            align_2.append('-')
            i -= 1
    
    # print alignments
    match = ["|"  if align_1[i]==align_2[i] else " " for i in range(len(align_1))]
    print("Alignment:")
    print('  '.join(align_1[::-1]))
    print('  '.join(match[::-1]))
    print('  '.join(align_2[::-1]))

    # print distance and percent  identity
    percent_id = percent_identity(align_1, align_2)
    print(f"Percent identiy is {percent_id}")

    hamming_dist = hamming_distance(sequence_1, sequence_2)
    print(f"Hamming distance is: {hamming_dist}")

#run 
local_alignment("PAWHEAE", "HDAGAWGHEQ")
