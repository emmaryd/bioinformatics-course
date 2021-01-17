import numpy as np
import pandas as pd

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

def find_paths(trace, 
               sequence_1, 
               sequence_2, 
               position,
               align_1 = [],
               align_2 = [],
               all_alignments=[], 
               stop_indices=[]):
    
    # define actions
    stop = 0
    up = 1
    left = 2
    diag = 3
    
    i, j = position
    action = trace[i][j]
    while action != stop:
        
        # check if there is more than one possible action
        if type(action) == list:
            alt_paths = action
            trace_cpy = trace.copy() # create copy of trace matrix

            # iterate over the additional alternatives (all except the first)
            for alt in alt_paths[1:]: 
                trace_cpy[i][j] = alt # replace mult option with one alternative
                
                # for each alternative action run find_path 
                all_alignments, stop_indices = find_paths(trace_cpy, 
                                                        sequence_1, 
                                                        sequence_2,
                                                        position = (i,j),
                                                        align_1= align_1.copy(),
                                                        align_2= align_2.copy(),
                                                        all_alignments=all_alignments, 
                                                        stop_indices=stop_indices)
            
            # change action to alternatice one
            action = alt_paths[0]

        # Now only one action is possible, take action and increment position (i,j)
        if action == diag:
            align_1.append(sequence_1[i-1])
            align_2.append(sequence_2[j-1])
            i -= 1
            j -= 1        
        elif action == left:
            align_1.append('-')
            align_2.append(sequence_2[j-1])
            j -= 1      
        elif action == up:
            align_1.append(sequence_1[i-1])
            align_2.append('-')
            i -= 1
        action = trace[i][j] # update action
    
    all_alignments.append((align_1, align_2))
    stop_indices.append((i,j))

    return all_alignments, stop_indices   

def global_alignment(sequence_1, sequence_2):
    # define scores
    match_score =  2
    mismatch_score = -1
    gap_penalty = 2

    #define actions,  same as in find_path
    stop = 0
    up = 1
    left = 2
    diag = 3

    len_seq_1 = len(sequence_1)
    len_seq_2 = len(sequence_2)
    #initiate matrices
    F = np.zeros(shape=(len_seq_1 +1, len_seq_2 +1))
    trace = np.zeros(shape=(len_seq_1  +1 ,  len_seq_2+1), dtype=object)
    
    for i in range(1, len_seq_1+1):
        F[i][0] = F[i-1][0] - gap_penalty

    for j in range(1, len_seq_2+1):
        F[0][j] = F[0][j-1] - gap_penalty

    # fill the matrices
    for i in range(1,len_seq_1+1):
        for j in range(1, len_seq_2+1):
            alternatives = []
            if sequence_1[i-1] == sequence_2[j-1]:
                diag_score = F[i-1][j-1] + match_score
            else:
                diag_score = F[i-1][j-1] + mismatch_score
            trace[i][j] = diag

            go_up_score = F[i-1][j] - gap_penalty        
            if go_up_score >= diag_score:
                alternatives.append(up)
                trace[i][j] = up

            go_left_score = F[i][j-1] - gap_penalty
            if go_left_score >= diag_score and go_left_score >= go_up_score :
                alternatives.append(left)
                trace[i][j] = left  
            F[i][j] = max(diag_score, go_up_score, go_left_score)
            
            #Check for alternative paths
            if go_up_score == diag_score:
                alternatives.append(diag)
                trace[i][j] = alternatives
            if go_left_score == diag_score:
                alternatives.append(diag)
                trace[i][j] = alternatives
    #print F matrix
    cols = [' '] + [i  for i in sequence_2]
    rows = [' '] + [i  for i in sequence_1]

    F_df = pd.DataFrame(F, index = rows, columns=cols, dtype=int)
    print("F matrix: global alignment")
    print(F_df)

    #print trace matrix
    trace_df= pd.DataFrame(trace, index = rows, columns=cols)
    print("trace matrix: global alignment")
    print(trace_df)

    #trace back from lower right corner
    i = len_seq_1 
    j = len_seq_2 
    all_alignments, stop_indices = find_paths(trace, 
                                            sequence_1, 
                                            sequence_2,
                                            position = (i,j),
                                            align_1=[],
                                            align_2=[],
                                            all_alignments=[])

    alignment_counter = 0
    for alignment, index in zip(all_alignments, stop_indices):
        alignment_counter += 1
        align_1,  align_2 = alignment
        i,j = index

        #unaligned begining
        while i > 0:
            align_1.append(sequence_1[i-1])
            align_2.append('-')
            i -= 1
        while j > 0:
            align_1.append('-')
            align_2.append(sequence_2[j-1])
            j -= 1

        # print alignments
        match = ["|"  if align_1[i]==align_2[i] else " " for i in range(len(align_1))]
        print(f"Alignment: {alignment_counter}")
        print('  '.join(align_1[::-1]))
        print('  '.join(match[::-1]))
        print('  '.join(align_2[::-1]))
        
        # print distance and percent  identity
        hamming_dist = hamming_distance(sequence_1, sequence_2)
        print(f"Hamming distance is: {hamming_dist}")

        percent_id = percent_identity(align_1, align_2)
        print(f"Percent identiy is {percent_id}")
    print(f"There is {alignment_counter} possible alignments")

#run 
global_alignment("ATTA", "ATTTTA")
