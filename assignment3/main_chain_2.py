#EMMA RYDHOLM

import numpy as np
import pandas as pd
import math 

UPPER_DIST = 3.95
LOWER_DIST = 3.65

UPPER_ANGLE = 1.65
LOWER_ANGLE = 0.55

def read_file(filepath):
    """ Reads the positions from file
    Args: 
        filepath: string that contains the filepath
    """
    positions = []

    with open(filepath) as f:
        lines = f.readlines()
        for line in lines:
            line = line.split()
            atom_position = (float(line[1]),float(line[2]),float(line[3]),)
            positions.append(atom_position)
    
    return positions

def distance(pos_1,pos_2):
    """Calculate the distance between 2 points.
    Args: 
        pos_1: describes the position 1 in 3 dimensions
        pos_2: describes the position 2 in 3 dimensions
    Returns: 
        distance is the euclidean distance between pos_1 and pos_2
    """
    x_1, y_1, z_1 = pos_1
    x_2, y_2, z_2 = pos_2

    distance = ( (x_2 - x_1)**2 + (y_2 - y_1)**2 + (z_2 - z_1)**2 ) **0.5

    return distance

def distance_between_all(positions):
    """Calculates the distances between all points given in positions
    Args:   
        positions: an itetable containing the positions for each atom
    Returns:
        distances: an array with the distances between all atoms
    """
    n_points = len(positions)
    distances = np.zeros((n_points, n_points))
    for i in range(n_points):
        for j in range(i, n_points):
            pos_1, pos_2 = positions[i], positions[j]
            dist = distance(pos_1, pos_2)
            distances[i][j] = dist
            distances[j][i] = dist

    return distances

def pseudo_valenc_angle(p1, p2, p3):
    """Calculate the pseudo_valenc_angle between 3 positions. 
    Args:
        p1: position 1
        p2: position 2 (middle point)
        p3: position 3 
    Returns:
        angle: the angle betweens
    """
    va = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]
    vb = [p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]]
   
    norm_a = ( va[0]**2 + va[1]**2 + va[2]**2 )**0.5
    norm_b = ( vb[0]**2 + vb[1]**2 + vb[2]**2 )**0.5

    angle = math.acos(np.dot(va,vb) / (norm_a*norm_b))

    return angle

def search_path(path, all_paths, positions, distances):
    """Finds a path that satisfies the distanse and angle conditions set by UPPER_DIST, LOWER_DIST and UPPER_ANGLE, LOWER_ANGLE
    Args:
        path: Should be a list, containing at least 2 atoms but could be longer
        all_path: A list where all valid paths are stored
        positions: All atom positions in order
        distances: A numpy array zero on all positions that does not satisfy the distance condition and the distance on 
                   all positions that do satisfy the condition 
    Returns: an uppdated all_path list 
    """
    prev, current = path[-2], path[-1]
    end_of_sequence = False
    
    while end_of_sequence == False:
            potential_next = []
            potential_next = np.nonzero(distances[current])[0]
            potential_next = [i for i in potential_next if i != prev]  

            potential_next_ = []

            # Sort out the atoms that does not satisfy the angle requirement
            for atom in potential_next:
                angle = pseudo_valenc_angle(positions[prev], positions[current], positions[atom])
                if angle > LOWER_ANGLE and angle < UPPER_ANGLE:
                    potential_next_.append(atom)   
            
            if len(potential_next_) == 1:
                path.append(potential_next_[0])

            elif len(potential_next_) > 1:  
                # Search the path for all next atoms in potential next              
                for atom in potential_next_[1:]:
                    opt_path = path.copy()
                    opt_path.append(atom)
                    all_path = search_path(opt_path, all_paths, positions, distances)
                path.append(potential_next_[0])                  
                        
            elif len(potential_next_) == 0:
                end_of_sequence = True # Set to true when the chain has no more connected atom
            
            # Check for any loops, if the lastly added atom is in any other place in tha chain, the search ends
            if path[-1] in path[0:-2]:
                end_of_sequence = True
                path.pop()
            current = path[-1]
            prev = path[-2]
    
    # append to all_path only if the chain has more than 2 atoms
    if len(path) > 2:
        all_paths.append(path)

    return all_paths

def find_path(positions):
    """Finds a potential path that describes how the alpha-carbon atoms are arranged, 
    based on the assumption that all alpha-carbon atoms has distance at approximately 3.8 Å of each other
    and that the pseudo valence angle have a certain value specified by UPPER_ANGLE and LOWER_ANGLE
    Args:
        distances: a matrix with the distances between the atom of interest
    Returns:
        path: The path represent to the longest probable path of the atoms that makes the alpha-carbon chain
    """
    distances = distance_between_all(positions)

    distances = pd.DataFrame(distances)
    distances = distances[distances > LOWER_DIST]
    distances = distances[distances < UPPER_DIST]
    distances = distances.fillna(0)
    distances = np.array(distances)

    potential_start = []
    for row in range(distances.shape[0]):
        if np.sum(distances[row]) > LOWER_DIST  and np.sum(distances[row]) < UPPER_DIST :
            nonzero_index = np.nonzero(distances[row])
            prev = row
            current = nonzero_index[0][0]
            potential_start.append([prev, current])

    all_paths = []
    for atom_pair in potential_start:
        path = atom_pair

        all_paths = search_path(path, all_paths, positions, distances)

    all_ca =[]
    longest_path_len = 0
    for path in all_paths:
        all_ca.extend(path)
        if len(path) > longest_path_len:
            longest_path = path
            longest_path_len = len(path)
        
    #fix correct index
    longest_path = [x + 1 for x in longest_path] # add one to all indices, since in the file the indices starts at 1
    print("Path:")
    for i in longest_path:
        print(i)
    print("Number of alpha-carbon atoms in main chain: ", longest_path_len)
    return longest_path

#run    
positions = read_file('data_q2.txt')
path = find_path(positions)

# How to run the program: python3 main_chain_2.py
