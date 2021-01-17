#EMMA RYDHOLM

import numpy as np
import math

UPPER_DIST = 3.95
LOWER_DIST = 3.65

def read_file(filepath):
    """ Reads the positions from file, assuming that the atom numbers starts at 1 and is aranged in order
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
    x_1,y_1,z_1 = pos_1
    x_2,y_2,z_2 = pos_2

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
        for j in range(n_points):
            pos_1, pos_2 = positions[i], positions[j]
            dist = distance(pos_1, pos_2)
            distances[i][j] = dist

    return distances

def find_path(distances):
    """Finds the best path that describes how the alpha-carbon atoms are arranged, 
    based on the assumption that all alpha-carbon atoms has distance at approximately 3.8Å of each other
    Args:
        distances: a matrix with the distances between the atom of interest
    Returns:
        path: The path represent to the most probable path of the atoms that makes the alpha-carbon chain
    """
    path = []
    for row in range(distances.shape[0]):
        if np.sum(row) == 1:
            current = row
            path.append(row)
            prev = None
            break

    end_of_sequence = False
    while end_of_sequence == False:
        potential_next = []
        for i, dist in enumerate(distances[current]):
            if dist > LOWER_DIST and dist < UPPER_DIST:
                if i != current and i != prev:
                    potential_next.append(i)
        if len(potential_next) == 1:
            path.append(potential_next[0])
        elif len(potential_next) > 1:
            most_probable = None
            least_deviation = 1
            for j in potential_next:
                deviation = np.absolute(distances[current][j] - 3.8)
                if deviation < least_deviation:
                    most_probable = j
                    least_deviation = deviation
            path.append(most_probable)

        elif len(potential_next) == 0:
            end_of_sequence = True # Set to true when the chain has no more connected atom
    
        current = path[-1]
        prev = path[-2]
    
    path = [x + 1 for x in path] # add one to all indices, since in the file the indices starts at 1
    print("Path:")
    for i in path:
        print(i)
    print("Number of alpha-carbon atoms in main chain: ", len(path))

    return path

#run    
files = ['data_q1.txt', 'test_q1.txt']
for f in files:
    print(f)
    positions = read_file(f)
    distances = distance_between_all(positions)
    path = find_path(distances)

# How to run the program: python3 main_chain.py