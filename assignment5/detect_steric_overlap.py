# EMMA RYDHOLM

# The approach taken to this problem is to compare the distance between 2 atoms. The atom radius is assumed to be 2 Å, 
# therefore an overlap is detected if two atoms are closer than 4 Å. The detected overlaps are stored as a list. 
# Since it is enought to find only 1 overlap "break" is used to end the inner loop if one  overlap is  detected. 
# This is to reduce the number of distance calculations.

import sys
import numpy as np

def read_pdb(filepath):
    """ Reads  the positions (x,y,z) from all atoms in a pdb file.
       Returns a list with tuples containing the atom number + information 
       and the positions of the atoms"""
    with open(filepath) as f:
    
        lines = f.readlines()
        relevant_lines = []

        for line in lines:
            line = line.split()
            if line[0] == "ATOM" or line[0] == "HETATM":
                relevant_lines.append(line)
    
    positions = []
    for line in relevant_lines:
        atom_number = (line[1], line[3], line[5], line[2])
        position = [float(line[6]), float(line[7]), float(line[8])]
        positions.append((atom_number, position))

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

def write_atoms_to_file(lines, file):
    """Writes line to the provided file. 
    Args: 
        lines: An iterable that contains the content that should be written to file. 
               Each value of the iterable will be printed on a new line and the values in  each value is separated  by a  tab.
        file: The file to where the lines  will be written.
    """
    f = open(file, "w")
    for line in lines:
        for i in line:
            f.write(i + "\t")
        f.write("\n")
    f.close()

def main(file_1, file_2, result_file):
    """ Finds the sterical overlaps and prints the atoms from file_1 that overlaps with 
    atoms from file_2 to result_file.
    Args:
        file_1: The file that gets compared with file_2. Should be a .pdb file.
        file_2: Reference file. Should be a .pdb file.
        result_file: The file to where the results is written  to, should be a .txt file
    """

    DIST_LIMIT = 4 # Limit for how close 2 atoms can be without any overlap

    positions_1 = read_pdb(file_1)
    positions_2 = read_pdb(file_2)

    counter_comparisions = 0
    counter_overlaps = 0
    # The atoms that overlaps are stored in a list
    overlaps = []

    for atom_1, pos_1 in positions_1:
        for atom_2, pos_2 in positions_2:
            dist = distance(pos_1, pos_2)
            counter_comparisions +=1
            if dist < DIST_LIMIT: # There is a overlap if the distance between the atoms is less than 2 atom radii
                counter_overlaps += 1
                overlaps.append(atom_1)
                # to reduce the number of comparisions the inner loop breaks if an overlap is detected
                # since it is enough to find overlap with one atom
                break

    write_atoms_to_file(overlaps,result_file)

    results_counting  = [f"Number of comparisions: {counter_comparisions}", 
                         f"Number of atom overlaps overlaps: {counter_overlaps}"]

    f = open(result_file, "a")
    for result in results_counting:
        f.write(result + "\n")
    f.close()

if __name__ == "__main__":
   main(*sys.argv[1:])

# HOW TO RUN:    python3 detect_steric_overlap.py file_1 file_2 result_file
# Compare 2csn.pdb with 1cdh.pdb:
#               python3 detect_steric_overlap.py 2csn.pdb 1cdh.pdb 2csn-1cdh.txt
# Compare  1cdh.pdb with 2csn.pdb:
#               python3 detect_steric_overlap.py 1cdh.pdb 2csn.pdb 1cdh-2csn.txt
