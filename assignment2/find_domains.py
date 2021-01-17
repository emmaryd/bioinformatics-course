# EMMA RYDHOLM ASSIGNMENT 2

import numpy as np
import matplotlib.pyplot as plt
import random

def read_pdb(filepath):
    """ Reads  the positions (x,y,z) from all carbon alpha atoms in a pdb file.
       Returns a list with the positions. """
    with open(filepath) as f:
    
        lines = f.readlines()
        relevant_lines = []
        coordinates = []

        for line in lines:
            line = line.split()
            if line[0] == "ATOM" and line[2] == "CA":
                relevant_lines.append(line)
    
    positions = []
    for line in relevant_lines:
        positions.append([float(line[6]), float(line[7]), float(line[8])])

    return positions

def calculate_distance(positions):
    """ Calculates distances between all positions in positions and return a matrix with these distances. """

    distance_array = np.zeros(shape=(len(positions), len(positions)))

    for i, pos_i in enumerate(positions):
        for j, pos_j in enumerate(positions):
            if i > j:
                distance_array[i][j] = distance_array[j][i] # becaus it is symmetric
            else:
                dist = ((pos_i[0] - pos_j[0])**2 + (pos_i[1] - pos_j[1])**2 + (pos_i[2] - pos_j[2])**2 )**0.5
                distance_array[i][j] = dist

    return distance_array

def plot_domains(distance_array, name, threshold, partion_index=None):
    """ Plot the distance dot plot and any domains """
    separation_array = distance_array < threshold

    dots_x = []
    dots_y = []
    dims = distance_array.shape

    for i in range(dims[0]):
        for j in range(dims[1]):
            if distance_array[i][j] < threshold:
                dots_x.append(i)
                dots_y.append(j)
    
    fig2 = plt.figure()
    ax = plt.axes()
    plt.scatter(dots_x, dots_y, c="b", marker='D', s=1.5, label="CA atoms")
    
    c=['r', 'g', 'cyan','purple']
    if partion_index:
        ax.set_title(f"Domains find by DOMAK on protein {name}")
        for num, i in enumerate(partion_index):
            linespace_1 = np.linspace(0, i)
            linespace_2 = np.linspace(i, dims[0])
            index_vector = np.ones(50)*i
            plt.plot(linespace_1, index_vector, c=c[num], linewidth=2, label=f'Split at residue {i}')
            plt.plot(index_vector, linespace_1, c=c[num], linewidth=2)

            plt.plot(linespace_2, index_vector, c=c[num], linewidth=2)
            plt.plot(index_vector, linespace_2, c=c[num], linewidth=2)

    else:
        ax.set_title(f"Distance map for protein {name}")

    plt.legend()
    plt.show()

def find_domains_DOMAK(distance_array, threshold=7):
    """ Find domains using DOMAK. """
    separation_array = distance_array < threshold
    dims = separation_array.shape
    max_score = 0
    best_partion = None
    for i in range(1, dims[0]-1):
        internal_a = np.sum(separation_array[0:i,0:i])
        internal_b = np.sum(separation_array[i:,i:])
        external = np.sum(separation_array[0:i,i:]) + np.sum(separation_array[i:,0:i])
        score = (internal_a / external) * (internal_b / external )
        if score > max_score:
            max_score = score
            best_partion = i

    return [best_partion]

def find_multiple_domains_DOMAK(distance_array, real_index_start, threshold=7, partion_index=[]):
    """ Find domains using DOMAK, this funktion can find more than just 2 domains compared to  find_domains_DOMAK"""

    min_domain_length = 40
    separation_array = distance_array < threshold
    dims = separation_array.shape
    max_score = 0
    best_partion = None
    np.fill_diagonal(separation_array, 0)
    for i in range(1, dims[0]-1):
        internal_a = np.sum(separation_array[0:i,0:i])
        internal_b = np.sum(separation_array[i:,i:])
        external = np.sum(separation_array[0:i,i:]) + np.sum(separation_array[i:,0:i])
        score = (internal_a / external) * (internal_b / external)
        if score > max_score:
            max_score = score
            best_partion = i
    
    domain_lengths = (best_partion + 1, dims[0] - best_partion +1)
    start_index = [real_index_start, real_index_start + best_partion]

    # Only is the domains found are longer than min_domain_lengt is the partioning valid
    if domain_lengths[0] > min_domain_length and domain_lengths[1] > min_domain_length:
        partion_index.append(best_partion + real_index_start)
        distance_arrays = [distance_array[0:best_partion, 0:best_partion],  distance_array[best_partion:, best_partion:]]
    
        for array, index in zip(distance_arrays, start_index):
            partion_index, N = find_multiple_domains_DOMAK(array, 
                                                           real_index_start=index,
                                                           threshold=threshold,
                                                           partion_index=partion_index, 
                                                           )

    n_domains = len(partion_index) + 1

    return partion_index, n_domains

def find_domains_STRUDL(distance_array, threshold=7, k=None):
    """ Find domains using STRDL and the Kernighan-Lin heuristic.  """

    separation_array = distance_array < threshold
    dims = separation_array.shape
    
    if k == None:
        k = int(dims[0]/2)

    # Initiate U and V and choose random u as starting point (u_seed)
    all_residues = [x for x in range(dims[0])]
    u_seed = random.choice(all_residues)
    U = [u_seed]
    all_residues.remove(u_seed)
    V = all_residues

    # add residues from V to U until len(U) == k
    for u in range(1, k):
        min_contacts = dims[0]
        for v in V:
            contacts = np.sum(separation_array[v,V])
            if contacts < min_contacts:
                min_contacts = contacts
                best_v = v
        U.append(best_v)
        V.remove(best_v)
    
    for _ in range(6): # Repeate 6 times
        best_U_V = None
        best_U_V_score = dims[0]**3 # set to a high number
        marked = []
        improvement = True
        
        # swap u and v until no improvements can be made or until all v is marked
        while improvement:
            max_links = 0
            u_swap = None
            # Chose the u that have most links to V
            for u in U:
                links_to_V = np.sum(separation_array[u,V])
                if links_to_V > max_links:
                    max_links = links_to_V
                    u_swap = u
            
            max_links = 0 
            v_swap = None   
            # Chose the v that have most links to U
            for v in V:
                if v not in marked:
                    links_to_U = np.sum(separation_array[U,v])
                    if links_to_U > max_links:
                        max_links = links_to_U
                        v_swap = v
                
            if u_swap == None or v_swap == None:
                improvement = False
            else:
                U.remove(u_swap)
                V.append(u_swap)

                V.remove(v_swap)
                U.append(v_swap)
                marked.append(u_swap)

            #calc U-V partioning score
            links_extern = 0
            for u in U:
                links_extern += np.sum(separation_array[u,V]) + np.sum(separation_array[V,u]) 
            if links_extern < best_U_V_score:
                best_U_V_score = links_extern
                best_U_V = (U.copy(), V.copy())
        U, V = best_U_V

    return best_U_V
  
"""   
# run tests
two_domain_files = ["1cdh.pdb", "2csn.pdb"]
for protein in two_domain_files:
    print(f"Protein = {protein}")
    positions = read_pdb(protein)
    distance_array = calculate_distance(positions)
    plot_domains(distance_array, name=protein, threshold=7)

    index = find_domains_DOMAK(distance_array, threshold=7)
    print(f"The best partition index using DOMAK is: {index}")
    U,V = find_domains_STRUDL(distance_array, threshold=7)
    print(f"The best partition index using STRUDL is: U = {np.sort(U)} and \n V = {np.sort(V)}")
    plot_domains(distance_array, name=protein, threshold=7, partion_index=index)

multiple_domain_files = ["4GAF_B.pdb", "1HZH_H.pdb"]
for protein in multiple_domain_files:
    print(f"Protein = {protein}")
    positions = read_pdb(protein)
    distance_array = calculate_distance(positions)
    index, N = find_multiple_domains_DOMAK(distance_array, real_index_start=0, threshold=7, partion_index=[])
    print(f"The best partition index using DOMAK is: {index} \n and there is {N} domains")
    plot_domains(distance_array, name=protein, threshold=7, partion_index=index)
"""
positions = read_pdb("2csn.pdb")
distance_array = calculate_distance(positions)
#plot_domains(distance_array, name=protein, threshold=7)

####
max_angle = 0
min_angle = 360
import math
angles = []
for protein in ["1cdh.pdb"]:#, "2csn.pdb"]:#,"4GAF_B.pdb", "1HZH_H.pdb"]:
    positions = read_pdb(protein)
    path = range(len(positions))
    for p,q,r in zip(path[0:-2], path[1:-1], path[2:]):
        #print(p,q,r)
        p = p-1
        q = q-1
        r = r-1

        p1 = positions[p]
        p2 = positions[q]
        p3 = positions[r]

        va = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]
        vb = [p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]]

        norm_a = ( va[0]**2 + va[1]**2 + va[2]**2 )**0.5
        norm_b = ( vb[0]**2 + vb[1]**2 + vb[2]**2 )**0.5

        angle = math.acos(np.dot(va,vb) / (norm_a*norm_b))
        if angle > max_angle:
            max_angle = angle
        elif angle < min_angle:
            min_angle = angle
        angles.append(angle)
print("mean angle: ",  np.sum(angles) / len(angles))
fig = plt.figure()
plt.hist(angles)
plt.show()

print('angle in range', min_angle, max_angle)
    