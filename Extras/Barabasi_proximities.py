# Functions to calculate Barabasi proximities
# Refer: https://www.nature.com/articles/ncomms10331
# Defined: AG


import numpy as np
import pandas as pd
import networkx as nx

# Barabasi proximity: closest
def proximity_closest(network, geneSet1, geneSet2): 
    # Check if the genes in geneSet1 and geneSet2 are present in the network
    geneSet1 = geneSet1 & network.nodes
    geneSet2 = geneSet2 & network.nodes
    
    # Calculate all distances between geneset1 and geneset2
    distance_matrix = pd.DataFrame(data = np.nan, index = geneSet1, columns  = geneSet2)
    for gene1 in geneSet1:
        for gene2 in geneSet2:
            if gene1 != gene2:
                distance_matrix.loc[gene1, gene2] = nx.shortest_path_length(network, gene1, gene2)

    # Find minimum distance of each gene in geneSet1 to geneSet2
    geneSet1_min_dist = pd.DataFrame()
    for gene in geneSet1:
        min_dist = min(list(distance_matrix.loc[gene,:].dropna()))
        tmp = pd.DataFrame({"Gene" : gene, "min_dist" : min_dist}, index = [0])
        geneSet1_min_dist = pd.concat([geneSet1_min_dist, tmp]).reset_index(drop = True)

    # Calculate the average of the ninimum distances
    closest_dist = geneSet1_min_dist["min_dist"].mean()
    return closest_dist


# Barabasi proximity: shortest
def proximity_shortest(network, geneSet1, geneSet2):
    # Check if the genes in geneSet1 and geneSet2 are present in the network
    geneSet1 = geneSet1 & network.nodes
    geneSet2 = geneSet2 & network.nodes

    # Calculate all distances between geneset1 and geneset2
    distance_matrix = pd.DataFrame(data = np.nan, index = geneSet1, columns  = geneSet2)
    for gene1 in geneSet1:
        for gene2 in geneSet2:
            if gene1 != gene2:
                distance_matrix.loc[gene1, gene2] = nx.shortest_path_length(network, gene1, gene2)

    # Find average distance of each gene in geneSet1 to geneSet2
    geneSet1_avg_dist = pd.DataFrame()
    for gene in geneSet1:
        avg_dist = distance_matrix.loc[gene,:].dropna().mean()
        tmp = pd.DataFrame({"Gene" : gene, "avg_dist" : avg_dist}, index = [0])
        geneSet1_avg_dist = pd.concat([geneSet1_avg_dist, tmp]).reset_index(drop = True)

    # Calculate the average of the average distances
    average_dist = geneSet1_avg_dist["avg_dist"].mean()
    return average_dist


# Get topological centre
def get_topological_centre(network, geneSet):
    # Check if the genes in geneSet are present in the network
    geneSet = geneSet & network.nodes

    # Calculate all distances within  geneset1 
    distance_matrix = pd.DataFrame(data = np.nan, index = geneSet, columns  = geneSet)
    for gene1 in geneSet:
        for gene2 in geneSet:
            if gene1 != gene2:
                distance_matrix.loc[gene1, gene2] = nx.shortest_path_length(network, gene1, gene2)

    # Find average distance within each gene in geneSet
    geneSet_sum_dist = pd.DataFrame()
    for gene in geneSet:
        sum_dist = distance_matrix.loc[gene,:].dropna().sum()
        tmp = pd.DataFrame({"Gene" : gene, "sum_dist" : sum_dist}, index = [0])
        geneSet_sum_dist = pd.concat([geneSet_sum_dist, tmp]).reset_index(drop = True)

    # Get the central genes
    min_value = geneSet_sum_dist.loc[geneSet_sum_dist["sum_dist"].argmin(), "sum_dist"]
    topological_centre = set(geneSet_sum_dist.loc[geneSet_sum_dist["sum_dist"] == min_value, "Gene"])
    return topological_centre

# Barabasi proximity: centre
def proximity_centre(network, topological_centre, geneSet2):
    
    # Check if the genes in topological_centre and geneSet2 are present in the network
    geneSet1 = topological_centre & network.nodes
    geneSet2 = geneSet2 & network.nodes

    # Calculate all distances between geneset1 and geneset2
    distance_matrix = pd.DataFrame(data = np.nan, index = geneSet1, columns  = geneSet2)
    for gene1 in geneSet1:
        for gene2 in geneSet2:
            if gene1 != gene2:
                distance_matrix.loc[gene1, gene2] = nx.shortest_path_length(network, gene1, gene2)

    # Find average distance of each gene in geneSet1 to geneSet2
    geneSet1_avg_dist = pd.DataFrame()
    for gene in geneSet1:
        avg_dist = distance_matrix.loc[gene,:].dropna().mean()
        tmp = pd.DataFrame({"Gene" : gene, "avg_dist" : avg_dist}, index = [0])
        geneSet1_avg_dist = pd.concat([geneSet1_avg_dist, tmp]).reset_index(drop = True)
    
    dist_to_centre = geneSet1_avg_dist["avg_dist"].mean()
    return dist_to_centre


# Barabasi proximity: kernel
def proximity_kernel(network, geneSet1, geneSet2):
    # Check if the genes in geneSet1 and geneSet2 are present in the network
    geneSet1 = geneSet1 & network.nodes
    geneSet2 = geneSet2 & network.nodes

    # Calculate all distances between geneset1 and geneset2
    distance_matrix = pd.DataFrame(data = np.nan, index = geneSet1, columns  = geneSet2)
    for gene1 in geneSet1:
        for gene2 in geneSet2:
            if gene1 != gene2:
                distance_matrix.loc[gene1, gene2] = nx.shortest_path_length(network, gene1, gene2)

    import math
    def transform_df(x):
        return (math.e**(-(x+1)))/len(geneSet1)

    # Apply kernel transformation
    transformed_matrix = distance_matrix.apply(np.vectorize(transform_df))

    # 
    geneSet1_transformed_dist = pd.DataFrame()
    for gene in geneSet1:
        transformed_dist = math.log(transformed_matrix.loc[gene,:].dropna().sum())
        tmp = pd.DataFrame({"Gene" : gene, "transformed_dist" : transformed_dist}, index = [0])
        geneSet1_transformed_dist = pd.concat([geneSet1_transformed_dist, tmp]).reset_index(drop = True)

    # Calculate the average of the transformed values
    kernel_dist = -geneSet1_transformed_dist["transformed_dist"].mean()
    return kernel_dist