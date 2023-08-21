# Functions for calculating Barabasi network proximities for drug combinations
# Refer: https://doi.org/10.1038/s41467-019-09186-x





calc_single_set_distance <- function(gene_network, geneSet){
  
  # Check if the input geneSet is in the network
  geneSet <- geneSet[geneSet %in% V(gene_network)$name]
  if(length(geneSet) == 0){
    print("Input geneSet not in network")
    return(NA)
  }else{
    # distances() creates a matrix of shortest path lengths with 'v' in rows and 'to' in columns
    distance_matrix <- distances(graph = gene_network,
                                 v = V(gene_network)[V(gene_network)$name %in% geneSet],
                                 to = V(gene_network)[V(gene_network)$name %in% geneSet],
                                 algorithm = "dijkstra")
    
    if(length(unique(geneSet)) == 1){ 
      return(0) # Return 0 or else will be returned NA
    }else{
      distance_matrix[is.infinite(distance_matrix)] <- NA 
      distance_matrix[distance_matrix == 0] <- NA # Since distance from the same gene will be 0
      
      # minimum distances from each gene in geneset
      geneSet_min_distance <- data.frame()
      for(gene in row.names(distance_matrix)){
        tmp <- data.frame(gene = gene,
                          min_distance = min((distance_matrix[gene,]), na.rm = TRUE))
        geneSet_min_distance <- rbind(geneSet_min_distance, tmp)
      }
      geneSet_min_distance[is.infinite(geneSet_min_distance$min_distance), "min_distance"] <- NA
      
      # Average of minimum distances from each gene in geneset
      geneSet_mean_distance <- mean(geneSet_min_distance$min_distance, na.rm = TRUE)
      return(geneSet_mean_distance)
    }
  }
}





calc_set_pair_distances <- function(gene_network, geneSet1, geneSet2){
  
  # Check if the input geneSet is in the network
  geneSet1 <- geneSet1[geneSet1 %in% V(gene_network)$name]
  geneSet2 <- geneSet2[geneSet2 %in% V(gene_network)$name]
  
  if(length(geneSet1) == 0 || length(geneSet2) == 0){
    print("One of input geneSet not in network")
    return(NA)
  }else{
    # distances() creates a matrix of shortest path lengths with geneSet1 in rows and geneSet2 in columns
    distance_matrix <- distances(graph = gene_network,
                                 v = V(gene_network)[V(gene_network)$name %in% geneSet1],
                                 to = V(gene_network)[V(gene_network)$name %in% geneSet2],
                                 algorithm = "dijkstra")
    
    distance_matrix[is.infinite(distance_matrix)] <- NA
    distance_matrix[distance_matrix == 0] <- NA # Since distance from the same gene will be 0
    
    # minimum distances from each gene in geneSet1
    geneSet_min_distance <- data.frame()
    for(gene in row.names(distance_matrix)){
      tmp <- data.frame(gene = gene,
                        min_distance = min((distance_matrix[gene,]), na.rm = TRUE))
      geneSet_min_distance <- rbind(geneSet_min_distance, tmp)
    }
    
    # minimum distances from each gene in geneSet2
    for(gene in colnames(distance_matrix)){
      tmp <- data.frame(gene = gene,
                        min_distance = min((distance_matrix[,gene]), na.rm = TRUE))
      geneSet_min_distance <- rbind(geneSet_min_distance, tmp)
    }
    geneSet_min_distance[is.infinite(geneSet_min_distance$min_distance), "min_distance"] <- NA
    
    # Average of minimum distances from each gene in geneSet1 and geneSet2
    geneSet_mean_distance <- mean(geneSet_min_distance$min_distance, na.rm = TRUE)
    return(geneSet_mean_distance)
  }
}





Barabasi_proximity_separation <- function(gene_network, geneSet1, geneSet2){
  
  # Check if input genesets are valid
  if(length(geneSet1) == 0 || length(geneSet2) == 0){
    print("One geneSet is empty")
    return(NA)
  }else{
    # Check if the input geneSets are in the network
    geneSet1 <- geneSet1[geneSet1 %in% V(gene_network)$name]
    geneSet2 <- geneSet2[geneSet2 %in% V(gene_network)$name]
    
    if(length(geneSet1) == 0 || length(geneSet2) == 0){
      print("One geneSet is not in the network")
      return(NA)
    }else{
      
      dist_geneSet1 <- calc_single_set_distance(gene_network, geneSet1)
      dist_geneSet2 <- calc_single_set_distance(gene_network, geneSet2)
      dist_geneSet1_geneSet2 <- calc_set_pair_distances(gene_network, geneSet1, geneSet2)
      separation <-  dist_geneSet1_geneSet2  - (dist_geneSet1 + dist_geneSet2)/2
      
      return(separation)
    }
  }
}





Barabasi_proximity_closest <- function(gene_network, geneSet1, geneSet2){
  
  # Check if input genesets are valid
  if(length(geneSet1) == 0 || length(geneSet2) == 0){
    print("One geneSet is empty")
    return(NA)
  }else{
    # Check if the input geneSets are in the network
    geneSet1 <- geneSet1[geneSet1 %in% V(gene_network)$name]
    geneSet2 <- geneSet2[geneSet2 %in% V(gene_network)$name]
    
    if(length(geneSet1) == 0 || length(geneSet2) == 0){
      print("One geneSet is not in the network")
      return(NA)
    }else{
      # distances() creates a matrix of shortest path lengths with geneSet1 in rows and geneSet2 in columns
      distance_matrix <- distances(graph = gene_network,
                                   v = V(gene_network)[V(gene_network)$name %in% geneSet1],
                                   to = V(gene_network)[V(gene_network)$name %in% geneSet2],
                                   algorithm = "dijkstra")
      
      
      distance_matrix[is.infinite(distance_matrix)] <- NA
      distance_matrix[distance_matrix == 0] <- NA # Since distance from the same gene will be 0
      
      
      # sum of minimum distances from each geneSet1 genes to the geneSet2 genes
      geneSet1_min_distance <- data.frame()
      for(gene in row.names(distance_matrix)){
        tmp <- data.frame(gene = gene,
                          min_distance = min(distance_matrix[gene, ], na.rm = TRUE))
        geneSet1_min_distance <- rbind(geneSet1_min_distance, tmp)
      }
      geneSet1_min_distance_sum <- sum(geneSet1_min_distance[, "min_distance"])
      
      
      
      # sum of minimum distances from each geneSet2 genes to the geneSet1 genes
      geneSet2_min_distance <- data.frame()
      for(gene in colnames(distance_matrix)){
        tmp <- data.frame(gene = gene,
                          min_distance = min(distance_matrix[,gene], na.rm = TRUE))
        geneSet2_min_distance <- rbind(geneSet2_min_distance, tmp)
      }
      geneSet2_min_distance_sum <- sum(geneSet2_min_distance[, "min_distance"])
      
      proximity_closest <- (geneSet1_min_distance_sum + geneSet2_min_distance_sum)/(length(geneSet1) + length(geneSet2))
      return(proximity_closest)
    }
  }
}





Barabasi_proximity_shortest <- function(gene_network, geneSet1, geneSet2){
  
  # Check if input genesets are valid
  if(length(geneSet1) == 0 || length(geneSet2) == 0){
    print("One geneSet is empty")
    return(NA)
  }else{
    # Check if the input geneSets are in the network
    geneSet1 <- geneSet1[geneSet1 %in% V(gene_network)$name]
    geneSet2 <- geneSet2[geneSet2 %in% V(gene_network)$name]
    
    if(length(geneSet1) == 0 || length(geneSet2) == 0){
      print("One geneSet is not in the network")
      return(NA)
    }else{
      
      # distances() creates a matrix of shortest path lengths with geneSet1 in rows and geneSet2 in columns
      distance_matrix <- distances(graph = gene_network,
                                   v = V(gene_network)[V(gene_network)$name %in% geneSet1],
                                   to = V(gene_network)[V(gene_network)$name %in% geneSet2],
                                   algorithm = "dijkstra")
      
      distance_matrix[is.infinite(distance_matrix)] <- NA
      # distance_matrix[distance_matrix == 0] <- NA # Since distance from the same gene will be 0
      
      proximity_shortest <- sum(distance_matrix, na.rm = TRUE)/(length(geneSet1) * length(geneSet2))
      
      return(proximity_shortest)
    }
  }
}





get_topological_centre <- function(gene_network, geneSet){
  
  if(length(geneSet) == 0){
    print("Input geneSet is empty")
    return(NA)
  }else{
    # Check if the input geneSet are in the network
    geneSet <- geneSet[geneSet %in% V(gene_network)$name]
    
    if(length(geneSet) == 0){
      print("Input geneSet is not in the network")
      return(NA)
    }else{
      # distances() creates a matrix of shortest path lengths with geneSet1 in rows and geneSet2 in columns
      distance_matrix <- distances(graph = gene_network,
                                   v = V(gene_network)[V(gene_network)$name %in% geneSet],
                                   to = V(gene_network)[V(gene_network)$name %in% geneSet],
                                   algorithm = "dijkstra")
      
      
      # Calculate sum of distances of each gene in geneSet to all other genes in the geneSet
      geneSet_sum_dist <- rowSums(distance_matrix, na.rm = TRUE)
      
      topological_centre <- names(geneSet_sum_dist[geneSet_sum_dist == min(geneSet_sum_dist)])
      
      return(topological_centre)
    }
  }
}





Barabasi_proximity_centre <- function(gene_network, geneSet1, geneSet2){
  
  # Check if input genesets are valid
  if(length(geneSet1) == 0 || length(geneSet2) == 0){
    print("One geneSet is empty")
    return(NA)
  }else{
    # Check if the input geneSets are in the network
    geneSet1 <- geneSet1[geneSet1 %in% V(gene_network)$name]
    geneSet2 <- geneSet2[geneSet2 %in% V(gene_network)$name]
    
    if(length(geneSet1) == 0 || length(geneSet2) == 0){
      print("One geneSet is not in the network")
      return(NA)
    }else{
      geneSet1_centre <- get_topological_centre(gene_network, geneSet1)
      geneSet2_centre <- get_topological_centre(gene_network, geneSet2)
      
      distance_matrix <- distances(graph = gene_network,
                                   v = V(gene_network)[V(gene_network)$name %in% geneSet1_centre],
                                   to = V(gene_network)[V(gene_network)$name %in% geneSet2_centre],
                                   algorithm = "dijkstra")
      
      if(nrow(distance_matrix) > 1 && ncol(distance_matrix) == 1){
        proximity_centre <- unname(colMeans(distance_matrix, na.rm = TRUE))
      }else if(nrow(distance_matrix) == 1 && ncol(distance_matrix) > 1){
        proximity_centre <- unname(rowMeans(distance_matrix, na.rm = TRUE))
      }else if(nrow(distance_matrix) > 1 && ncol(distance_matrix) > 1){
        proximity_centre <- sum(distance_matrix, na.rm = TRUE)/(nrow(distance_matrix) + ncol(distance_matrix))
      }else if(nrow(distance_matrix) == 1 && ncol(distance_matrix) == 1){
        proximity_centre <- distance_matrix[1,1]
      }
      return(proximity_centre)
    }
  }
}





Barabasi_proximity_kernel <- function(gene_network, geneSet1, geneSet2){
  
  # Check if input genesets are valid
  if(length(geneSet1) == 0 || length(geneSet2) == 0){
    print("One geneSet is empty")
    return(NA)
  }else{
    # Check if the input geneSets are in the network
    geneSet1 <- geneSet1[geneSet1 %in% V(gene_network)$name]
    geneSet2 <- geneSet2[geneSet2 %in% V(gene_network)$name]
    
    if(length(geneSet1) == 0 || length(geneSet2) == 0){
      print("One geneSet is not in the network")
      return(NA)
    }else{
      
      # distances() creates a matrix of shortest path lengths with geneSet1 in rows and geneSet2 in columns
      distance_matrix <- distances(graph = gene_network,
                                   v = V(gene_network)[V(gene_network)$name %in% geneSet1],
                                   to = V(gene_network)[V(gene_network)$name %in% geneSet2],
                                   algorithm = "dijkstra")
      
      
      distance_matrix[is.infinite(distance_matrix)] <- NA
      # distance_matrix[distance_matrix == 0] <- NA # Since distance from the same gene will be 0
      
      
      transform_matrix <- function(x, size){
        y = exp(-(x+1))/size
        return(y)
      }
      
      geneSet1_transformed_dist <- data.frame()
      distance_matrix_transformed <- apply(distance_matrix, c(1,2), function(z)transform_matrix(z, length(geneSet1)))
      for(gene in row.names(distance_matrix_transformed)){
        transformed_dist <- log(sum(distance_matrix_transformed[gene,], na.rm = TRUE))
        tmp <- data.frame(gene = gene, transformed_dist = transformed_dist)
        geneSet1_transformed_dist <- rbind(geneSet1_transformed_dist, tmp)
      } 
      geneSet1_transformed_dist_sum <- sum(geneSet1_transformed_dist[, "transformed_dist"])
      
      geneSet2_transformed_dist <- data.frame()
      distance_matrix_transformed <- apply(distance_matrix, c(1,2), function(z)transform_matrix(z, length(geneSet2)))
      for(gene in colnames(distance_matrix_transformed)){
        transformed_dist <- log(sum(distance_matrix_transformed[,gene], na.rm = TRUE))
        tmp <- data.frame(gene = gene, transformed_dist = transformed_dist)
        geneSet2_transformed_dist <- rbind(geneSet2_transformed_dist, tmp)
      } 
      geneSet2_transformed_dist_sum <- sum(geneSet2_transformed_dist[, "transformed_dist"])
      
      proximity_kernel <- (-(geneSet1_transformed_dist_sum + geneSet2_transformed_dist_sum)/(length(geneSet1) + length(geneSet2)))
      
      return(proximity_kernel)
    }
  }
}
