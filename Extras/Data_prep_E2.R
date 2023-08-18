library(igraph)
library(tidyr)


# Prepare the network for exercise
# Read network
Barabasi_Net <- read.csv("Data/Barabasi_CoV2_PPI.csv")
Barabasi_Net <- Barabasi_Net <- graph_from_data_frame(Barabasi_Net, directed = FALSE)


# Read list of human genes interacting with SARS-Cov-2
SarsCov_targets <- read.csv("Data/SARSCoV2_Targets.csv")

# Create subnetwork with the covid targets
Barabasi_SarsCov_Net <- induced_subgraph(Barabasi_Net, V(Barabasi_Net)[V(Barabasi_Net)$name %in% SarsCov_targets$EntrezID])
Barabasi_SarsCov_Net <- simplify(Barabasi_SarsCov_Net, remove.loops = TRUE) # remove loops
isolated <- names(which(degree(Barabasi_SarsCov_Net)==0)) # remove isolated nodes
Barabasi_SarsCov_Net <- delete.vertices(Barabasi_SarsCov_Net, isolated)
Barabasi_SarsCov_Net <- set_vertex_attr(Barabasi_SarsCov_Net, 
                                        name = "clustering_coefficient", 
                                        value = transitivity(Barabasi_SarsCov_Net, type = "local", isolates = "zero"))

Barabasi_SarsCov_Net <- induced_subgraph(Barabasi_SarsCov_Net, V(Barabasi_SarsCov_Net)[V(Barabasi_SarsCov_Net)$clustering_coefficient > 0.3])
selected_SarsCov_genes <- V(Barabasi_SarsCov_Net)$name


# Read disease genes association
disease_gene_links <- data.frame()
for(line in readLines("Data/Guney2016_GenesDisease.tsv")){
  tmp1 <- strsplit(line, "\t")
  tmp2 <- data.frame(disease = tmp1[[1]][2],
                     genes = paste(tmp1[[1]][3:length(tmp1[[1]])], collapse = ";"))
  disease_gene_links <- rbind(disease_gene_links, tmp2)
}
disease_gene_links <- separate_rows(disease_gene_links, "genes", sep = ";")

# Get disease genes from selected distant diseases from SARS-Cov-2
distant_diseases = c("peroxisomal disorders", "cardiomyopathy, hypertrophic", "anemia", "sarcoma")
dis_genes_select <- unique(disease_gene_links$genes[disease_gene_links$disease %in%  distant_diseases])

# Create network for exercise
Exercise_Net <- induced_subgraph(Barabasi_Net, 
                                 V(Barabasi_Net)[V(Barabasi_Net)$name %in% c(selected_SarsCov_genes, dis_genes_select)])
Exercise_Net <- as_data_frame(Exercise_Net, what = "edges")
write.csv(Exercise_Net, "Data/Exercise_PPI_Net.csv", row.names = FALSE)
