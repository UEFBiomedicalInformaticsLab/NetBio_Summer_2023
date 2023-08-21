set.seed(5081)



# Data preparation (Exercise 1)


# Load libraries
library(tidyverse)
library(igraph)


# Extract the human drug targets from DrugBank
DrugBank_Targets_idMap <- read.csv("Data_prep/targets_polypeptides_ext_id.csv", header = TRUE)
DrugBank_Targets_idMap <- reshape(as.data.frame(DrugBank_Targets_idMap), idvar = "parent_key", timevar = "resource", direction = "wide")
colnames(DrugBank_Targets_idMap) <- gsub("identifier.", "", colnames(DrugBank_Targets_idMap))
DrugBank_Targets_idMap <- DrugBank_Targets_idMap[grep("_HUMAN$", x = DrugBank_Targets_idMap$`UniProt Accession`),]


# Process STRING DB

# Read protein Id mappings
String_proteins <- read.table("Data_prep/9606.protein.aliases.v11.5.txt.gz", fill = TRUE, sep = "\t")
colnames(String_proteins) <- c("string_protein_id", "alias", "source")
String_proteins <- String_proteins[String_proteins$source == "Ensembl_EntrezGene",]

# Read all PPI in String DB along with scores
String_ppi <- read.table("Data_prep/9606.protein.links.detailed.v11.5.txt.gz", header = TRUE)


# Filter all interactions with scores greater than 0 and sourced from database/experimets
String_ppi <- String_ppi[(String_ppi$combined_score > quantile(String_ppi$combined_score, 0.95)), ] #

String_ppi$from_geneSymbol <- String_proteins$alias[match(String_ppi$protein1, String_proteins$string_protein_id)]
String_ppi$to_geneSymbol <- String_proteins$alias[match(String_ppi$protein2, String_proteins$string_protein_id)]

String_ppi_Net <- na.exclude(String_ppi[, c("from_geneSymbol", "to_geneSymbol", "combined_score")])


drug_target_Net <- String_ppi_Net[(String_ppi_Net$from_geneSymbol %in% DrugBank_Targets_idMap$GenAtlas &
                                         String_ppi_Net$to_geneSymbol %in% DrugBank_Targets_idMap$GenAtlas),]

write.csv(drug_target_Net, "Data/Human_drug_target_net.csv", row.names = FALSE)




# Create drugs

# Extract all drug target interactions in human
DrugBank_Drug_Target <- read.csv("Data_prep/targets.csv", header = TRUE)
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$organism == "Humans", ]
colnames(DrugBank_Drug_Target)[colnames(DrugBank_Drug_Target) == "parent_key"] <- "drug_id"
DrugBank_Drug_Target$drug_target <- DrugBank_Targets_idMap$GenAtlas[match(DrugBank_Drug_Target$id, DrugBank_Targets_idMap$parent_key)]


# Keep only approved drugs
DrugBank_drug_groups <- read.csv("Data_prep/drug_groups.csv", header = TRUE, check.names = FALSE)
DrugBank_drug_groups <- DrugBank_drug_groups[DrugBank_drug_groups$group == "approved", ]
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$drug_id %in% DrugBank_drug_groups$`drugbank-id`, ]

# Keep only drugs whose targets are in the network
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$drug_target %in% c(drug_target_Net$from_geneSymbol, drug_target_Net$to_geneSymbol),]


# Add the drug names
DrugBank_Drugs <- read.csv("Data_prep/drug.csv", header = TRUE)
DrugBank_Drug_Target$drug_name <- DrugBank_Drugs$name[match(DrugBank_Drug_Target$drug_id, DrugBank_Drugs$primary_key)]
DrugBank_Drugs <- DrugBank_Drugs[DrugBank_Drugs$type == "small molecule", ]
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$drug_id %in% DrugBank_Drugs$primary_key, ]


DrugBank_Drug_Target <- DrugBank_Drug_Target[, c("drug_id", "drug_name", "drug_target", "name")]
colnames(DrugBank_Drug_Target) <- c("drug_id", "drug_name", "drug_target", "drug_target_name")

# Remove drugs that directly target ACE2
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$drug_target != "ACE2", ]

# Get drugs with only certain categories
DrugBank_drug_categories <- read.csv("Data_prep/drug_categories.csv", header = TRUE)
DrugBank_drug_categories <- DrugBank_drug_categories[DrugBank_drug_categories$category %in% c("Antiviral Agents", "Antiprotozoals", "Antiparasitic Agents"), ]
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$drug_id %in% DrugBank_drug_categories$parent_key, ]


# Filter to keep only the most distant targets from ACE2
# If there are multiple targets with same distance, they are removed
# Doing this to reduce the number of drugs for exercise
network <- graph_from_data_frame(drug_target_Net, directed = FALSE)

cc <- components(network)
largest_components <- V(network)[cc$membership == which.max(cc$csize)] # remove connected components
network <- induced_subgraph(network, largest_components)

distance_table <- as.data.frame(distances(network, v = "ACE2", to = V(network)$name))


final_df <- data.frame()
for(drug in unique(DrugBank_Drug_Target$drug_id)){
  tmp1 <- DrugBank_Drug_Target[DrugBank_Drug_Target$drug_id == drug, ]
  
  if(nrow(tmp1) == 1){
    if(tmp1$drug_target %in% V(network)$name){
      final_df <- rbind(final_df, tmp1)
    }
  }else{
    tmp2 <- distance_table[names(distance_table) %in% tmp1$drug_target]
    tmp1 <- tmp1[tmp1$drug_target %in% names(tmp2[which(tmp2 == max(tmp2))]),]
    if(nrow(tmp1) == 1){
      final_df <- rbind(final_df, tmp1)
    }
  }
}



write.csv(final_df, "Data/Drugs.csv", row.names = FALSE)

# "Antiviral Agents", "Antiprotozoals", "Antiparasitic Agents", "Antivirals for Systemic Use",
#"Experimental Unapproved Treatments for COVID-19", "Approved Treatments for COVID-19"


print(warnings())