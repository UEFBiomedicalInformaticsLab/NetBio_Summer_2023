set.seed(5081)


# Script to prepare the data for exercise 2



# Load libraries
library(igraph)
library(tidyverse)
library(sparklyr)
library(sparklyr.nested)
source("Extras/Functions_Barabasi_metrics.R")





# Create the PPI

# Read protein Id mappings
String_proteins <- read.table("Data_prep/9606.protein.aliases.v11.5.txt.gz", fill = TRUE, sep = "\t")
colnames(String_proteins) <- c("string_protein_id", "alias", "source")
String_proteins <- String_proteins[String_proteins$source == "Ensembl_EntrezGene",]

# Read all PPI in String DB along with scores
String_ppi <- read.table("Data_prep/9606.protein.links.detailed.v11.5.txt.gz", header = TRUE)

# Filter all interactions 
String_ppi <- String_ppi[(String_ppi$combined_score > quantile(String_ppi$combined_score, 0.75)), ] #

String_ppi$from_geneSymbol <- String_proteins$alias[match(String_ppi$protein1, String_proteins$string_protein_id)]
String_ppi$to_geneSymbol <- String_proteins$alias[match(String_ppi$protein2, String_proteins$string_protein_id)]

String_ppi_Net <- na.exclude(String_ppi[, c("from_geneSymbol", "to_geneSymbol", "combined_score")])
String_ppi_Net <- graph_from_data_frame(d = String_ppi_Net, directed = FALSE)

# Get disease genes
if(!dir.exists("Data_prep/")){dir.create("Data_prep/", recursive = TRUE)}
if(!dir.exists("Data_prep/associationByDatatypeDirect")){
  system("wget --recursive --no-parent --no-host-directories -P Data_prep/ --cut-dirs 8 ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/21.04/output/etl/parquet/associationByDatatypeDirect
")
}

sc <- spark_connect(master = "local")
OpenTargets_Target_Disease <- spark_read_parquet(sc, "OpenTargets_Target_Disease", "Data_prep/associationByDatatypeDirect/")

# List columns in the tibble
columns <- OpenTargets_Target_Disease %>%
  sdf_schema() %>%
  lapply(function(x) do.call(tibble, x)) %>%
  bind_rows()

OpenTargets_Target_Disease_lit <- OpenTargets_Target_Disease %>%
  filter(datatypeId == "literature")

OpenTargets_Target_Disease_lit %>%
  collect()  
OpenTargets_Target_Disease_lit <- as.data.frame(OpenTargets_Target_Disease_lit)
OpenTargets_Target_Disease_lit <- OpenTargets_Target_Disease_lit[OpenTargets_Target_Disease_lit$datatypeHarmonicScore > quantile(OpenTargets_Target_Disease_lit$datatypeHarmonicScore, 0.75),]

OpenTargets_Target_Disease_lit <- OpenTargets_Target_Disease_lit[, c("targetSymbol", "diseaseId", "diseaseLabel")]
diseaseID_2_diseaseLabel <- unique(OpenTargets_Target_Disease_lit[,c("diseaseId", "diseaseLabel")])

OpenTargets_Disease2Gene_lit_lib <- list()
for(i in unique(OpenTargets_Target_Disease_lit$diseaseId)) {
  OpenTargets_Disease2Gene_lit_lib[[i]] <- OpenTargets_Target_Disease_lit[OpenTargets_Target_Disease_lit$diseaseId == i,]$targetSymbol
}
for(i in 1:length(OpenTargets_Disease2Gene_lit_lib)){
  diseaseName <- paste0(diseaseID_2_diseaseLabel[diseaseID_2_diseaseLabel$diseaseId == names(OpenTargets_Disease2Gene_lit_lib)[i],]$diseaseLabel, " (", 
                        diseaseID_2_diseaseLabel[diseaseID_2_diseaseLabel$diseaseId == names(OpenTargets_Disease2Gene_lit_lib)[i],]$diseaseId, ")")
  names(OpenTargets_Disease2Gene_lit_lib)[i] <- diseaseName
}

# Filter to keep selected diseases
OpenTargets_Disease2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[lengths(OpenTargets_Disease2Gene_lit_lib) < 500]
OpenTargets_Disease2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[lengths(OpenTargets_Disease2Gene_lit_lib) > 20]

covid_genes <- OpenTargets_Disease2Gene_lit_lib[[grep("MONDO_0100096", names(OpenTargets_Disease2Gene_lit_lib))]]
write.csv(covid_genes, "Data/Exercise_SARS-Cov_genes.csv", row.names = FALSE)

res <- data.frame()
for(disease in names(OpenTargets_Disease2Gene_lit_lib)){
  disease_genes <- OpenTargets_Disease2Gene_lit_lib[[disease]]
  disease_net <- induced_subgraph(String_ppi_Net, 
                                  V(String_ppi_Net)[V(String_ppi_Net)$name %in% disease_genes])
  disease_net_clustcoef <- transitivity(disease_net, type = "global", isolates = "zero")
  
  covid_net <- induced_subgraph(String_ppi_Net, 
                                V(String_ppi_Net)[V(String_ppi_Net)$name %in% covid_genes])
  
  BBSI_sep <- Barabasi_proximity_separation(gene_network = String_ppi_Net, geneSet1 = covid_genes, geneSet2 = disease_genes)
  
  res <- rbind(res, data.frame("Disease" = disease, 
                               "Separation" = BBSI_sep,
                               "Disease_cc" = disease_net_clustcoef, 
                               "Disease_genes" = vcount(disease_net)))
  
}

write.csv(res, "Data_prep/Covid_vs_dis.csv", row.names = FALSE)

select_disease <- res[(res$Separation >= 0.5 &  # quantile(res$Separation, 0.75)
                         res$Disease_cc >= 0.5 &
                         res$Disease_genes >= 0), "Disease"]

select_disease <- c(select_disease, "COVID-19 (MONDO_0100096)")


select_disease 

select_disease_genes <- unique(unlist(OpenTargets_Disease2Gene_lit_lib[select_disease], recursive = FALSE, use.names = FALSE))

select_disease_net <- induced_subgraph(String_ppi_Net, 
                                       V(String_ppi_Net)[V(String_ppi_Net)$name %in% select_disease_genes])
Exercise_Net <- igraph::as_data_frame(select_disease_net, "edges")

write.csv(Exercise_Net, "Data/Exercise_PPI_Net.csv", row.names = FALSE)



# Create drugs
# Extract the human drug targets from DrugBank
DrugBank_Targets_idMap <- read.csv("Data_prep/targets_polypeptides_ext_id.csv", header = TRUE)
DrugBank_Targets_idMap <- reshape(as.data.frame(DrugBank_Targets_idMap), idvar = "parent_key", timevar = "resource", direction = "wide")
colnames(DrugBank_Targets_idMap) <- gsub("identifier.", "", colnames(DrugBank_Targets_idMap))
DrugBank_Targets_idMap <- DrugBank_Targets_idMap[grep("_HUMAN$", x = DrugBank_Targets_idMap$`UniProt Accession`),]

# Extract all drug target interactions in human
DrugBank_Drug_Target <- read.csv("Data_prep/targets.csv", header = TRUE)
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$organism == "Humans", ]
colnames(DrugBank_Drug_Target)[colnames(DrugBank_Drug_Target) == "parent_key"] <- "drug_id"
DrugBank_Drug_Target$drug_target <- DrugBank_Targets_idMap$GenAtlas[match(DrugBank_Drug_Target$id, DrugBank_Targets_idMap$parent_key)]

# Keep only approved drugs, experimental and investigational
DrugBank_drug_groups <- read.csv("Data_prep/drug_groups.csv", header = TRUE, check.names = FALSE)
DrugBank_drug_groups <- DrugBank_drug_groups[DrugBank_drug_groups$group %in% c("approved", "experimental", "investigational"), ]
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$drug_id %in% DrugBank_drug_groups$`drugbank-id`, ]

# Keep only drugs whose targets are in the network
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$drug_target %in% c(Exercise_Net$from, Exercise_Net$to),]

# Add the drug names
DrugBank_Drugs <- read.csv("Data_prep/drug.csv", header = TRUE)
DrugBank_Drug_Target$drug_name <- DrugBank_Drugs$name[match(DrugBank_Drug_Target$drug_id, DrugBank_Drugs$primary_key)]
DrugBank_Drugs <- DrugBank_Drugs[DrugBank_Drugs$type == "small molecule", ]
DrugBank_Drug_Target <- DrugBank_Drug_Target[DrugBank_Drug_Target$drug_id %in% DrugBank_Drugs$primary_key, ]

DrugBank_Drug_Target <- DrugBank_Drug_Target[, c("drug_id", "drug_name", "drug_target", "name")]
colnames(DrugBank_Drug_Target) <- c("drug_id", "drug_name", "drug_target", "drug_target_name")

# Remove drugs with certain categories
DrugBank_drug_categories <- read.csv("Data_prep/drug_categories.csv", header = TRUE)
DrugBank_drug_categories <- DrugBank_drug_categories[DrugBank_drug_categories$parent_key %in% DrugBank_Drug_Target$drug_id, ]

categories_remove <- c("Amino Acids", "Anions", "Cosmetics", "Food", "Ions", "Metals", "Metal cations", "Nucleotides")
remove_drugs <- c()

for(cat_select in categories_remove){
  tmp1 <- DrugBank_drug_categories[DrugBank_drug_categories$category == cat_select, ]
  remove_drugs <- c(remove_drugs, tmp1$parent_key)
}

DrugBank_Drug_Target <- DrugBank_Drug_Target[!DrugBank_Drug_Target$drug_id %in% remove_drugs, ]


write.csv(DrugBank_Drug_Target, "Data/Exercise_Drugs.csv", row.names = FALSE)


print(warnings())