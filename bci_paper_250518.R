
############################################################################################
############################################################################################
#Damla Cinoglu | March 24, 2025
############################################################################################
############################################################################################

library(dplyr)
library(tidyr)
library(tidyverse)
library(tibble)
library(openxlsx)
library(ggplot2)
library(lme4)
library(nlme)
library(lmerTest)
library(cowplot)
library(phylosignal)
library(phyloTop)      
library(MASS)
library(phylobase)
library(adephylo)
library(vegan)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(brms)
library(glmnet)
library(tidyr)library(dplyr)library(here)library(data.table)library(splitstackshape)library(rinchi)
library(ChemmineR)library(chemodiv)
library(ggtree)
library(tidybayes)
library(ggplot2)
library(dplyr)
library(ggridges)  
library(caper)
library(phytools)
library(geiger)
library(tidybayes)
library(arm)

setwd("/Users/damlacinoglu/Desktop/code and data")

load("DC_ALL_Metabmastertable_240920.RData")
pca_scores = read.csv("sp_pca_scores.csv", header = TRUE)
meta = read.xlsx("panama_seedlings_adults_metab_metadata_2024_v4.xlsx")
adult_trait = read.csv("BCITRAITS_20101220.csv")
sp_list = read.csv("species_list.csv", header = T)
Tree<-ape::read.tree("bci_phylogeny/BCI_Trees_Lianas_3gene.1.tre")
setwd("/Users/damlacinoglu/Desktop/final_panama_data")
all_individual_data = read.table("panama_data_analysis_ready_20241216.txt")

setwd("/Users/damlacinoglu/Desktop/code and data")
num_ind = read.csv("Copy of sp_list_for_brian_DC.csv")

############################################################################################
############################################################################################
#MORE DATA ORGANIZATION AND WRANGLING
############################################################################################
############################################################################################

#################################################
#Compile species level pools
#################################################

#new.master.real = master.real[,c(1:18)]
#new.master.real = cbind(new.master.real, master.real[,grepl("BCI", colnames(master.real))])

#meta $genu_species = paste0(meta $ATTRIBUTE_Genus,"_", meta $ATTRIBUTE_Species)
#meta $genu_species = gsub(' _',"_", meta $genu_species)
#meta $genu_species = gsub('_ ',"_", meta $genu_species)
#meta_seedling = meta[1:502,]

#colnames(master.real) = gsub(".Peak.area", "", colnames(master.real))
#colnames(master.real) = gsub("X", "", colnames(master.real))
#colnames(master.real) = gsub("mzML", "mzXML", colnames(master.real))
#for(i in 1:ncol(master.real)) {
#  if(colnames(master.real)[i] == "01.mzXML") {colnames(master.real)[i] =  "1.mzXML"}
#  if(colnames(master.real)[i] == "02.mzXML") {colnames(master.real)[i] =  "2.mzXML"}
#  if(colnames(master.real)[i] == "03.mzXML") {colnames(master.real)[i] =  "3.mzXML"}
#  if(colnames(master.real)[i] == "04.mzXML") {colnames(master.real)[i] =  "4.mzXML"}
#  if(colnames(master.real)[i] == "05.mzXML") {colnames(master.real)[i] =  "5.mzXML"}
#  if(colnames(master.real)[i] == "06.mzXML") {colnames(master.real)[i] =  "6.mzXML"}
#  if(colnames(master.real)[i] == "07.mzXML") {colnames(master.real)[i] =  "7.mzXML"}
#  if(colnames(master.real)[i] == "08.mzXML") {colnames(master.real)[i] =  "8.mzXML"}
#  if(colnames(master.real)[i] == "09.mzXML") {colnames(master.real)[i] =  "9.mzXML"}
#}

#for(i in unique(meta_seedling$genu_species)) {
  #files <- meta_seedling[meta_seedling[,"genu_species"] == i & meta_seedling[,"ATTRIBUTE_LEAF"] == #"YOUNG", "filename"]
#   files <- meta_seedling[meta_seedling[,"genu_species"] == i , "filename"]
 # new.master.real[[i]] <- rowMeans(master.real[, colnames(master.real) %in% files], na.rm = TRUE)
#}

#save(new.master.real, file = "DC_ALL_Metabmastertable_sp_250416_ALL.RData")

#################################################
#Individual level traits
#################################################

load("DC_ALL_Metabmastertable_240920.RData")

class = "secondary" 
if(class == "secondary"){
  heat.class = master.real[which(master.real$NPC.pathway %in% c("Alkaloids","Amino acids and Peptides","Polyketides","Shikimates and Phenylpropanoids", "Terpenoids")),]
}

nrow(heat.class)

colnames(heat.class) = gsub(".Peak.area", "", colnames(heat.class))
colnames(heat.class) = gsub("X", "", colnames(heat.class))
colnames(heat.class) = gsub("mzML", "mzXML", colnames(heat.class))
for(i in 1:ncol(heat.class)) {
  if(colnames(heat.class)[i] == "01.mzXML") {colnames(heat.class)[i] =  "1.mzXML"}
  if(colnames(heat.class)[i] == "02.mzXML") {colnames(heat.class)[i] =  "2.mzXML"}
  if(colnames(heat.class)[i] == "03.mzXML") {colnames(heat.class)[i] =  "3.mzXML"}
  if(colnames(heat.class)[i] == "04.mzXML") {colnames(heat.class)[i] =  "4.mzXML"}
  if(colnames(heat.class)[i] == "05.mzXML") {colnames(heat.class)[i] =  "5.mzXML"}
  if(colnames(heat.class)[i] == "06.mzXML") {colnames(heat.class)[i] =  "6.mzXML"}
  if(colnames(heat.class)[i] == "07.mzXML") {colnames(heat.class)[i] =  "7.mzXML"}
  if(colnames(heat.class)[i] == "08.mzXML") {colnames(heat.class)[i] =  "8.mzXML"}
  if(colnames(heat.class)[i] == "09.mzXML") {colnames(heat.class)[i] =  "9.mzXML"}
}

#meta = meta[,1:38]

unique(heat.class$NPC.pathway)

alk = heat.class[which(heat.class$NPC.pathway %in% c("Alkaloids")),]
shiki = heat.class[which(heat.class$NPC.pathway %in% c("Shikimates and Phenylpropanoids")),]
terp = heat.class[which(heat.class$NPC.pathway %in% c("Terpenoids")),]

for(i in 1:nrow(meta)) {
  vector = heat.class[,colnames(heat.class) == meta[i,"filename"]]
  meta[i,"richness"] = sum( vector != 0)
  
  vector1 = alk[,colnames(alk) == meta[i,"filename"]]
  meta[i,"Alkaloids_richness"] = sum( vector1 != 0)
  
  vector2 = shiki[,colnames(shiki) == meta[i,"filename"]]
  meta[i,"Shikimates_richness"] = sum( vector2 != 0)
  
  vector3 = terp[,colnames(terp) == meta[i,"filename"]]
  meta[i,"Terpenoids_richness"] = sum( vector3 != 0)
}

#To make up species level compound richness:
load("DC_ALL_Metabmastertable_sp_250412_YOUNGONLY.RData")

class = "secondary" 
if(class == "secondary"){
  heat.class = new.master.real[which(new.master.real $NPC.pathway %in% c("Alkaloids","Amino acids and Peptides","Polyketides","Shikimates and Phenylpropanoids", "Terpenoids")),]
}

nrow(heat.class)

alk = heat.class[which(heat.class$NPC.pathway %in% c("Alkaloids")),]
shiki = heat.class[which(heat.class$NPC.pathway %in% c("Shikimates and Phenylpropanoids")),]
terp = heat.class[which(heat.class$NPC.pathway %in% c("Terpenoids")),]

meta[,"genu_species"] = paste0(meta[,"ATTRIBUTE_Genus"],"_",meta[,"ATTRIBUTE_Species"] )
meta[,"genu_species"] = gsub(" ","", meta[,"genu_species"] )

for(i in 1:nrow(meta)) {
	ont = meta[i,"ATTRIBUTE_Ontogeny"]
	if(ont == "Seedling") {
		sp = meta[i,"genu_species"]
  vector = heat.class[[sp]]
  meta[i,"richness_sp"] = sum( vector != 0)
  
  vector1 = alk[,colnames(alk) == meta[i,"genu_species"]]
  meta[i,"Alkaloids_richness_sp"] = sum( vector1 != 0)
  
  vector2 = shiki[,colnames(shiki) == meta[i,"genu_species"]]
  meta[i,"Shikimates_richness_sp"] = sum( vector2 != 0)
  
  vector3 = terp[,colnames(terp) == meta[i,"genu_species"]]
  meta[i,"Terpenoids_richness_sp"] = sum( vector3 != 0)
	}
	}

tail(meta)

 meta[, "ATTRIBUTE_SPECIES"] = gsub(" ","", meta[, "ATTRIBUTE_SPECIES"] )
for(i in 1:nrow(meta)) {
		ont = meta[i,"ATTRIBUTE_Ontogeny"]
	if(ont == "Seedling") {
	sp = meta[i, "ATTRIBUTE_SPECIES"]
	meta[i,"num_ind"] = nrow(meta[meta[,"ATTRIBUTE_SPECIES"] == sp & meta[,"ATTRIBUTE_LEAF"] == "YOUNG", ])	
}
if(ont == "Adult") {
	sp = meta[i, "ATTRIBUTE_SPECIES"]
	if(sp != "BLANK") {
			meta[i,"num_ind"] = num_ind[num_ind[,1] == sp, "IndividualsInPool"] 
	}}
}

#meta[meta[,"num_ind"] == 0,]
	
meta[meta[,"ATTRIBUTE_Ontogeny"] == "Adult","richness_sp"] = meta[meta[,"ATTRIBUTE_Ontogeny"] == "Adult","richness"]
meta[meta[,"ATTRIBUTE_Ontogeny"] == "Adult","Alkaloids_richness_sp"] = meta[meta[,"ATTRIBUTE_Ontogeny"] == "Adult","Alkaloids_richness"]
meta[meta[,"ATTRIBUTE_Ontogeny"] == "Adult","Shikimates_richness_sp"] = meta[meta[,"ATTRIBUTE_Ontogeny"] == "Adult","Shikimates_richness"]
meta[meta[,"ATTRIBUTE_Ontogeny"] == "Adult","Terpenoids_richness_sp"] = meta[meta[,"ATTRIBUTE_Ontogeny"] == "Adult","Terpenoids_richness"]

#################################################
#Calculating chemo div and including those plots later
#################################################
#master.real has compounds as rows and samples as columns 
#load("DC_ALL_Metabmastertable_240920.RData")
#head(master.real)

load("DC_ALL_Metabmastertable_sp_250412_YOUNGONLY.RData")
master.real = new.master.real
head(master.real)
row.names(master.real) = master.real$row.ID
colnames(master.real)
master.real = master.real[,19:ncol(master.real)]
alpinaSampData = t(master.real)

#NOT THE ABSOLUTE VALUES, WE NEED RELATIVE ABUNDANCES
for(i in 1:nrow(alpinaSampData)) {
	alpinaSampData[i,] = alpinaSampData[i,] / sum(alpinaSampData[i,] )
}

#load("DC_ALL_Metabmastertable_sp_250412_YOUNGONLY.RData")
#master.real = new.master.real
#load("DC_ALL_Metabmastertable_240920.RData")
#head(master.real)

#initial = master.real[,c(1,4)]
#colnames(initial) = c("compound","smiles")
#head(initial)

#master.real[master.real[,"row.ID"] == 267,"smiles"]

#unique_structures = initial
#unique_structures$compound <- paste0("compound_",1:nrow(unique_structures))
#unique_structures <- unique_structures[unique_structures$smiles !=" ",]
#unique_structures <- unique_structures %>% filter(smiles != " ",smiles != "c")

#convert smiles to inchi-key
#unique_structures <- cbind(unique_structures,data.frame(inchikey=sapply(unique_structures$smiles, get.inchi.key),inchi=sapply(unique_structures$smiles, get.inchi)))
#use compDis to determine calculate the differences between smiles based on the PubChemFingerprints

# Filter the data frame to find rows where the 'smiles' column is not empty
#non_empty_smiles <-  initial[initial $smiles !=" ",]
#non_empty_smiles <- non_empty_smiles %>% filter(smiles != " ",smiles != "c")

# Extract the corresponding compound names into a vector
#compound_names <- non_empty_smiles$compound

#alpinaCompData = unique_structures
#alpinaCompData[,1] = compound_names

#alpinaCompDis <- compDis(compoundData = alpinaCompData,
 #                        type = "PubChemFingerprint")
#write.table(alpinaCompDis, file = "alpinaCompDis_MAY16")

#ONLY INCLUDE SECONDARY METABOLOME!
load("DC_ALL_Metabmastertable_sp_250412_YOUNGONLY.RData")
master.real = new.master.real
#load("DC_ALL_Metabmastertable_240920.RData")
#head(master.real)

secondary = unique(master.real[master.real[,"NPC.pathway"] %notin% c("Carbohydrates", "Fatty acids") & !is.na(master.real[,"NPC.pathway"]), "row.ID"])

secondary_compounds = intersect(secondary, compound_names)

# Assuming compound_names is a vector of column names to keep
alpinaSampData_subset <- alpinaSampData[, colnames(alpinaSampData) %in% secondary_compounds]

#alpinaCompDis$fingerDisMat <- alpinaCompDis$fingerDisMat[as.character(secondary_compounds), as.character(secondary_compounds)]

#colnames(alpinaCompDis$fingerDisMat) = secondary_compounds
#row.names(alpinaCompDis$fingerDisMat) = secondary_compounds

for(i in 1:nrow(alpinaSampData_subset)) {
	alpinaSampData_subset[i,] = alpinaSampData_subset[i,] / sum(alpinaSampData_subset[i,] )
}

alpinaDiv <- calcDiv(sampleData = alpinaSampData_subset, 
                     compDisMat = alpinaCompDis $fingerDisMat,
                     type = "FuncHillDiv",
                     q = 0)

alpinaDiv = cbind(alpinaDiv, row.names(alpinaSampData_subset))

alpinaDiv[,1] = alpinaDiv[,1]/2

write.table(alpinaDiv, file = "alpinaDiv_MAY16_SP_ALL")

#name = names(alpinaSampData_subset[1, alpinaSampData_subset[2, ] != 0])
#subset_fingerDisMat <- alpinaCompDis2$fingerDisMat[name, name]
#sum(subset_fingerDisMat[1:length(name),1:length(name)])
#sum_dissim <- sum(subset_fingerDisMat[upper.tri(subset_fingerDisMat)])
#sum(subset_fingerDisMat)

#data("alpinaSampData")
#data("alpinaCompData")
#alpinaCompDis5 <- compDis(compoundData = alpinaCompData,
#                         type = "PubChemFingerprint")
#alpinaDiv2 <- calcDiv(sampleData = alpinaSampData, 
 #                    compDisMat = alpinaCompDis5$fingerDisMat,
 #                    type = "FuncHillDiv",
 #                    q = 0)

#name = names(alpinaSampData[1, alpinaSampData[1, ] != 0])
#subset_fingerDisMat <- alpinaCompDis5 $fingerDisMat[name, name]
#sum(subset_fingerDisMat[1:length(name),1:length(name)])
#subset_fingerDisMat

##########################################################
##########################################################
##########################################################

#merge functional hill diversity here
hill_div = read.table(file = "alpinaDiv_MAY16_IND_ALL")
head(hill_div)
head(meta)
colnames(hill_div)[2] = "filename"

hill_div$filename = gsub(".Peak.area", "",hill_div$filename )
hill_div$filename  = gsub("X", "", hill_div$filename )
hill_div$filename  = gsub("mzML", "mzXML", hill_div$filename )
for(i in 1:nrow(hill_div)) {
  if(hill_div[i,"filename"] == "01.mzXML") {hill_div[i,"filename"] =  "1.mzXML"}
  if(hill_div[i,"filename"]  == "02.mzXML") {hill_div[i,"filename"]=  "2.mzXML"}
  if(hill_div[i,"filename"]  == "03.mzXML") {hill_div[i,"filename"] =  "3.mzXML"}
  if(hill_div[i,"filename"]  == "04.mzXML") {hill_div[i,"filename"] =  "4.mzXML"}
  if(hill_div[i,"filename"]  == "05.mzXML") {hill_div[i,"filename"] =  "5.mzXML"}
  if(hill_div[i,"filename"]  == "06.mzXML") {hill_div[i,"filename"] =  "6.mzXML"}
  if(hill_div[i,"filename"] == "07.mzXML") {hill_div[i,"filename"] =  "7.mzXML"}
  if(hill_div[i,"filename"] == "08.mzXML") {hill_div[i,"filename"] =  "8.mzXML"}
  if(hill_div[i,"filename"]== "09.mzXML") {hill_div[i,"filename"] =  "9.mzXML"}
}

meta = cbind(meta, hill_div = rep(NA, nrow(meta)))
for(i in 1:nrow(meta)) {
  if(length(hill_div[hill_div[,"filename"] == meta[i,"filename"], "FuncHillDiv"]) == 1) {
    meta[i, "hill_div"] = hill_div[hill_div[,"filename"] == meta[i,"filename"], "FuncHillDiv"]
  }
}
#NOTICE THERE ARE NO HILL DIV NUMBERS FOR THE BLANKS.

hill_div = read.table(file = "alpinaDiv_MAY16_SP_ALL")
head(hill_div)
head(meta)
colnames(hill_div)[2] = "filename"

hill_div$filename = gsub(".Peak.area", "",hill_div$filename )
hill_div$filename  = gsub("X", "", hill_div$filename )
hill_div$filename  = gsub("mzML", "mzXML", hill_div$filename )

hill_div[hill_div[,2] == "ylopia_macrantha",2]  = "Xylopia_macrantha"
meta = meta[meta[,"genu_species"] != "BLANK_BLANK",]
meta = meta[meta[,"genu_species"] != "NA_NA",]

meta = cbind(meta, hill_div_sp = rep(NA, nrow(meta)))
for(i in 1:nrow(meta)) {
genu = meta[i,"genu_species"]
ont = meta[i,"ATTRIBUTE_Ontogeny"]
if(ont == "Seedling") {
	if(genu %in%hill_div[,"filename"]) {
	    meta[i, "hill_div_sp"] = hill_div[hill_div[,"filename"] == genu, "FuncHillDiv"]
}
  }
  if(ont == "Adult") {
  if( meta[i,"filename"] %in% hill_div[,"filename"]){
	print("t")
	    meta[i, "hill_div_sp"] = hill_div[hill_div[,"filename"] == meta[i,"filename"], "FuncHillDiv"]	
}}
 
}

setwd("/Users/damlacinoglu/Desktop/final_panama_data")
all_individual_data = read.table("panama_data_analysis_ready_20241216.txt")
all_individual_data = all_individual_data[all_individual_data[,"Site"] %in% c("PEARSON", "BARBOUR"),]

all_individual_data$X..total.herbivory.per.individual <- as.numeric(all_individual_data$X..total.herbivory.per.individual)

summary_data <- all_individual_data %>%
 # filter(X..total.herbivory.per.individual > 0) %>%
  dplyr::group_by(SP) %>%
  dplyr::summarise(
    avg_herbivory = mean(X..total.herbivory.per.individual, na.rm = TRUE),
    se_herbivory = sd(X..total.herbivory.per.individual, na.rm = TRUE) / sqrt(n()),
    n = n()
  )
  
meta = cbind(meta, avg_herbivory = rep(NA, nrow(meta)), avg_sd = rep(NA, nrow(meta)))
for(i in 1:nrow(meta)) {
	ont = meta[i,"ATTRIBUTE_Ontogeny"]
	  sp =  gsub(" ", "", meta[i,"ATTRIBUTE_SPECIES"]) 
if(ont == "Seedling") {
  if(sp == "XYLIMA") {sp = "XYL1MA"}
  if(sp != "BLANK") {
    meta[i,"avg_herbivory"] = summary_data[summary_data[,"SP"] == sp,"avg_herbivory"]
    meta[i,"avg_sd"] = summary_data[summary_data[,"SP"] == sp,"se_herbivory"]
    }
  }
 if(ont == "Adult") { 
  if(nrow(pca_scores[pca_scores[,"sp"] == sp,]) == 1) {
    meta[i,"ATTRIBUTE_Group"] = pca_scores[pca_scores[,"sp"] == sp,"PFT_2axes"]
    meta[i,"ATTRIBUTE_PC2score"] = pca_scores[pca_scores[,"sp"] == sp,"PC2score"]
    meta[i,"ATTRIBUTE_PC1score"] = pca_scores[pca_scores[,"sp"] == sp,"PC1score"] 
  }
}
}
meta_seedling = meta[meta[,"ATTRIBUTE_Ontogeny"] == "Seedling",]
meta_adult = meta[meta[,"ATTRIBUTE_Ontogeny"] == "Adult",]

meta$new_ind = paste0(meta[,"ATTRIBUTE_SPECIES"],"_",meta[,"ATTRIBUTE_INDIVIDUAL"])

meta[meta[,"ATTRIBUTE_SPECIES"] == "COCCMA", "ATTRIBUTE_Group"] = 1

for(i in 1:nrow(meta)) {
	sp = meta[i,"ATTRIBUTE_SPECIES"] 
	if(!is.na(meta[i,"ATTRIBUTE_Group"])) {
		if(as.numeric(meta[i,"ATTRIBUTE_Group"]) != pca_scores[pca_scores[,"sp"] == sp,"PFT_2axes"]) {
		print(i)
	}		
	}
	}
	
#Check if hill diversity is possible?
store = NULL
for(i in 1:nrow(meta)) {
	if(meta[i,"richness_sp"] > 0) { 	if((((meta[i,"richness_sp"]^2) - meta[i,"richness_sp"])/2) < meta[i,"hill_div_sp"]) {print(i)}
} else { store = append(store, i)}
}
	
############################################################################################
############################################################################################
#FIGURES AND ANALYSES
############################################################################################
############################################################################################

#################################################
#Number of individuals figure
#################################################
 species_new <- c(
  "BROSAL", "ACALDI", "ANNOSP", "CROTBI", "COU2CU", "CASECO", "CAPPFR", "FARAOC", "GUSTSU", "EUGEGA",
  "AMAICO", "GUATDU", "HEISCO", "EUGENE", "EUGECO", "SOROAF", "PSYCDE", "COCCMA", "ALSEBL", "TET2PA",
  "HYBAPR", "EUGEOE", "MICOAR", "SIMAAM", "PSYCLI", "INGACO", "POUTST", "NECTPU", "POCHSE", "INGAPE",
  "PALIGU", "OURALU", "PSYCHO", "PROTTE", "XYL1MA", "MOURMY", "PSYCMA", "PROTPA", "LUEHSE", "PIPECO",
  "LACIAG", "LAETTH"
) 

  species_vector <- c(
  "Acalypha_diversifolia", "Alseis_blackiana", "Annona_spraguei",
  "Coussarea_curvigemmia","Coccoloba_manzinellensis",
  "Croton_billbergianus", "Eugenia_galalonensis", "Eugenia_nesiotica",
  "Eugenia_oerstediana", "Faramea_occidentalis", "Guatteria_dumetorum",
  "Gustavia_superba", "Heisteria_concinna", "Hybanthus_prunifolius",
  "Inga_cocleensis", "Lacistema_aggregatum", "Laetia_thamnia",
  "Luehea_seemannii", "Miconia_argentea", "Mouriri_myrtilloides",
  "Ouratea_lucens", "Palicourea_guianensis", "Piper_cordulatum",
  "Protium_panamense", "Protium_tenuifolium", "Psychotria_deflexa",
  "Psychotria_horizontalis", "Psychotria_limonensis", "Psychotria_marginata",
  "Simarouba_amara", "Sorocea_affinis", "Tetragastris_panamensis",
  "Xylopia_macrantha"
)
length(species_vector)
nrow(meta_unique[meta_unique[,"ATTRIBUTE_Ontogeny"] == "Seedling",])
nrow(meta_unique[meta_unique[,"ATTRIBUTE_Ontogeny"] == "Adult",])

# "Brosimum_alicastrum",

meta_newish = meta[meta[,"genu_species"] %in% species_vector,]
meta_unique <- meta_newish %>%
  distinct(ATTRIBUTE_SPECIES, ATTRIBUTE_Ontogeny, genu_species, .keep_all = TRUE)
  
model = lmer(richness_sp ~ num_ind +(1| ATTRIBUTE_Ontogeny), data = meta_unique)
summary(model)
model = lm(richness_sp ~ num_ind , data = meta_unique[meta_unique[,"ATTRIBUTE_Ontogeny"] == "Seedling",])

summary(model)
fig <- ggplot() +
  theme_bw() +
  theme(axis.ticks.x = element_blank()) +
  labs(y = "Secondary metabolomic compound richness\n(Number of unique compounds)", 
       x = "Number of individuals in species pool") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  geom_point(data =  meta_unique, aes(x = as.numeric(num_ind), y = as.numeric(richness_sp)), color = "gray") +
  geom_smooth(data =  meta_unique, aes(x = as.numeric(num_ind), y = as.numeric(richness_sp)), method = "lm", color = "blue", size = 1) +
  theme(text = element_text(size = 11))
  
model = lmer(hill_div_sp ~ num_ind +(1| ATTRIBUTE_Ontogeny), data = meta_unique)
summary(model)
model = lm(hill_div_sp ~ num_ind , data = meta_unique[meta_unique[,"ATTRIBUTE_Ontogeny"] == "Seedling",])

fig2 <- ggplot() +
  theme_bw() +
  theme(axis.ticks.x = element_blank()) +
  labs(y = "Functional hill diversity\n(Sum of the pairwise dissimilarities\n in the compound dissimilarity matrix", 
       x = "Number of individuals in species pool") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) + 
  geom_point(data =  meta_unique, aes(x = as.numeric(num_ind), y = as.numeric(hill_div_sp)), color = "gray") +
  geom_smooth(data =  meta_unique, aes(x = as.numeric(num_ind), y = as.numeric(hill_div_sp)), method = "lm", color = "blue", size = 1) +
  theme(text = element_text(size = 11))

final_plot <- plot_grid(fig, fig2, labels = c('A', 'B'), label_size = 12, ncol = 2)

ggsave("/Users/damlacinoglu/Desktop/code and data/figures/num_ind_fig_plot.png", plot = final_plot, width = 7, height = 4, dpi = 300)

#################################################
#Table comparing ontogeny, life history, phylogeny
#################################################

meta_newish $ATTRIBUTE_PC1score = as.numeric(meta_newish $ATTRIBUTE_PC1score)
meta_newish $ATTRIBUTE_PC2score = as.numeric(meta_newish $ATTRIBUTE_PC2score)
meta_newish$ATTRIBUTE_Group = as.numeric(meta_newish$ATTRIBUTE_Group)

summary_data <- meta_newish %>%
  dplyr::group_by(genu_species, ATTRIBUTE_Ontogeny) %>%
  dplyr::summarize(
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  as.data.frame()

# Add 'name' column: Extract everything before the second underscore
summary_data$name <- sapply(strsplit(summary_data$genu_species, "_"), function(x) paste(x[1], x[2], sep = "_"))

phylo_tree <-ape::read.tree("bci_phylogeny/BCI_Trees_Lianas_3gene.1.tre")

summary_data$name <- gsub(" ", "", summary_data$name)
species_in_data <- unique(summary_data$name)
trimmed_phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, species_in_data))

summary_data_final = summary_data[summary_data[,"name"] %in% trimmed_phylo_tree$tip,]
unique(summary_data_final$name)

library(nlme)
summary_data_final_clean <- summary_data_final[complete.cases(summary_data_final$richness_sp,
                                                              summary_data_final$ATTRIBUTE_PC1score,
                                                              summary_data_final$ATTRIBUTE_PC2score,
                                                              summary_data_final$ATTRIBUTE_Group,
                                                              summary_data_final$ATTRIBUTE_Ontogeny), ]
       
       summary_data_final_clean = data.frame(summary_data_final_clean)                                                       
 summary_data_final_clean[,3:21] = as.numeric( summary_data_final_clean[,3:21] )

species_in_data <- unique(summary_data_final$name)
trimmed_phylo_tree <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, species_in_data))

A <- vcv.phylo(phylo_tree)
# Extract species names from your data
species_in_data <- unique(summary_data_final_clean$name)

# Match the species in your data to those in the phylogenetic tree
phylo_tree_species <- phylo_tree$tip.label
species_in_tree <- intersect(species_in_data, phylo_tree_species)

# Subset the phylogenetic covariance matrix to match the species in the data
A_trimmed <- A[species_in_tree, species_in_tree]
A_normalized <- A_trimmed / sqrt(diag(A_trimmed) %*% t(diag(A_trimmed)))
diag(A_normalized) = diag(A_normalized)+ 0.0000001

hist(log(summary_data_final_clean$richness_sp))
summary_data_final_clean[,"name"] = paste0(summary_data_final_clean[,"genu_species"],"_",summary_data_final_clean[,"ATTRIBUTE_Ontogeny"]  )

model_repeat1 <- brm(
  as.numeric(log(richness_sp))~ as.factor(ATTRIBUTE_Group) + ATTRIBUTE_Ontogeny + 
    (1 | gr(genu_species, cov = A)) +(1 | name) +  (1 | num_ind),
  data = summary_data_final_clean,
  family = gaussian(),
  data2 = list(A = A_normalized),  # Pass the normalized A matrix in data2
  prior = c(
    prior(normal(0, 10), "b"),
    prior(normal(0, 10), "Intercept"),
    prior(student_t(3, 0, 20), "sd"),
    prior(student_t(3, 0, 20), "sigma")
  ),
  sample_prior = TRUE, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

#There's a 95% probability that seedlings have higher richness than the reference level.
posterior <- as_draws_df(model_repeat1)
mean(posterior$ b_ATTRIBUTE_OntogenySeedling  > 0)
posterior <- model_repeat1 %>%
  spread_draws(
    b_ATTRIBUTE_OntogenySeedling,
    sd_genu_species__Intercept,
    b_as.factorATTRIBUTE_Group2,
    b_as.factorATTRIBUTE_Group3,
    b_as.factorATTRIBUTE_Group4
  ) %>%
  mutate(
    b_ATTRIBUTE_Group_Mean = rowMeans(
      cbind(0, as.matrix(across(starts_with("b_as.factorATTRIBUTE_Group"))))
    )
  )

posterior_long <- posterior %>%
 dplyr:: select(
    b_ATTRIBUTE_OntogenySeedling,
    sd_genu_species__Intercept,
    b_ATTRIBUTE_Group_Mean
  ) %>%
 pivot_longer(cols = everything(), names_to = "term", values_to = "estimate") %>%
  dplyr:: mutate(term = case_when(
    term == "b_ATTRIBUTE_OntogenySeedling" ~ "Ontogeny",
    term == "sd_genu_species__Intercept" ~ "Phylogeny",
    term == "b_ATTRIBUTE_Group_Mean" ~ "Demographic Group"
  ))
  
summary_ontogeny <- posterior_long %>%
  filter(term == "Ontogeny") %>%
  summarise(
    posterior_mean = mean(estimate),
    lower_95_CrI = quantile(estimate, 0.025),
    upper_95_CrI = quantile(estimate, 0.975)
  )

summary_phylogeny <- posterior_long %>%
  filter(term == "Phylogeny") %>%
  summarise(
    posterior_mean = mean(estimate),
    lower_95_CrI = quantile(estimate, 0.025),
    upper_95_CrI = quantile(estimate, 0.975)
  )

summary_demographic <- posterior_long %>%
  filter(term == "Demographic Group") %>%
  summarise(
    posterior_mean = mean(estimate),
    lower_95_CrI = quantile(estimate, 0.025),
    upper_95_CrI = quantile(estimate, 0.975)
  )

plot_1 =  ggplot(
  posterior_long,
  aes(x = estimate, y = term, fill = term)
) +
  stat_halfeye(.width = c(0.66, 0.95)) +  # Removed the point_interval argument
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Posterior estimate",
    y = NULL,
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")

#SPECIES LEVEL HILL DIVERSITY

model_repeat1 <- brm(
  as.numeric(log(hill_div_sp))~ as.factor(ATTRIBUTE_Group) + ATTRIBUTE_Ontogeny + 
    (1 | gr(genu_species, cov = A)) +(1 | name) +  (1 | num_ind),
  data = summary_data_final_clean,
  family = gaussian(),
  data2 = list(A = A_normalized),  # Pass the normalized A matrix in data2
  prior = c(
    prior(normal(0, 10), "b"),
    prior(normal(0, 10), "Intercept"),
    prior(student_t(3, 0, 20), "sd"),
    prior(student_t(3, 0, 20), "sigma")
  ),
  sample_prior = TRUE, chains = 2, cores = 2,
  iter = 4000, warmup = 1000
)

#There's a 95% probability that seedlings have higher richness than the reference level.
posterior <- as_draws_df(model_repeat1)
mean(posterior$ b_ATTRIBUTE_OntogenySeedling  > 0)
posterior <- model_repeat1 %>%
  spread_draws(
    b_ATTRIBUTE_OntogenySeedling,
    sd_genu_species__Intercept,
    b_as.factorATTRIBUTE_Group2,
    b_as.factorATTRIBUTE_Group3,
    b_as.factorATTRIBUTE_Group4
  ) %>%
  mutate(
    b_ATTRIBUTE_Group_Mean = rowMeans(
      cbind(0, as.matrix(across(starts_with("b_as.factorATTRIBUTE_Group"))))
    )
  )

posterior_long <- posterior %>%
 dplyr:: select(
    b_ATTRIBUTE_OntogenySeedling,
    sd_genu_species__Intercept,
    b_ATTRIBUTE_Group_Mean
  ) %>%
 pivot_longer(cols = everything(), names_to = "term", values_to = "estimate") %>%
  dplyr:: mutate(term = case_when(
    term == "b_ATTRIBUTE_OntogenySeedling" ~ "Ontogeny",
    term == "sd_genu_species__Intercept" ~ "Phylogeny",
    term == "b_ATTRIBUTE_Group_Mean" ~ "Demographic Group"
  ))

summary_ontogeny <- posterior_long %>%
  filter(term == "Ontogeny") %>%
  summarise(
    posterior_mean = mean(estimate),
    lower_95_CrI = quantile(estimate, 0.025),
    upper_95_CrI = quantile(estimate, 0.975)
  )

summary_phylogeny <- posterior_long %>%
  filter(term == "Phylogeny") %>%
  summarise(
    posterior_mean = mean(estimate),
    lower_95_CrI = quantile(estimate, 0.025),
    upper_95_CrI = quantile(estimate, 0.975)
  )

summary_demographic <- posterior_long %>%
  filter(term == "Demographic Group") %>%
  summarise(
    posterior_mean = mean(estimate),
    lower_95_CrI = quantile(estimate, 0.025),
    upper_95_CrI = quantile(estimate, 0.975)
  )
  
plot_2 = ggplot(
  posterior_long,
  aes(x = estimate, y = term, fill = term)
) +
  stat_halfeye(.width = c(0.66, 0.95)) +  # Removed the point_interval argument
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Posterior estimate",
    y = NULL,
    title = ""
  ) +
  theme_bw() +
  theme(legend.position = "none")

final_plot <- plot_grid(plot_1, plot_2,  labels = c('A', 'B'), ncol = 2, label_size = 12)
ggsave("/Users/damlacinoglu/Desktop/code and data/figures/posterior_plot.png", final_plot, width = 8, height = 4)

#############################################
#COMPARING SEEDLING AND ADULT METABOLOMIC SPACE
#############################################
sp_cscs = read.table("cscs_APR12_young.txt")
cscs = sp_cscs
head(cscs)

nmds.samps = isoMDS(as.dist(1-cscs), k=2)
nmds.samps.points = nmds.samps$points

meta_newish$ATTRIBUTE_Group = as.numeric(meta_newish $ATTRIBUTE_Group)  # Ensure it's numeric
valid_species = meta_newish $genu_species[meta_newish $ATTRIBUTE_Group %in% c(1, 2, 3, 4)] 

seedlings = nmds.samps.points[row.names(nmds.samps.points) %in% species_vector, ]
#seedlings = nmds.samps.points[row.names(nmds.samps.points)%in% unique(meta[meta["ATTRIBUTE_Ontogeny"] == "Seedling", "genu_species"]), ]
adults = meta_newish[meta_newish[,"ATTRIBUTE_Ontogeny"] == "Adult",]
files = adults [adults[,"genu_species"] %in% row.names(seedlings), "filename"]
row.names(nmds.samps.points) = gsub(".Peak.area","",row.names(nmds.samps.points))
adults_nmds = nmds.samps.points[row.names(nmds.samps.points) %in% files, ]

final_nmds = rbind(seedlings, adults_nmds)

for(i in 1:nrow(final_nmds)) {
	 sp = row.names(final_nmds)[i]
if(grepl( "BCI", sp)) {
	ColorCode[i] = colors[2]
}	
}

ColorCode <- ifelse(grepl("BCI", rownames(final_nmds)), "black", "red")

final_nmds <- as.data.frame(final_nmds)
final_nmds = cbind(final_nmds, ColorCode)
final_nmds $Stage <- ifelse(final_nmds $ColorCode == "black", "Adult", "Seedling")

space_plot = ggplot(final_nmds, aes(x = V1, y = V2, color = Stage, fill = Stage)) +
  stat_ellipse(geom = "polygon", alpha = 0.2) +  
  geom_point(size = 3, shape = 16) +  
  scale_color_manual(values = c("Seedling" = "red", "Adult" = "black")) + 
    scale_fill_manual(values = c("Seedling" = "red", "Adult" = "black")) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Ontogeny", fill = "Ontogeny") +
  theme_bw() +
  theme(legend.position = "right") +   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

species_vector %in% unique(meta[meta[,"ATTRIBUTE_SPECIES"] %in% species_new,"genu_species"])  
# Wilcoxon test for NMDS Axis 1 and 2
wilcox.test(final_nmds[,1] ~ final_nmds[,"Stage"])
wilcox.test(final_nmds[,2] ~ final_nmds[,"Stage"])

for(i in 1:nrow(final_nmds)) {
	if(final_nmds[i,"Stage"] == "Seedling") {
		final_nmds[i,"num_ind"] = meta[meta[,"genu_species"] == rownames(final_nmds)[i] & 
	meta[,"ATTRIBUTE_Ontogeny"] == final_nmds[i,"Stage"], "num_ind"] [1]
	}
	if(final_nmds[i,"Stage"] == "Adult") {
		final_nmds[i,"num_ind"] = meta[meta[,"filename"] == rownames(final_nmds)[i] & 
	meta[,"ATTRIBUTE_Ontogeny"] == final_nmds[i,"Stage"], "num_ind"] [1]
	}		
}
model1 <- lmer(final_nmds[,1] ~ Stage + (1 | num_ind), data = final_nmds)
summary(model1)
model1 <- lmer(final_nmds[,2] ~ Stage + (1 | num_ind), data = final_nmds)
summary(model1)

simulation_output <- simulateResiduals(fittedModel = model1, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)

#Start of the euclidiean space analyses
final_nmds = cbind(final_nmds, Group = rep(NA, nrow(final_nmds)))
for(i in 1:nrow(final_nmds)) {
  sp = row.names(final_nmds)[i]
     if(final_nmds[i,"Stage"] == "Seedling") {
  	 group = as.numeric(meta_newish[meta_newish[,"genu_species"] == sp, "ATTRIBUTE_Group"][1])
 if(!is.na(group)) {
  	  final_nmds[i, "Group"] = group  }
  }
  
  	    if(final_nmds[i,"Stage"] == "Adult") {
  group = as.numeric(meta_newish[meta_newish[,"filename"] == sp, "ATTRIBUTE_Group"][1])
 if(!is.na(group)) {
  	  final_nmds[i, "Group"] = group  }
  }
}

filtered_data <- final_nmds[final_nmds[, "Group"] %in% c(1, 2, 3, 4), ]

# Step 2: Extract seedling and adult data
seedlings <- filtered_data[filtered_data[,"Stage"] == "Seedling", ]
adults <-  filtered_data[filtered_data[,"Stage"] == "Adult", ]

# Step 3: Match species between seedlings and adults by their names before "_Seedling" and "_Adult"
# Extract the species names by removing the "_Seedling" or "_Adult" suffix
#species_names <- gsub("_Seedling|_Adult", "", rownames(seedlings))

# Ensure we have matching species between Seedling and Adult
# Step 1: Extract the species names from seedlings and adults by removing the "_Seedling" or "_Adult" suffix
species_names_seedlings <- rownames(seedlings)
species_names_adults <- meta_newish[meta_newish [,"filename"] %in% rownames(adults),"genu_species"]

# Step 2: Match species names between seedlings and adults
matched_seedlings <- seedlings[species_names_seedlings %in% species_names_adults, ]
matched_seedlings_sorted <- matched_seedlings[order(rownames(matched_seedlings)), ]

matched_adults <- adults[species_names_adults %in% species_names_seedlings, ]
matched_adults_sorted <- matched_adults[order(rownames(matched_adults)), ]

# Step 4: Calculate Euclidean distance for each species between Seedling and Adult
distances <- apply(cbind(matched_seedlings_sorted[, 1:2], matched_adults_sorted[, 1:2]), 1, function(x) {
  sqrt((x[1] - x[3])^2 + (x[2] - x[4])^2)  # Euclidean distance calculation
})

# Add distances to the dataset
matched_data <- data.frame(Species = gsub("_Adult","",row.names(matched_adults_sorted)), Distance = distances, Group = matched_adults_sorted[, "Group"], PC1_score = NA, PC2_score = NA)

for(i in 1:nrow(matched_data)) {
	matched_data[i, "PC1_score"] = as.numeric(meta_newish[meta_newish[,"genu_species"] == row.names(matched_data)[i], "ATTRIBUTE_PC1score"    ][1])
	matched_data[i, "PC2_score"] = as.numeric(meta_newish[meta_newish[,"genu_species"] == row.names(matched_data)[i], "ATTRIBUTE_PC2score"    ][1])
}

matched_data = cbind(matched_data, num_ind = rep(NA, nrow(matched_data)) )
for(i in 1:nrow(matched_data)) {
	matched_data[i, "num_ind"] = as.numeric(meta_newish[meta_newish[,"genu_species"] == row.names(matched_data)[i], "num_ind"][1]) + as.numeric(meta_newish[meta_newish[,"filename"] == matched_data[i,"Species"], "num_ind"][1])
}

#lmm_model <- lmer(Distance ~ scale(PC1_score) *scale(PC2_score) + (1|num_ind), data = matched_data)
lmm_model <- lmer(Distance ~ as.factor(Group) + (1|num_ind), data = matched_data)
summary(lmm_model)

emm1.1 <- emmeans(lmm_model, specs = pairwise ~ as.factor(Group), type = "response", adjust = "tukey")
 multcomp::cld(emm1.1$emmeans, Letters = letters, adjust = "tukey")
 
# Step 2: Extract the residuals from the model
residuals_model2 <- residuals(lmm_model)

# Step 3: Assess phylogenetic signal in residuals
# Read the phylogenetic tree
library(ape)
Tree <- read.tree("bci_phylogeny/BCI_Trees_Lianas_3gene.1.tre")

#residuals_data <- data.frame(species = meta_adult $genu_species, residuals = residuals_model2)

Tree_filtered <- drop.tip(Tree, setdiff(Tree$tip.label, names(residuals_model2)))
residuals_data2 = residuals_model2[names(residuals_model2) %in% Tree_filtered$tip.label]
names(residuals_data2)

#names(residuals_model2) = residuals_data2[,"species"]
lambda <- phylosig(Tree_filtered, residuals_model2,method="K", test=TRUE, nsim=999)
print(lambda)

pairwise_results <- emmeans(lmm_model, pairwise ~ Group, adjust = "none")
summary(pairwise_results)

new_meta = meta[meta[,"genu_species"] %in% rownames(matched_data),]
new_meta_unique <- new_meta %>%
  dplyr::distinct(genu_species, .keep_all = TRUE)
#new_meta_unique = new_meta_unique[,-c(ncol(new_meta_unique))]

new_meta_adult = new_meta[new_meta[,"ATTRIBUTE_Ontogeny"] == "Adult",]
#new_meta_adult = new_meta_adult[,-c(ncol(new_meta_adult))]

all = rbind(new_meta_adult, new_meta_unique[new_meta_unique[,"ATTRIBUTE_Ontogeny"]== "Seedling",])

model = lmer(richness_sp ~ ATTRIBUTE_Ontogeny + (1|genu_species)  + (1|num_ind), data = all)
model = lmer(log(richness_sp) ~ ATTRIBUTE_Ontogeny + (1|genu_species)  + (1|num_ind), data = all)

summary(model)
simulation_output <- simulateResiduals(fittedModel = model1, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)

distance_plot = ggplot(matched_data, aes(x = factor(Group), y = Distance)) +
  geom_jitter(color = "gray", size = 2, width = 0.2, show.legend = FALSE) +
  stat_summary(
    fun = mean,
    fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)),
    fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
    geom = "pointrange",
    aes(group = Group, color = factor(Group)),
    position = position_dodge(width = 0.5),
    size = 0.8
  ) +
  scale_color_manual(values = c("1" = "black", "2" = "black", "3" = "black", "4" = "black")) +
  scale_x_discrete(labels = c("1" = "Slow", "2" = "Fast", "3" = "Long-lived\n pioneer", "4" = "Short-lived\n breeder")) +
  theme_bw() +
  labs(
    x = "",
    y = "Euclidean distance in metabolomic space"
  ) +
  theme(legend.position = "none")+   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

richness_plot <- ggplot(all, aes(x = ATTRIBUTE_Ontogeny, y = richness_sp)) +
  geom_jitter(aes(color = ATTRIBUTE_Ontogeny), size = 2, show.legend = FALSE, color ="gray") +
    stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)),
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = "pointrange",
               aes(group = ATTRIBUTE_Ontogeny,
                   color = ATTRIBUTE_Ontogeny),
               position = position_dodge(width = 0.8), size = 0.8) +
    scale_color_manual(values = c("Seedling" = "red", "Adult" = "black")) +
  theme_bw() +
  labs(
    x = "",
    y = "Secondary metabolomic compound richness\n(Number of unique compounds)"
  ) +
  theme(legend.position = "none")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

model = lmer(hill_div_sp ~ ATTRIBUTE_Ontogeny + (1|genu_species)  + (1|num_ind), data = all)

simulation_output <- simulateResiduals(fittedModel = model, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)

summary(model)

simulation_output <- simulateResiduals(fittedModel = model1, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)

hill_plot <- ggplot(all, aes(x = ATTRIBUTE_Ontogeny, y = hill_div_sp)) +
  geom_jitter(aes(color = ATTRIBUTE_Ontogeny), size = 2, show.legend = FALSE, color ="gray") +
    stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x)/sqrt(length(x)),
               fun.max = function(x) mean(x) + sd(x)/sqrt(length(x)),
               geom = "pointrange",
               aes(group = ATTRIBUTE_Ontogeny,
                   color = ATTRIBUTE_Ontogeny),
               position = position_dodge(width = 0.8), size = 0.8) +
    scale_color_manual(values = c("Seedling" = "red", "Adult" = "black")) +
  theme_bw() +
  labs(
    x = "",
    y = "Functional hill diversity\n(Sum of the pairwise dissimilarities\n in the compound dissimilarity matrix)"
  ) +
  theme(legend.position = "none") +   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    
final_plot <- plot_grid(space_plot, distance_plot, richness_plot, hill_plot, labels = c('A', 'B', "C", "D"), ncol = 2, label_size = 12)

ggsave("/Users/damlacinoglu/Desktop/code and data/figures/fig_2.png", final_plot, width = 7, height = 8)

#################################################
#Herbivory distribution bar plot
#################################################

group_colors <- c(`1` = "pink", `2` = "orange", `3` = "green", `4` = "blue", `5` = "red")
group_labels <- c(`1` = "Slow", `2` = "Fast", 
                  `3` = "Long-lived pioneers", `4` = "Short-lived breeders", `5` = "Intermediate")

summary_data <- meta[meta["ATTRIBUTE_Ontogeny"] == "Seedling", ] %>%
  dplyr::group_by(ATTRIBUTE_SPECIES) %>%
  dplyr::summarise(
    avg_herbivory = mean(avg_herbivory, na.rm = TRUE),
    avg_sd = mean(avg_sd, na.rm = TRUE),
    ATTRIBUTE_Group = dplyr::first(ATTRIBUTE_Group),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(avg_herbivory))

summary_data = data.frame(summary_data)
summary_data[,"ATTRIBUTE_Group"] = as.factor(summary_data[,"ATTRIBUTE_Group"] )
summary_data = summary_data[summary_data[,"ATTRIBUTE_Group"] %in% c('1','2','3','4','5'),]

summary_data$ATTRIBUTE_SPECIES <- factor(summary_data$ATTRIBUTE_SPECIES, 
                                         levels = summary_data$ATTRIBUTE_SPECIES[order(summary_data$avg_herbivory, decreasing = TRUE)])

summary_data[,1] %in% species_new

# Then plot
plot <- ggplot(summary_data, aes(x = ATTRIBUTE_SPECIES, y = avg_herbivory, color = as.character(ATTRIBUTE_Group))) +
  geom_point(size = 3) +  # Dots for mean herbivory
  geom_errorbar(aes(ymin = avg_herbivory - avg_sd, ymax = avg_herbivory + avg_sd), width = 0.2) +  
  scale_color_manual(values = group_colors, labels = group_labels, name = "Demographic group") +
  labs(x = "Species", y = "Total herbivory per individual (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("/Users/damlacinoglu/Desktop/code and data/figures/herbivory_plot.png", plot = plot, width = 8, height = 4, dpi = 300)

#################################################
#PCA figure with the species
#################################################
#Plot all PC scores of all species for the background.
spdata = pca_scores[,c(1,4,5)]; spdata[,1] = tolower(spdata[,1]); colnames(spdata) = c("sp","PC1","PC2")
spdata$dg = NA

sp_cords_seedlings = unique(meta[meta[,"ATTRIBUTE_Ontogeny"] == "Seedling" ,"ATTRIBUTE_SPECIES"])

new_coordinates = data.frame(sp = sp_cords_seedlings, dg =rep(NA,length(sp_cords_seedlings)), PC1 = rep(NA,length(sp_cords_seedlings)), PC2 = rep(NA,length(sp_cords_seedlings)))

for(i in 1:nrow(new_coordinates)) {
	if(new_coordinates[i,"sp"] %in% pca_scores[,"sp"]) {
		new_coordinates[i,"dg"] =  pca_scores[pca_scores[,"sp"] == new_coordinates[i,"sp"],"PFT_2axes"]
  new_coordinates[i,"PC1"] = pca_scores[pca_scores[,"sp"] == new_coordinates[i,"sp"],"PC1score"]
  new_coordinates[i,"PC2"] = pca_scores[pca_scores[,"sp"] == new_coordinates[i,"sp"],"PC2score"]		
	} else {new_coordinates[i,"dg"] =  6
  	}
  }

new_coordinates$dg <- factor(new_coordinates$dg, levels = c(1, 2, 3, 4, 5, 6), 
                             labels = c("Slow", "Fast", "Long-lived pioneer", "Short-lived breeder", "Intermediate", "Not assigned"))
#new_coordinates = new_coordinates[!is.na(new_coordinates[,"dg"]),]

new_coordinates = new_coordinates[new_coordinates[,1] %in% species_new,]

pca_fig = ggplot(spdata, aes(x = PC1, y = PC2)) + 
  theme_bw() + 
  theme(axis.ticks.x = element_blank()) +
  labs(y = "Stature-recruitment tradeoff", x = "Growth-survival tradeoff") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_point(color = "gray") +
  theme(text = element_text(size = 11)) + 
  geom_vline(xintercept = 0, color = "black") + 
  geom_hline(yintercept = 0, color = "black") +
  geom_point(data = new_coordinates, aes(x = PC1, y = PC2, color = dg), size = 2) +  
    theme(legend.position = "none") + 
   scale_color_manual(values = c("Slow" = "pink", 
                                "Fast" = "orange", 
                               "Long-lived pioneer" = "green", 
                                "Short-lived breeder" = "blue", 
                                 "Intermediate" = "red")) +
                               labs(color="Demographic group") +
  theme(legend.position = c(0.20, 0.17)) 

ggsave("/Users/damlacinoglu/Desktop/code and data/figures/fig_1A.png", pca_fig, width = 5, height = 6, dpi = 300)

################################
#phylogenetic analysis
################################

meta_phylo = meta[!meta[,"genu_species"] == "Miconia_NA",]
meta_phylo = meta_phylo[meta_phylo[,"ATTRIBUTE_Ontogeny"] == "Seedling",]

# Define the tip labels you want to keep
unique(meta_phylo$genu_species)

meta_phylo = meta_phylo[meta_phylo$ATTRIBUTE_SPECIES %in% new_coordinates[,1],]
selected_labels = unique(meta_phylo $genu_species[meta_phylo $genu_species%in%Tree$tip.label])
# Remove specific items from selected_labels
#selected_labels <- selected_labels[!selected_labels %in% c("Amaioua_corymbosa", "Eugenia_coloradoensis")]

# Subset the tree to keep only the selected tip labels
sub_tree <- ape::drop.tip(Tree, tip = Tree$tip.label[!(Tree$tip.label %in% selected_labels),drop = T])

# Create a named vector of colors
color_map <- setNames(c("pink", "orange", "green", "blue", "red", "gray"), c(1, 2, 3, 4, 5, 6))

# Match the `ATTRIBUTE_Group` directly to `sub_tree$tip.label`
tip_groups <- meta_phylo $ATTRIBUTE_Group[match(sub_tree$tip.label, meta_phylo $genu_species)]

# Assign colors, default to grey if missing
colors <- color_map[as.character(tip_groups)]
colors[is.na(colors)] <- "grey"  # Assign grey to missing values

# Plot
quartz()
plot(sub_tree, cex = 0.8, tip.color = colors)

# Convert the phylogenetic tree into a ggtree plot
tree_plot <- ggtree(sub_tree) + 
  geom_tiplab(size = 3, color = colors) + 
  xlim(c(0,1.2)) # Add tip labels with assigned colors
  theme_tree2()  # Better tree visualization

# Arrange the PCA plot (pca_fig) and the tree (tree_plot) side by side
combined_plot <- pca_fig + tree_plot + 
  plot_layout(ncol = 2, widths = c(1, 1))  # Adjust width if needed

# Label and arrange the plots side by side
final_plot <- plot_grid(pca_fig, tree_plot, labels = c('A', 'B'), label_size = 12)

# Save the plot to a folder (adjust the file path as needed)
ggsave("/Users/damlacinoglu/Desktop/code and data/figures/fig_1B.png", final_plot, width = 10, height = 6)

#################################################
#Make table
#################################################
final_dat = new_coordinates
colnames(final_dat)[1] = "ATTRIBUTE_SPECIES"

final_dat = cbind(final_dat, genu_species = rep(NA, nrow(final_dat)), num_ind_seedling = rep(NA, nrow(final_dat)), num_ind_adult = rep(NA, nrow(final_dat)))
for(i in 1:nrow (final_dat)) {
	final_dat[i,"genu_species"] = meta[meta[,"ATTRIBUTE_SPECIES"] == final_dat[i,1],"genu_species"][1]
	final_dat[i,"num_ind_seedling"] = meta[meta[,"ATTRIBUTE_SPECIES"] == final_dat[i,1] &meta[,"ATTRIBUTE_Ontogeny"] == "Seedling" ,"num_ind"][1]
	final_dat[i,"num_ind_adult"] = meta[meta[,"ATTRIBUTE_SPECIES"] == final_dat[i,1] &meta[,"ATTRIBUTE_Ontogeny"] == "Adult" ,"num_ind"][1]
}

final_dat = cbind(final_dat, num_ind_herbivory = rep(NA, nrow(final_dat)))
for(i in 1:nrow (final_dat)) {
	final_dat[i,"num_ind_herbivory"] = summary_data[summary_data[,"SP"] == final_dat[i,1], "n"]
}

sla_data = read.csv("SLA_data_2023.csv", header = T)
final_dat = cbind(final_dat, num_ind_seedling_LMA= rep(NA, nrow(final_dat)), num_leaf_seedling_LMA = rep(NA, nrow(final_dat)))
for(i in 1:nrow (final_dat)) {
	final_dat[i,"num_ind_seedling_LMA"] = nrow(sla_data[sla_data[,"SPECIES"] == final_dat[i,1], ])
		final_dat[i,"num_leaf_seedling_LMA"] = length(unique(sla_data[sla_data[,"SPECIES"] == final_dat[i,1], "INDIVIDUAL"]))
}

writexl ::write_xlsx(final_dat, "final_dat.xlsx")

####################################################################################
####################################################################################
####################################################################################
#ADULTS
####################################################################################
####################################################################################
####################################################################################

meta_adult = meta_adult[!is.na(meta_adult $richness)& meta_adult $richness > 0 &!is.na(meta_adult $ATTRIBUTE_PC1score)&
!is.na(meta_adult $ATTRIBUTE_PC2score)&
!is.na(meta_adult $genu_species)
,]

meta_adult = meta_adult[meta_adult[,"ATTRIBUTE_SPECIES"] %in% species_new , ]
length(unique(meta_adult[,"ATTRIBUTE_SPECIES"]) %in% species_new)

model2 <- lmer(log(richness) ~ as.factor(ATTRIBUTE_Group)+ (1|num_ind), 
             data = meta_adult)
             
emm1.1 <- emmeans(model2, specs = pairwise ~ as.factor(ATTRIBUTE_Group), type = "response", adjust = "tukey")
 multcomp::cld(emm1.1$emmeans, Letters = letters, adjust = "tukey")
 
simulation_output <- simulateResiduals(fittedModel = model2, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)

# Step 2: Extract the residuals from the model
residuals_model2 <- residuals(model2)

# Step 3: Assess phylogenetic signal in residuals
# Read the phylogenetic tree
library(ape)
Tree <- read.tree("bci_phylogeny/BCI_Trees_Lianas_3gene.1.tre")

residuals_data <- data.frame(species = meta_adult $genu_species, residuals = residuals_model2)

Tree_filtered <- drop.tip(Tree, setdiff(Tree$tip.label, unique(residuals_data[,"species"])))
residuals_data2 = residuals_data[residuals_data[,"species"] %in% Tree_filtered$tip.label,]
unique(residuals_data2$species)

names(residuals_model2) = residuals_data[,"species"]
lambda <- phylosig(Tree_filtered, residuals_model2,method="K", test=TRUE, nsim=999)
print(lambda)
  
   summary_data <- meta_adult %>%
  dplyr::group_by(ATTRIBUTE_Group) %>%
  dplyr::summarize(
    mean_richness = mean(richness, na.rm = TRUE),
    se_richness = sd(richness, na.rm = TRUE) / sqrt(n())
  )
  
fig1 <- ggplot() +  
  theme_bw() + 
  theme(axis.ticks.x = element_blank()) +
  labs(y = "Secondary metabolomic compound\nrichness\n(Number of unique compounds)", 
       x = "Demographic group") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) + 
#  geom_point(data = meta_adult, aes(x = as.numeric(ATTRIBUTE_Group), y = richness), color = "gray", position = position_jitter(width = 0.1)) +  
  geom_pointrange(
    data = summary_data, 
    aes(x = as.numeric(ATTRIBUTE_Group), y = mean_richness, ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
    color = "black", size = 0.4
  ) +
  scale_x_continuous(
    breaks = 1:length(unique(meta_new$ATTRIBUTE_Group)),
    labels = c("Slow", "Fast", "LLP","SLB","Int")  # <-- customize these!
  )  + coord_cartesian(ylim = c(90, 150)) 
  
#Metabolomic richness

model1 = lmer(hill_div ~ as.factor(ATTRIBUTE_Group) + (1|num_ind) , data = meta_adult)
emm1.1 <- emmeans(model1, specs = pairwise ~ as.factor(ATTRIBUTE_Group), type = "response", adjust = "tukey")
 multcomp::cld(emm1.1$emmeans, Letters = letters, adjust = "tukey")

simulation_output <- simulateResiduals(fittedModel = model1, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)   

# Step 3: Assess phylogenetic signal in residuals
# Read the phylogenetic tree
Tree <- read.tree("bci_phylogeny/BCI_Trees_Lianas_3gene.1.tre")

residuals_data <- data.frame(species = meta_adult $genu_species, residuals = residuals_model2)

Tree_filtered <- drop.tip(Tree, setdiff(Tree$tip.label, unique(residuals_data[,"species"])))
residuals_data2 = residuals_data[residuals_data[,"species"] %in% Tree_filtered$tip.label,]
unique(residuals_data2$species)

names(residuals_model2) = residuals_data[,"species"]
lambda <- phylosig(Tree_filtered, residuals_model2,method="K", test=TRUE, nsim=999)
print(lambda)
 
    summary_data <- meta_adult %>%
  dplyr::group_by(ATTRIBUTE_Group) %>%
  dplyr::summarize(
    mean_richness = mean(hill_div, na.rm = TRUE),
    se_richness = sd(hill_div, na.rm = TRUE) / sqrt(n())
  )
  
fig2 <- ggplot() +  
  theme_bw() + 
  theme(axis.ticks.x = element_blank()) +
  labs(y = "Functional hill diversity\n(Sum of the pairwise dissimilarities\n in the compound dissimilarity matrix)", 
       x = "Demographic group") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) + 
 #geom_point(data = meta_adult, aes(x = as.numeric(ATTRIBUTE_Group), y = hill_div), color = "gray", position = position_jitter(width = 0.1)) +  
  geom_pointrange(
    data = summary_data, 
    aes(x = as.numeric(ATTRIBUTE_Group), y = mean_richness, ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
    color = "black", size = 0.4
  ) +
  scale_x_continuous(
    breaks = 1:length(unique(meta_new$ATTRIBUTE_Group)),
    labels = c("Slow", "Fast", "LLP","SLB","Int")  
  )  + coord_cartesian(ylim = c(2500, 7000)) 

####################################################################################
####################################################################################
####################################################################################

sla_data = read.csv("lma_data.csv", header = T)

sla_data$genu_species = paste0(sla_data$GENUS. , "_", sla_data$SPECIES.)

# Filter SLA data to species present in meta_adult
mean_lma_by_species_indiv <- sla_data %>%
  dplyr::filter(genu_species %in% unique(meta_adult$genu_species))

# Get distinct PC scores per species from meta
meta_pc_scores <- meta %>%
  dplyr::select(genu_species, ATTRIBUTE_PC1score, ATTRIBUTE_PC2score, ATTRIBUTE_Group) %>%
  dplyr::distinct()

# Join the two
mean_lma_annotated <- mean_lma_by_species_indiv %>%
  dplyr::left_join(meta_pc_scores, by = "genu_species")

mean_lma_unique <- mean_lma_annotated %>%
  dplyr::distinct(genu_species, .keep_all = TRUE)
  
 unique(mean_lma_unique[,"genu_species"]) %in% unique(meta_new[ meta_new[,"ATTRIBUTE_SPECIES"] %in% species_new , "genu_species"])
  
model2 = lm(LMA/1000 ~ as.factor(ATTRIBUTE_Group) , data = mean_lma_unique[mean_lma_unique[,"genu_species"]%in% unique(meta_new[ meta_new[,"ATTRIBUTE_SPECIES"] %in% species_new , "genu_species"]),])
emm1.1 <- emmeans(model2, specs = pairwise ~ as.factor(ATTRIBUTE_Group), type = "link", adjust = "tukey")
 multcomp::cld(emm1.1$emmeans, Letters = letters, adjust = "tukey")

simulation_output <- simulateResiduals(fittedModel = model2, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)

residuals_model2 <- residuals(model2)

# Step 3: Assess phylogenetic signal in residuals
# Read the phylogenetic tree
library(ape)
Tree <- read.tree("bci_phylogeny/BCI_Trees_Lianas_3gene.1.tre")

residuals_data <- data.frame(species = mean_lma_unique $genu_species, residuals = residuals_model2)

Tree_filtered <- drop.tip(Tree, setdiff(Tree$tip.label, unique(residuals_data[,"species"])))
residuals_data2 = residuals_data[residuals_data[,"species"] %in% Tree_filtered$tip.label,]
unique(residuals_data2$species)

names(residuals_model2) = residuals_data[,"species"]
lambda <- phylosig(Tree_filtered, residuals_model2,method="K", test=TRUE, nsim=999)
print(lambda)

   summary_data <- mean_lma_unique %>%
 dplyr:: group_by(ATTRIBUTE_Group) %>%
  dplyr:: summarize(
    mean_richness = mean(LMA/1000, na.rm = TRUE),
    se_richness = sd(LMA/1000, na.rm = TRUE) / sqrt(n())
  )
  
fig3 <- ggplot() +  
  theme_bw() + 
  theme(axis.ticks.x = element_blank()) +
  labs(y = "Leaf mass per unit area (kg/m2)", 
       x = "Demographic group") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) + 
 # geom_point(data = mean_lma_unique, aes(x = as.numeric(ATTRIBUTE_Group), y = LMA), color = "gray", position = position_jitter(width = 0.1)) +  
  # Mean and error bars
  geom_pointrange(
    data = summary_data, 
    aes(x = as.numeric(ATTRIBUTE_Group), y = mean_richness, ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
    color = "black", size = 0.4
  ) +
  scale_x_continuous(
    breaks = 1:length(unique(meta_new$ATTRIBUTE_Group)),
    labels = c("Slow", "Fast", "LLP","SLB","Int")  # <-- customize these!
  )  + coord_cartesian(ylim = c(0.040, 0.08)) 

combined_plot <- plot_grid(fig1, fig2, 
                           fig3,
                           ncol = 3, align = "hv", labels = c('A', 'B',"C"), label_size = 12)

ggsave("/Users/damlacinoglu/Desktop/code and data/figures/fig_S1**.png", combined_plot, width = 10, height = 3.5)


####################################################################################
####################################################################################
####################################################################################
#SEEDLINGS
####################################################################################
####################################################################################
####################################################################################

meta_new = meta[!is.na(meta$richness)& meta$richness > 0 &!is.na(meta$ATTRIBUTE_PC1score)&
!is.na(meta$ATTRIBUTE_PC2score)&
!is.na(meta$genu_species)&
!is.na(meta$ATTRIBUTE_LEAF)&
!is.na(meta$ATTRIBUTE_INDIVIDUAL),]
#meta_new[,"genu_species"] = gsub("_Seedling","", meta_new[,"genu_species"])
 #meta_new = meta_new[,-49]
 
#meta_new = meta_new[meta_new[,"genu_species"] %in% species_to_keep,]
unique(meta_new[,"genu_species"])

unique(meta_new[,"ATTRIBUTE_SPECIES"]) %in% species_new

model4 <- lmer(richness ~ as.factor(ATTRIBUTE_Group) + ATTRIBUTE_LEAF  +  (1 | ATTRIBUTE_SPECIES:ATTRIBUTE_INDIVIDUAL), 
             data = meta_new)
             
           model =   aov(richness ~ as.factor(ATTRIBUTE_Group) + ATTRIBUTE_SPECIES, meta_new)
emm1.1 <- emmeans(model4, specs = pairwise ~ as.factor(ATTRIBUTE_Group), type = "response", adjust = "tukey")
 multcomp::cld(emm1.1$emmeans, Letters = letters, adjust = "tukey")

simulation_output <- simulateResiduals(fittedModel = model2, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)

# Step 2: Extract the residuals from the model
residuals_model2 <- residuals(model2)

# Step 3: Assess phylogenetic signal in residuals
# Read the phylogenetic tree
library(ape)
Tree <- read.tree("bci_phylogeny/BCI_Trees_Lianas_3gene.1.tre")

residuals_data <- data.frame(species = meta_new $genu_species, residuals = residuals_model2)

Tree_filtered <- drop.tip(Tree, setdiff(Tree$tip.label, unique(residuals_data[,"species"])))
residuals_data2 = residuals_data[residuals_data[,"species"] %in% Tree_filtered$tip.label,]
unique(residuals_data2$species)

names(residuals_model2) = residuals_data[,"species"]
lambda <- phylosig(Tree_filtered, residuals_model2,method="K", test=TRUE, nsim=999)
print(lambda)
  
 summary_data <- meta_new %>%
  dplyr::group_by(ATTRIBUTE_Group) %>%
  dplyr::summarize(
    mean_richness = mean(richness, na.rm = TRUE),
    se_richness = sd(richness, na.rm = TRUE) / sqrt(n())
  )
  
fig1 <- ggplot() +  
  theme_bw() + 
  theme(axis.ticks.x = element_blank()) +
  labs(y = "Secondary metabolomic compound\nrichness\n(Number of unique compounds)", 
       x = "Demographic group") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) + 
  #geom_point(data = meta_new, aes(x = as.numeric(ATTRIBUTE_Group), y = richness), color = "gray", position = position_jitter(width = 0.1)) +  
  geom_pointrange(
    data = summary_data, 
    aes(x = as.numeric(ATTRIBUTE_Group), y = mean_richness, ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
    color = "black", size = 0.4
  ) +
  scale_x_continuous(
    breaks = 1:length(unique(meta_new$ATTRIBUTE_Group)),
    labels = c("Slow", "Fast", "LLP","SLB","Int")  # <-- customize these!
  ) + coord_cartesian(ylim = c(100, 180)) 

####################################################################################
####################################################################################
####################################################################################

meta_new = meta[!is.na(meta$hill_div)& meta$hill_div > 0 &!is.na(meta$ATTRIBUTE_PC1score)&
!is.na(meta$ATTRIBUTE_PC2score)&
!is.na(meta$genu_species)&
!is.na(meta$ATTRIBUTE_LEAF)&
!is.na(meta$ATTRIBUTE_INDIVIDUAL),]

unique(meta_new[,"ATTRIBUTE_SPECIES"]) %in% species_new

model3 <- lmer(hill_div ~ as.factor(ATTRIBUTE_Group) + ATTRIBUTE_LEAF + 
                 (1 | ATTRIBUTE_SPECIES:ATTRIBUTE_INDIVIDUAL) , 
             data = meta_new)

 simulation_output <- simulateResiduals(fittedModel = model3, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)

 emm1.1 <- emmeans(model3, specs = pairwise ~ as.factor(ATTRIBUTE_Group), type = "response", adjust = "Tukey")
 multcomp::cld(emm1.1$emmeans, Letters = letters, adjust = "Tukey")

# Step 2: Extract the residuals from the model
residuals_model2 <- residuals(model2)

# Step 3: Assess phylogenetic signal in residuals
# Read the phylogenetic tree
library(ape)
Tree <- read.tree("bci_phylogeny/BCI_Trees_Lianas_3gene.1.tre")

residuals_data <- data.frame(species = meta_new $genu_species, residuals = residuals_model2)

Tree_filtered <- drop.tip(Tree, setdiff(Tree$tip.label, unique(residuals_data[,"species"])))
residuals_data2 = residuals_data[residuals_data[,"species"] %in% Tree_filtered$tip.label,]
unique(residuals_data2$species)

names(residuals_model2) = residuals_data[,"species"]
lambda <- phylosig(Tree_filtered, residuals_model2,method="K", test=TRUE, nsim=999)
print(lambda)
  
 summary_data <- meta_new %>%
 dplyr:: group_by(ATTRIBUTE_Group) %>%
 dplyr:: summarize(
    mean_richness = mean(hill_div, na.rm = TRUE),
    se_richness = sd(hill_div, na.rm = TRUE) / sqrt(n())
  )
  
fig2 <- ggplot() +  
  theme_bw() + 
  theme(axis.ticks.x = element_blank()) +
  labs(y = "Functional hill diversity\n(Sum of the pairwise dissimilarities\n in the compound dissimilarity matrix)", 
       x = "Demographic group") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) + 
#  geom_point(data = meta_new, aes(x = as.numeric(ATTRIBUTE_Group), y = hill_div), color = "gray", position = position_jitter(width = 0.1)) +  
  geom_pointrange(
    data = summary_data, 
    aes(x = as.numeric(ATTRIBUTE_Group), y = mean_richness, ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
    color = "black", size = 0.4
  ) +
  scale_x_continuous(
    breaks = 1:length(unique(meta_new$ATTRIBUTE_Group)),
    labels = c("Slow", "Fast", "LLP","SLB","Int")    )  + coord_cartesian(ylim = c(4200, 11500)) 
      
####################################################################################
####################################################################################
####################################################################################

sla_data = read.csv("SLA_data_2023.csv", header = T)

sla_data = cbind(sla_data, mean_LMA_DAMLA_kg_m2 = rep(NA, nrow(sla_data)))
sla_data$mean_LMA_DAMLA_kg_m2 = as.numeric(sla_data$LMA_DAMLA) *10
sla_data_annotated <- sla_data %>%
  dplyr::left_join(meta, by = c("SPECIES" = "ATTRIBUTE_SPECIES"))
sla_data_annotated = sla_data_annotated[sla_data_annotated[,"SPECIES"] %in% species_new,]

#40 SPECIES INCLUDED!
unique(sla_data_annotated[,"SPECIES"])

model2 <- lmer(log(mean_LMA_DAMLA_kg_m2) ~ as.factor(ATTRIBUTE_Group) + (1|SPECIES: INDIVIDUAL ), 
             data = sla_data_annotated)
emm1.1 <- emmeans(model2, specs = pairwise ~ as.factor(ATTRIBUTE_Group), type = "response", adjust = "tukey")
 multcomp::cld(emm1.1$emmeans, Letters = letters, adjust = "tukey")
 
simulation_output <- simulateResiduals(fittedModel = model2, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)
   
# Step 2: Extract the residuals from the model
residuals_model2 <- residuals(model2)

# Step 3: Assess phylogenetic signal in residuals
# Read the phylogenetic tree
library(ape)
Tree <- read.tree("bci_phylogeny/BCI_Trees_Lianas_3gene.1.tre")

residuals_data <- data.frame(species = mean_lma_filtered $genu_species, residuals = residuals_model2)

Tree_filtered <- drop.tip(Tree, setdiff(Tree$tip.label, unique(residuals_data[,"species"])))
residuals_data2 = residuals_data[residuals_data[,"species"] %in% Tree_filtered$tip.label,]
unique(residuals_data2$species)

names(residuals_model2) = residuals_data[,"species"]
lambda <- phylosig(Tree_filtered, residuals_model2,method="K", test=TRUE, nsim=999)
print(lambda)

     summary_data <- sla_data_annotated %>%
 dplyr:: group_by(ATTRIBUTE_Group) %>%
 dplyr:: summarize(
    mean_richness = mean(mean_LMA_DAMLA_kg_m2, na.rm = TRUE),
    se_richness = sd(mean_LMA_DAMLA_kg_m2, na.rm = TRUE) / sqrt(n())
  )
  
fig3 <- ggplot() +  
  theme_bw() + 
  theme(axis.ticks.x = element_blank()) +
  labs(y = "Leaf mass per unit area (kg/m2)", 
       x = "Demographic group") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) + 
#  geom_point(data = mean_lma_annotated, aes(x = as.numeric(ATTRIBUTE_Group), y = mean_LMA_DAMLA), color = "gray", position = position_jitter(width = 0.1)) +  
  geom_pointrange(
    data = summary_data, 
    aes(x = as.numeric(ATTRIBUTE_Group), y = mean_richness, ymin = mean_richness - se_richness, ymax = mean_richness + se_richness),
    color = "black", size = 0.4
  ) +
  scale_x_continuous(
    breaks = 1:length(unique(meta_new$ATTRIBUTE_Group)),
    labels = c("Slow", "Fast", "LLP","SLB","Int")   )  + coord_cartesian(ylim = c(0.035, 0.055)) 
      

combined_plot <- plot_grid(fig1, fig2, 
                           fig3,
                           ncol = 3, align = "hv", labels = c('A', 'B',"C"), label_size = 12)

ggsave("/Users/damlacinoglu/Desktop/code and data/figures/fig_3.png", combined_plot, width = 10, height = 3.5)

#################################################
#Herbivory analyses
#################################################

meta_pc_scores <- meta %>%
  filter(ATTRIBUTE_Ontogeny == "Seedling") %>%
  distinct(avg_herbivory, .keep_all = TRUE)
hist(log(meta_pc_scores$avg_herbivory))
hist(meta_pc_scores$hill_div_sp)

hist(meta[meta["ATTRIBUTE_Ontogeny"] == "Seedling" ,"avg_herbivory"])
hist(meta[meta["ATTRIBUTE_Ontogeny"] == "Seedling" ,"hill_div"])
hist(meta[meta["ATTRIBUTE_Ontogeny"] == "Seedling" ,"richness"])

 model3 <- lmer(
log(hill_div) ~  log(avg_herbivory) + ATTRIBUTE_LEAF  + (1| ATTRIBUTE_SPECIES:ATTRIBUTE_INDIVIDUAL) ,
    data = meta[meta[,"ATTRIBUTE_Ontogeny"] == "Seedling" & meta[,"richness"] > 0 & meta[,"ATTRIBUTE_SPECIES"] %in% species_new,]
)
 model5 <- lmer(
log(richness) ~  log(avg_herbivory) + ATTRIBUTE_LEAF + (1| ATTRIBUTE_SPECIES:ATTRIBUTE_INDIVIDUAL) ,
  data = meta[meta[,"ATTRIBUTE_Ontogeny"] == "Seedling" & meta[,"richness"] > 0 & meta[,"ATTRIBUTE_SPECIES"] %in% species_new,])
 
summary(model3)
summary(model5)

simulation_output <- simulateResiduals(fittedModel = v, plot = TRUE)
testUniformity(simulation_output)
testDispersion(simulation_output)
testOutliers(simulation_output)

###########################################################
# Subset and drop unused levels
meta_subset <- droplevels(meta[meta$ATTRIBUTE_Ontogeny == "Seedling" & meta$richness > 0, ])

# Fit the model using the cleaned subset
model4 <- lmer(
  log(richness) ~ log(avg_herbivory) + ATTRIBUTE_LEAF  + (1| ATTRIBUTE_SPECIES:ATTRIBUTE_INDIVIDUAL) ,
  data = meta_subset
)

# Create sequence of avg_herbivory values for prediction
herb_seq <- seq(
  from = min(meta_subset$avg_herbivory, na.rm = TRUE),
  to = max(meta_subset$avg_herbivory, na.rm = TRUE),
  length.out = 100
)

# Use only levels present in the subset
grid <- expand.grid(
  avg_herbivory = herb_seq,
  ATTRIBUTE_LEAF = unique(meta_subset$ATTRIBUTE_LEAF))

# Predictions
pred <- marginaleffects::predictions(model4, newdata = grid, re.form = NA)

# Summarize predictions
new_data_avg <- pred %>%
  dplyr::group_by(avg_herbivory) %>%
  dplyr::summarize(
    predicted_p = exp(mean(estimate, na.rm = TRUE)),
    predicted_p_lower = exp(mean(conf.low, na.rm = TRUE)),
    predicted_p_upper = exp(mean(conf.high, na.rm = TRUE)))

plot1 <- ggplot(meta[meta$ATTRIBUTE_Ontogeny == "Seedling" & meta$richness > 0, ], 
                aes(x = avg_herbivory, y = richness)) +
  geom_point(size = 1, color = "grey") +
  geom_line(data = new_data_avg, aes(x = avg_herbivory, y = predicted_p), color = "blue", size = 1, inherit.aes = FALSE) +
  geom_ribbon(data = new_data_avg, aes(x = avg_herbivory, ymin = predicted_p_lower, ymax = predicted_p_upper),
              fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
  labs(
    x = "Mean herbivory per species (%)",
    y = "Secondary metabolomic compound richness\n(Number of unique compounds)",
    title = " "
  ) +  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) + 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12))   + theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 
        
###########################################################      
# Fit the model using the cleaned subset
model4 <- lmer(
  log(hill_div) ~ log(avg_herbivory) + ATTRIBUTE_LEAF + (1| ATTRIBUTE_SPECIES:ATTRIBUTE_INDIVIDUAL),
  data = meta_subset[meta_subset[,"hill_div"] < 20000 & meta_subset[,"avg_herbivory"] > 1.5,]
)

# Create sequence of avg_herbivory values for prediction
herb_seq <- seq(
  from = 1.5,
  to = max(meta_subset$avg_herbivory, na.rm = TRUE),
  length.out = 100
)

# Use only levels present in the subset
grid <- expand.grid(
  avg_herbivory = herb_seq,
  ATTRIBUTE_LEAF = unique(meta_subset$ATTRIBUTE_LEAF)
)

# Predictions
pred <- marginaleffects::predictions(model4, newdata = grid, re.form = NA)

# Summarize predictions
new_data_avg <- pred %>%
  dplyr::group_by(avg_herbivory) %>%
  dplyr::summarize(
    predicted_p = exp(mean(estimate, na.rm = TRUE)),
    predicted_p_lower = exp(mean(conf.low, na.rm = TRUE)),
    predicted_p_upper = exp(mean(conf.high, na.rm = TRUE))
  )

plot2 <- ggplot(meta[meta$ATTRIBUTE_Ontogeny == "Seedling" & meta$hill_div <20000 & meta$avg_herbivory>1.5, ], 
                aes(x = avg_herbivory, y = hill_div)) +
  geom_point(size = 1, color = "grey") +
  geom_line(data = new_data_avg, aes(x = avg_herbivory, y = predicted_p), color = "blue", size = 1, inherit.aes = FALSE) +
  geom_ribbon(data = new_data_avg, aes(x = avg_herbivory, ymin = predicted_p_lower, ymax = predicted_p_upper),
              fill = "blue", alpha = 0.2, inherit.aes = FALSE) +
  labs(
    x = "Mean herbivory per species (%)",
    y = "Functional hill diversity\n(Sum of the pairwise dissimilarities\n in the compound dissimilarity matrix)",
    title = " "
  ) +  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) + 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12)
  ) +  coord_cartesian(ylim = c(0, 20000)) +   theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

###########################################################      
 model3 <- lmer(
log(hill_div) ~  log(avg_herbivory) + ATTRIBUTE_LEAF  + (1| ATTRIBUTE_SPECIES:ATTRIBUTE_INDIVIDUAL) ,
    data = meta[meta[,"ATTRIBUTE_Ontogeny"] == "Seedling" & meta[,"richness"] > 0 & meta[,"ATTRIBUTE_SPECIES"] %in% species_new,]
)
model3_slopes <- lmer(
  log(hill_div) ~ log(avg_herbivory) + ATTRIBUTE_LEAF + 
    (1 + log(avg_herbivory) | ATTRIBUTE_SPECIES),
  data = meta_subset
)
 model5 <- lmer(
log(richness) ~  log(avg_herbivory) + ATTRIBUTE_LEAF + (1| ATTRIBUTE_SPECIES:ATTRIBUTE_INDIVIDUAL) ,
  data = meta[meta[,"ATTRIBUTE_Ontogeny"] == "Seedling" & meta[,"richness"] > 0 & meta[,"ATTRIBUTE_SPECIES"] %in% species_new,])
 model3_slopes <- lmer(
  log(richness) ~ log(avg_herbivory) + ATTRIBUTE_LEAF + 
    (1 + log(avg_herbivory) | ATTRIBUTE_SPECIES),
  data = meta_subset
)
 
data = coef(model3)$ATTRIBUTE_SPECIES

species_group_lookup <- unique(meta[, c("ATTRIBUTE_SPECIES", "ATTRIBUTE_Group", "ATTRIBUTE_PC1score","ATTRIBUTE_PC2score")])
data$ATTRIBUTE_SPECIES <- rownames(data)
merged_data <- merge(data, species_group_lookup, by = "ATTRIBUTE_SPECIES", all.x = TRUE)
colnames(merged_data)[3] = "herb"

 model4 <- lm(herb  ~  ATTRIBUTE_Group , data = merged_data)
 summary(model4)
pairs(emmeans(model4, ~ ATTRIBUTE_Group), adjust = "Tukey")

species_coefs <- coef(model3)$ATTRIBUTE_SPECIES
species_se <- se.coef(model3)$ATTRIBUTE_SPECIES  # from arm package

df <- data.frame(
  ATTRIBUTE_SPECIES = rownames(species_coefs),
  herb = species_coefs[, "log(avg_herbivory)"],
  se = species_se[, "log(avg_herbivory)"]
)
df <- merge(df, species_group_lookup, by = "ATTRIBUTE_SPECIES")
df$weight <- 1 / (df$se^2)
model4_weighted <- lm(herb ~ ATTRIBUTE_Group, data = df, weights = weight)
summary(model4_weighted)
pairs(emmeans(model4_weighted, ~ ATTRIBUTE_Group), adjust = "Tukey")

meta$ATTRIBUTE_LEAF = as.factor(meta$ATTRIBUTE_LEAF )
meta$ATTRIBUTE_SPECIES = as.factor(meta$ATTRIBUTE_SPECIES )

plot1a <- ggplot(
  meta[meta$ATTRIBUTE_Ontogeny == "Seedling" & meta$richness > 0 & meta$ATTRIBUTE_Group %in% c("1","2","3","4","5"), ], 
  aes(x = avg_herbivory, y = richness, color = ATTRIBUTE_Group)
) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_manual(
    name = "Demographic Group",  # Legend title
    values = c(
      "1" = "pink", 
      "2" = "orange", 
      "3" = "green" ,
        "4" = "blue",
          "5" = "red"
    )
  ) +
  labs(
    x = "Mean herbivory per species (%)",
    y = "Secondary metabolomic compound richness\n(Number of unique compounds)",
    title = ""
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
 legend.position = "none")

plot2a <- ggplot(
  meta[meta$ATTRIBUTE_Ontogeny == "Seedling" & meta$hill_div < 20000 & meta$avg_herbivory > 1.5 & meta$ATTRIBUTE_Group %in% c("1", "2", "3", "4", "5"), ], 
  aes(x = avg_herbivory, y = hill_div, color = ATTRIBUTE_Group)
) +
  geom_point(size = 1) +
  geom_smooth(method = "lm", se = FALSE, size = 1) +
  scale_color_manual(
    name = "Demographic Group",
    values = c(
      "1" = "pink", 
      "2" = "orange", 
      "3" = "green",
      "4" = "blue",
      "5" = "red"
    ),
    labels = c(
      "1" = "Slow",
      "2" = "Fast",
      "3" = "LLP",
      "4" = "SLB",
      "5" = "Int"
    )
  ) +
  labs(
    x = "Mean herbivory per species (%)",
    y = "Functional hill diversity\n(Sum of the pairwise dissimilarities\n in the compound dissimilarity matrix",
    title = ""
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

final_plot <- plot_grid(plot1a, plot2a, labels = c('A', 'B'), label_size = 12, ncol = 2,rel_widths = c(1, 1.6))

ggsave("/Users/damlacinoglu/Desktop/code and data/figures/fig_herbivory_supp.png", final_plot, width = 9, height = 4)

#################################################
#LASSO regression
#################################################
#run cv.glmnet() once using 'herbivory' as the phenotype and then a few more times using Ruger group or PC1 or 2.
#x needs to be a matrix in which rows are tree species and columns are metabolites (or "traits" defined by "clades" in the chemical dendrogram); y for you will be a vector of herbivory or life history group or life history PC1 or 2.

#Using individual seedlings:

load("DC_ALL_Metabmastertable_240920.RData")

class = "secondary" 
if(class == "secondary"){
  heat.class = master.real[which(master.real$NPC.pathway %in% c("Alkaloids","Amino acids and Peptides","Polyketides","Shikimates and Phenylpropanoids", "Terpenoids")),]
}

nrow(heat.class)

colnames(heat.class) = gsub(".Peak.area", "", colnames(heat.class))
colnames(heat.class) = gsub("X", "", colnames(heat.class))
colnames(heat.class) = gsub("mzML", "mzXML", colnames(heat.class))
for(i in 1:ncol(heat.class)) {
  if(colnames(heat.class)[i] == "01.mzXML") {colnames(heat.class)[i] =  "1.mzXML"}
  if(colnames(heat.class)[i] == "02.mzXML") {colnames(heat.class)[i] =  "2.mzXML"}
  if(colnames(heat.class)[i] == "03.mzXML") {colnames(heat.class)[i] =  "3.mzXML"}
  if(colnames(heat.class)[i] == "04.mzXML") {colnames(heat.class)[i] =  "4.mzXML"}
  if(colnames(heat.class)[i] == "05.mzXML") {colnames(heat.class)[i] =  "5.mzXML"}
  if(colnames(heat.class)[i] == "06.mzXML") {colnames(heat.class)[i] =  "6.mzXML"}
  if(colnames(heat.class)[i] == "07.mzXML") {colnames(heat.class)[i] =  "7.mzXML"}
  if(colnames(heat.class)[i] == "08.mzXML") {colnames(heat.class)[i] =  "8.mzXML"}
  if(colnames(heat.class)[i] == "09.mzXML") {colnames(heat.class)[i] =  "9.mzXML"}
}

x_matrix = heat.class
rownames(heat.class) = heat.class[,"row.ID"]


x_matrix = x_matrix[,19:ncol(x_matrix)]

y_matrix = meta[!is.na(meta[,"avg_herbivory"]),]

#ONLY THE 42 SPECIES INCLUDED!!!!
y_matrix = y_matrix[y_matrix[,"ATTRIBUTE_SPECIES"] %in% species_new,]
unique(y_matrix[,"ATTRIBUTE_SPECIES"])

filtered_filenames <- y_matrix[!grepl("BCI", y_matrix[,"filename"]), ]
x_matrix = x_matrix[,colnames(x_matrix) %in% filtered_filenames[,"filename"],]

row.names(y_matrix) = y_matrix[,"filename"]
y_matrix = y_matrix[colnames(x_matrix),]

#these two should be equal to each other!
unique(length(colnames(x_matrix)))
nrow(y_matrix)

cvfit = cv.glmnet(x = t(x_matrix), y = t(y_matrix[,"avg_herbivory"]))

coef = coef(cvfit, s = "lambda.min")
compeffs = coef@Dimnames[[1]][(coef@i+1)]

coef_df = coef_df[rownames(coef_df) != "(Intercept)", ]

new_metab = heat.class[rownames(heat.class) %in% coef_df$Metabolite, ]
length(unique(compeffs)); nrow(new_metab)

coef_values = as.matrix(coef(cvfit, s = "lambda.min"))  # Convert coefficients to a matrix
coef_df = data.frame(Metabolite = rownames(coef_values), Coefficient = coef_values[,1])  
coef_df = coef_df[coef_df$Coefficient != 0, ]  # Keep only nonzero coefficients
positive_influencers = coef_df[coef_df$Coefficient > 0, ]
negative_influencers = coef_df[coef_df$Coefficient < 0, ]

nrow(positive_influencers)
nrow(negative_influencers)

positive_metab = heat.class[heat.class[,"row.ID"] %in% positive_influencers$Metabolite,]
negative_metab = heat.class[heat.class[,"row.ID"] %in% negative_influencers$Metabolite,]
positive_influencers = positive_influencers[positive_influencers[,1] %in% positive_metab[,"row.ID"],]
negative_influencers = negative_influencers[negative_influencers[,1] %in% negative_metab[,"row.ID"],]

coef_df_final = rbind(positive_influencers, negative_influencers) 
for(i in 1:nrow(coef_df_final)) {
	coef_df_final[i,"smile"] = heat.class[rownames(heat.class) == as.numeric(coef_df_final[i,1]) , "smiles"]
		coef_df_final[i,"pathway"] = heat.class[heat.class[,"row.ID"] == coef_df_final[i,1] , "NPC.pathway"     ]  
}

writexl ::write_xlsx(coef_df_final, "coef_df_final_final_dat2.xlsx")

pathway_counts <- as.data.frame(table(negative_metab[,"NPC.pathway"]))
for(i in 1:nrow(pathway_counts)) {
pathway_counts[i,3] = 	pathway_counts[i,2] / nrow(heat.class[heat.class[,"NPC.pathway"] == pathway_counts[i,1],])
}

negative_plot = ggplot(pathway_counts, aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_text(aes(label = Var1, y = 0.5), angle = 90, hjust = 0, vjust = 0.5, size = 6) +  # Start from x-axis
  theme_bw() +
  labs(title = "Pathways downregulating herbivory",
       x = NULL,  # Remove x-axis label
       y = "Count") +
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank()) +  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

pathway_counts <- as.data.frame(table(positive_metab[,"NPC.pathway"]))
for(i in 1:nrow(pathway_counts)) {
pathway_counts[i,3] = 	pathway_counts[i,2] / nrow(heat.class[heat.class[,"NPC.pathway"] == pathway_counts[i,1],])
}

positive_plot = ggplot(pathway_counts, aes(x = reorder(Var1, -Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_text(aes(label = Var1, y = 0.5), angle = 90, hjust = 0, vjust = 0.5, size = 6) +  # Start from x-axis
  theme_bw() +
  labs(title = "Pathways upregulating herbivory",
       x = NULL,  # Remove x-axis label
       y = "Count") +
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank()) +   theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

final_plot <- plot_grid(plot1, plot2, negative_plot, positive_plot, labels = c('A', 'B',"C", "D"), label_size = 12, ncol = 2)

ggsave("/Users/damlacinoglu/Desktop/code and data/figures/fig_herbivory.png", final_plot, width = 7, height = 10)


 
