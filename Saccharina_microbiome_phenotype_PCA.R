setwd("/Volumes/Artemis/Kelp Microbiome/Saccharina Sugar Kelp")
library(phyloseq)
library(dplyr)
library(tidyr)
library(microbiome)
library(vegan)
library(FactoMineR)
library(Factoshiny)

## This script takes the saccharina phenotype data, metaxa2 output (SILVA128 SSU classification), and performs exploratory factor analysis.

## Remember: Grab abundance table from /project/noujdine_61/mgosborn/Latissima/Biomass_Subset/level_7_output
## Command: metaxa2_dc -o AbundanceTable.txt -r "metaxa2_" -p "^[^.]+" *level_7.txt

otumat <- read.delim("AbundanceTable.txt", sep = "\t", row.names = 1, check.names = FALSE)
taxmat <- matrix(nrow = nrow(otumat), ncol = 0)
rownames(taxmat) <- rownames(otumat)
taxmat <- as.data.frame(cbind(Taxon = rownames(taxmat), taxmat))
levels <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
taxmat <- taxmat %>% separate(Taxon, levels, sep = ";")
taxmat <- as.matrix(taxmat)

#create phyloseq objects
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)

#Only keep Bacteria
physeq = subset_taxa(physeq, Domain=="Bacteria")

#### REMOVE SINGLETONS & DOUBLETONS
physeq <- prune_taxa(taxa_sums(physeq) > 2, physeq) 


#### CONSTRUCTING METADATA 
names <- as.data.frame(sample_names(physeq))
colnames(names)[1] <- "SampleName" #rename first column
names$fullID <- names$SampleName
levels <- c("SampleID","Index","Lane")
names <- names %>% separate(SampleName, levels, sep = "_") #Split Names by '_'

pheno <- read.csv("sub60_for_microbiome.csv", header = TRUE)
pheno <- pheno[ , -which(names(pheno) %in% c("X","Name_InGeno","Name_InPheno"))]
names(pheno)[names(pheno) == "Sample.ID_in_FastQ"] <- "SampleID"

metadata <- merge(pheno, names, by = "SampleID", no.dups = FALSE)
rownames(metadata) <- metadata$fullID
temp = metadata[ , -which(names(metadata) %in% c("Index","Lane","fullID"))]
metadata_merged = unique(temp)
rownames(metadata_merged) <- metadata_merged$SampleID

#add population info to physeq object
physeq@sam_data = sample_data(metadata)

#Standardize abundances to the median sequencing depth
total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
temp = transform_sample_counts(physeq,standf)

#Merge samples from multiple lanes
physeq_merged = merge_samples(temp, "SampleID",fun = mean)
# For samples with multiple files, divide by number of times sample occurs in dataset
temp_otu <- t(as.data.frame(physeq_merged@otu_table))
temp2 <- as.data.frame(table(metadata$SampleID))
temp2 <- subset(temp2, temp2$Freq>1)
row.names(temp2) <- temp2$Var1
for (z in temp2$Var1) {
  temp_otu[,z] <- round((temp_otu[,z])/temp2[z,"Freq"])
}
temp_otu <- t(temp_otu)

OTU = otu_table(temp_otu, taxa_are_rows = FALSE)
physeq_merged@otu_table = OTU
sample_data(physeq_merged) <- metadata_merged

#### REMOVE SINGLETONS & DOUBLETONS
physeq_merged <- prune_taxa(taxa_sums(physeq_merged) > 2, physeq_merged) 

### ALPHA DIVERSITY
dev.new()
plot_richness(physeq_merged, x="Select", measures=c("Chao1", "Shannon"))

## BAR PLOT
dev.new()
plot_bar(physeq_merged)


##### EXPLORATORY FACTOR ANALYSIS - FactoMineR
# center log transform
physeq_merged_transformed <- microbiome::transform(physeq_merged, 'clr')

# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(physeq_merged_transformed@otu_table)
sub_metadata <- subset(metadata_merged, select = c("Select","Location","Parent"))

# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,sub_metadata, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
facto_mat <- na.omit(facto_mat)

### FactoMineR
res.pca <- PCA(facto_mat, quali.sup=1168:1170, graph = FALSE)
res.shiny=PCAshiny(res.pca)


######### Now do the same for only top 50% most abundant OTUs
PM_top50 = prune_taxa(names(sort(taxa_sums(physeq_merged), TRUE))[1:round(ncol(physeq_merged@otu_table)/2)], physeq_merged)

# center log transform
PM_top50_transformed <- microbiome::transform(PM_top50, 'clr')

# get otu, tax table, and metadata from physeq
otumat <- as.data.frame(PM_top50_transformed@otu_table)

# combine otu and metadata - prepare for FactoMineR
facto_mat <- merge(otumat,sub_metadata, by = 0) # merge by row name
row.names(facto_mat) <- facto_mat$Row.names
facto_mat <- subset(facto_mat, select = -c(Row.names))
facto_mat <- na.omit(facto_mat)

### FactoMineR
PM_top50.pca <- PCA(facto_mat, quali.sup=585:587, graph = FALSE)
PM_top50.shiny=PCAshiny(PM_top50.pca)
