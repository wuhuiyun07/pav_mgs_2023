# Load necessary libraries
library(vegan)
library(ggforce)
library(ggplot2)
library(readxl)
library(tidyverse)
library(cluster)
library(dplyr)
library(tidyr)
library(factoextra)
library(phyloseq)
library(openxlsx)
library(reshape2)

# BiocManager::install(c("phyloseq"))
# BiocManager::install(c("ggforce"))
# BiocManager::install(c("readxl"))
# BiocManager::install(c("cluster"))
# BiocManager::install(c("factoextra"))
# BiocManager::install(c("openxlsx"))
# BiocManager::install(c("reshape2"))
# install.packages("cli")


##tutorial https://rpubs.com/lconteville/714853
##tutorial https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
getwd()
# Read the Excel file
bivalve_data <- read_excel("results/onr/bivalve_data.xlsx")

# Load sample data
sample_metadata <- read.xlsx("results/onr/bivalve_sample.xlsx")


# Identify unique species and assign OTU names
unique_species <- unique(bivalve_data$species)
otu_codes <- paste0("Otu", sprintf("%05d", 1:length(unique_species)))
species_to_otu <- data.frame(species = unique_species, OTU = otu_codes)




# Merge OTU codes with the bivalve data
bivalve_data <- bivalve_data %>%
  left_join(species_to_otu, by = "species")

# Extract abundance data for the OTU table
otu_data <- bivalve_data %>%
  group_by(sample.name, OTU) %>%
  summarise(abundance = n()) %>%
  pivot_wider(names_from = sample.name, values_from = abundance, values_fill = 0)


# Convert otu_data tibble to data frame
otu_table <- as.data.frame(otu_data)

# Extract taxonomic information for taxonomic table
tax_table <- bivalve_data %>%
  select(superkingdom, phylum, class, order, family, genus, species, OTU) %>%
  distinct()

# Convert tax_table tibble to data frame
tax_table <- as.data.frame(tax_table)


# Save OTU table as CSV
# write.csv(otu_mat, file = "C:/Users/katie/OneDrive/Desktop/R/ONR_bivalves/data/otu_table.csv", row.names = TRUE)

# Save taxonomic table as CSV
# write.csv(tax_mat, file = "C:/Users/katie/OneDrive/Desktop/R/ONR_bivalves/data/tax_table.csv", row.names = TRUE)



# ## Read in new OTU, tax, and sample tables
# otu_mat <- read.csv("C:/Users/katie/OneDrive/Desktop/R/ONR_bivalves/data/otu_table.csv")
# tax_mat <- read.csv("C:/Users/katie/OneDrive/Desktop/R/ONR_bivalves/data/tax_table.csv")
# samples_df <- read_excel("C:/Users/katie/OneDrive/Desktop/R/ONR_bivalves/data/bivalve_sample_1.xlsx")

otu_mat <- otu_table
tax_mat <- tax_table
samples_df<- sample_metadata

str(samples_df)
colnames(otu_mat)

## Define the row names from the OTU column
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("OTU")

## Remove duplicate rows from the taxonomic table
tax_mat <- unique(tax_mat)

# Convert taxonomic table to a matrix
tax_mat <- as.matrix(tax_mat)

# Assign OTU as row names
rownames(tax_mat) <- tax_mat[, "OTU"]

## Define row names for the sample table
samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample.name")

# Transform into matrices otu and tax tables (sample table can be left as a data frame)
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

## Transform to phyloseq objects
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(samples_df)

# Identify duplicated taxa names
duplicated_taxa <- duplicated(taxa_names(TAX))

# Check if there are any duplicated taxa names
if (any(duplicated_taxa)) {
  # Remove duplicated taxa from taxonomic table
  TAX <- TAX[!duplicated_taxa, ]
}

# Create phyloseq object
carbom <- phyloseq(OTU, TAX, samples)
carbom


##visualize data
sample_names(carbom)

rank_names(carbom)

sample_variables(carbom)

# Load required library

# Subset to marine mammal families
carbom_subset <- subset_taxa(carbom, family %in% c("Adenoviridae", "Anelloviridae", "Caliciviridae", 
                                                   "Coronaviridae", "Hepadnaviridae", "Orthomyxoviridae",
                                                   "Papillomaviridae", "Paramyxoviridae", "Parvoviridae",
                                                   "Picornaviridae", "Retroviridae", "Rhabdoviridae"))

# Convert counts to relative abundance
carbom_subset_rel <- transform_sample_counts(carbom_subset, function(x) 100 * x / sum(x))


# Plot bar chart with customized legend title and y-axis label
bivalve_barplot <- plot_bar(carbom_subset_rel, fill = "family") + 
  geom_bar(aes(color = family, fill = family), stat = "identity", position = "stack") +
  # Customize legend title
  labs(fill = "Viral Family") +
  # Customize y-axis label
  ylab("% Abundance") +
  # Adjust legend appearance
  theme(legend.key = element_rect(colour = NA)) +
  guides(color = FALSE)

# Save the plot as a high-resolution figure
ggsave("results/onr/bivalve_barplot.png", bivalve_barplot, dpi = 300, width = 8, height = 6, units = "in")




#########Stacked barplot with different colors###########

# Subset to marine mammal families
carbom_subset <- subset_taxa(carbom, family %in% c("Adenoviridae", "Anelloviridae", "Caliciviridae", 
                                                   "Coronaviridae", "Hepadnaviridae", "Orthomyxoviridae",
                                                   "Papillomaviridae", "Paramyxoviridae", "Parvoviridae",
                                                   "Picornaviridae", "Retroviridae", "Rhabdoviridae"))

# Convert counts to relative abundance
carbom_subset_rel <- transform_sample_counts(carbom_subset, function(x) 100 * x / sum(x))

# Define custom color palette
custom_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", 
                    "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ffbb78", "#98df8a")

# Plot bar chart with customized color palette
bivalve_barplot_colors <- plot_bar(carbom_subset_rel, fill = "family") + 
  geom_bar(aes(color = family, fill = family), stat = "identity", position = "stack") +
  # Use custom color palette
  scale_fill_manual(values = custom_palette) +
  scale_color_manual(values = custom_palette) +
  # Customize legend title
  labs(fill = "Viral Family") +
  # Customize y-axis label
  ylab("% Abundance") +
  # Adjust legend appearance
  theme(legend.key = element_rect(colour = NA)) +
  guides(color = FALSE) 
  # Save the plot as a high-resolution figure
  ggsave("results/onr/bivalve_barplot_colors.png", bivalve_barplot_colors, dpi = 300, width = 10, height = 6, units = "in")


##add sample date to the top of the graph
  
  # Subset to marine mammal families
  carbom_subset <- subset_taxa(carbom, family %in% c("Adenoviridae", "Anelloviridae", "Caliciviridae", 
                                                     "Coronaviridae", "Hepadnaviridae", "Orthomyxoviridae",
                                                     "Papillomaviridae", "Paramyxoviridae", "Parvoviridae",
                                                     "Picornaviridae", "Retroviridae", "Rhabdoviridae"))
  
  # Convert counts to relative abundance
  carbom_subset_rel <- transform_sample_counts(carbom_subset, function(x) 100 * x / sum(x))
  
  # Define custom color palette
  custom_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", 
                      "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#ffbb78", "#98df8a")
  
  # Plot bar chart with customized color palette and facet by collection date
  bivalve_barplot_colors <- plot_bar(carbom_subset_rel, fill = "family") + 
    geom_bar(aes(color = family, fill = family), stat = "identity", position = "stack") +
    # Use custom color palette
    scale_fill_manual(values = custom_palette) +
    scale_color_manual(values = custom_palette) +
    # Customize legend title
    labs(fill = "Viral Family") +
    # Customize y-axis label
    ylab("% Abundance") +
    # Adjust legend appearance
    theme(legend.key = element_rect(colour = NA)) +
    guides(color = FALSE) +
    # Facet by collection date
    facet_grid(. ~ collection.date, scales = "free_x", space = "free_x")
  
  # Save the plot as a high-resolution figure
  ggsave("results/onr/bivalve_barplot_colors.png", bivalve_barplot_colors, dpi = 300, width = 10, height = 6, units = "in")
  

##########HEATMAP#######
  # Load necessary library
  # Load necessary library
  # Load necessary library
  library(phyloseq)
  
  # Subset to marine mammal families
  carbom_subset <- subset_taxa(carbom, family %in% c("Adenoviridae", "Anelloviridae", "Caliciviridae", 
                                                     "Coronaviridae", "Hepadnaviridae", "Orthomyxoviridae",
                                                     "Papillomaviridae", "Paramyxoviridae", "Parvoviridae",
                                                     "Picornaviridae", "Retroviridae", "Rhabdoviridae"))
  
  # Create a list to store aggregated OTU tables
  aggregate_otu <- list()
  
  # Aggregate OTU counts within each viral family
  for (family in viral_families) {
    subset_family <- subset_taxa(carbom_subset, family == family)
    if (nrow(otu_table(subset_family)) > 0) {
      aggregate_otu[[family]] <- tax_glom(subset_family, taxrank = "family")
    }
  }
  
  # Combine aggregated OTU tables
  combined_otu <- do.call(combine, lapply(aggregate_otu, otu_table))
  
  # Append unique identifiers to taxon names
  taxa_names(combined_otu) <- make_unique(taxa_names(combined_otu))
  
  # Plot heatmap
  plot_heatmap(phyloseq(otu_table(combined_otu)), method = "NMDS", distance = "bray", 
               sample_label = "none", taxa_label = "row")
  
 

  ######beta-diversity plots########
 
  # Calculate abundance of each viral species
  viral_abundance <- tax_table(carbom)[, "species"]
  
  # Subset only the viral species
  viral_species <- subset_taxa(carbom, species != "NA")
  
  # Calculate abundance of each viral species
  viral_abundance <- tax_table(viral_species)[, "species"]
  
  # Calculate beta diversity
  beta_diversity <- vegdist(otu_table(viral_species), method = "bray")
  
  # Perform NMDS
  nmds <- metaMDS(beta_diversity)
  
  # Extract NMDS coordinates
  nmds_data <- data.frame(nmds$points)
  
  # Rename the columns
  colnames(nmds_data) <- c("NMDS1", "NMDS2")
  
  # Assuming carbom is your phyloseq object
  
  # Access the sample metadata
  sample_metadata <- sample_data(carbom)
  
  # Change sample site names from "p160" and "p161" to "p159"
  sample_metadata$sample.site.1[sample_metadata$sample.site.1 %in% c("p160", "p161")] <- "p159"
  
  # Update the sample metadata in the phyloseq object
  sample_data(carbom) <- sample_metadata
  
  
  # Extract sample site information from carbom dataset
  sample_sites <- sample_data(carbom)$sample.site.1
  
  # Repeat sample site information to match the number of rows in nmds_data
  sample_sites <- rep(sample_sites, length.out = nrow(nmds_data))
  
  # Merge sample site information with NMDS data
  nmds_data_with_sites <- cbind(nmds_data, sample_site = sample_sites)
  
  # Assuming you've already performed NMDS and have the nmds object
  
  # Plot NMDS with sample sites
  ggplot(nmds_data_with_sites, aes(x = NMDS1, y = NMDS2, color = sample_site)) +
    geom_point() +
    theme_minimal()
  
  
  ##bray curtis boxplot

  # Ensure sample_sites matches the length of beta_diversity
  sample_sites <- sample_sites[1:length(beta_diversity)]
  
  # Create a data frame with beta diversity and sample site information
  beta_data <- data.frame(Bray_Curtis = as.vector(beta_diversity), Sample_Site = sample_sites)
  
  # Remove NA values from beta_data
  beta_data <- na.omit(beta_data)
  
  # Plot boxplot
  ggplot(beta_data, aes(x = Sample_Site, y = Bray_Curtis)) +
    geom_boxplot(fill = "skyblue", color = "black") +
    labs(x = "Sample Site", y = "Bray-Curtis Dissimilarity") +
    theme_minimal()
  
 
  ###PCoA Principal Coordiante Analysis or Multidimensional Scaling (MDS) used to visualize and analyzing similaries
  ##or dissimilarities among a set of samples
  # Assuming you have calculated Bray-Curtis dissimilarities and stored them in beta_data
  
  # Perform PCoA
  pcoa <- cmdscale(as.dist(beta_data$Bray_Curtis), k = 2)
  
  # Create a data frame with PCoA coordinates
  pcoa_data <- data.frame(PCoA1 = pcoa[, 1], PCoA2 = pcoa[, 2], Sample_Site = beta_data$Sample_Site)
  
  # Plot PCoA
pcoa_plot<-ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, color = Sample_Site, shape = Sample_Site)) +
    geom_point(size = 1) +
    stat_ellipse(aes(group = Sample_Site)) +
    labs(color = "Sample Site", shape = "Sample Site")
  
  # Save as PNG with high resolution (300 dpi)
  ggsave("C:/Users/katie/OneDrive/Desktop/R/ONR_bivalves/results/pcoa_plot.png", pcoa_plot, width = 6, height = 4, dpi = 300)
  
  
  
  
  
  
  
  
  






