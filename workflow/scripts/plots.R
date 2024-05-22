

library(tidyr)
library(dplyr)
library(ggplot2)
library(rstatix)
library(tidyverse)



screenedFiles <- list.files("results/annotation2", full.name = TRUE)

screened<-read_csv(screenedFiles)
# write_csv(screened, "results/screenedFiles.csv")
colnames(screened)
#     colnames(screened)
#  [1] "sampleID"            "contig_id"           "gene.x"              "dsDNAphage"          "ssDNA"               "NCLDV"               "RNA"                
#  [8] "lavidaviridae"       "max_score"           "max_score_group"     "length.x"            "hallmark"            "viral"               "cellular"           
# [15] "contig_length"       "provirus"            "proviral_length"     "gene_count"          "viral_genes"         "host_genes"          "checkv_quality"     
# [22] "miuvig_quality"      "completeness"        "completeness_method" "contamination"       "kmer_freq"           "warnings"            "gene.y"             
# [29] "sseqid"              "pident"              "length.y"            "mismatch"            "evalue"              "bitscore"            "staxids"            
# [36] "sscinames"           "sskingdoms"          "skingdoms"           "sphylums"            "stitle"             


data1<- read_csv("results/bigtable2.csv")
# colnames(data1)
# [1] "staxids"      "superkingdom" "phylum"       "class"        "order"        "family"       "genus"       
# [8] "species"     
data2<-read_csv("results/screenedFiles.csv")
colnames(data2)
#  [1] "sampleID"            "contig_id"           "gene.x"              "dsDNAphage"          "ssDNA"               "NCLDV"               "RNA"                
#  [8] "lavidaviridae"       "max_score"           "max_score_group"     "length.x"            "hallmark"            "viral"               "cellular"           
# [15] "contig_length"       "provirus"            "proviral_length"     "gene_count"          "viral_genes"         "host_genes"          "checkv_quality"     
# [22] "miuvig_quality"      "completeness"        "completeness_method" "contamination"       "kmer_freq"           "warnings"            "gene.y"             
# [29] "sseqid"              "pident"              "length.y"            "mismatch"            "evalue"              "bitscore"            "staxids"            
# [36] "sscinames"           "sskingdoms"          "skingdoms"           "sphylums"            "stitle"       
# data<- full_join(data1, data3, by = "staxids")
data<- left_join(data2, data1, by = 'staxids')%>%
distinct(contig_id, .keep_all = TRUE)
write_csv(data, "results/bigdata.disctict.csv")

data <- read_csv("results/bigdata.disctict.csv")
colnames(data)
#  [1] "sampleID"            "contig_id"           "gene.x"              "dsDNAphage"          "ssDNA"               "NCLDV"               "RNA"                
#  [8] "lavidaviridae"       "max_score"           "max_score_group"     "length.x"            "hallmark"            "viral"               "cellular"           
# [15] "contig_length"       "provirus"            "proviral_length"     "gene_count"          "viral_genes"         "host_genes"          "checkv_quality"     
# [22] "miuvig_quality"      "completeness"        "completeness_method" "contamination"       "kmer_freq"           "warnings"            "gene.y"             
# [29] "sseqid"              "pident"              "length.y"            "mismatch"            "evalue"              "bitscore"            "staxids"            
# [36] "sscinames"           "sskingdoms"          "skingdoms"           "sphylums"            "stitle"              "superkingdom"        "phylum"             
# [43] "class"               "order"               "family"              "genus"               "species"            

viruses = data 

high_quality <- viruses %>%
    arrange(desc(pident))%>%
    filter(checkv_quality == "High-quality")%>%
    select(sampleID, contig_id, max_score, max_score_group, contig_length, sseqid, staxids, checkv_quality, length.x, length.y, completeness, pident, family, species, stitle)

write_csv(high_quality, "results/high_quality.csv")


# r$> colnames(viruses)
# [1] "staxids"             "superkingdom"        "phylum"              "class"               "order"              
# [6] "family"              "genus"               "species"             "contig_id"           "gene.x"             
# [11] "dsDNAphage"          "ssDNA"               "NCLDV"               "RNA"                 "lavidaviridae"      
# [16] "max_score"           "max_score_group"     "length.x"            "hallmark"            "viral"              
# [21] "cellular"            "contig_length"       "provirus"            "proviral_length"     "gene_count"         
# [26] "viral_genes"         "host_genes"          "checkv_quality"      "miuvig_quality"      "completeness"       
# [31] "completeness_method" "contamination"       "kmer_freq"           "warnings"            "gene.y"             
# [36] "sseqid"              "pident"              "length.y"            "mismatch"            "evalue"             
# [41] "bitscore"            "sscinames"           "sskingdoms"          "skingdoms"           "sphylums"           
# [46] "stitle"     

ggplot(viruses) + 
    geom_point(
        aes(x=contig_length,y=pident),
        alpha=0.2) + 
    facet_wrap(~family)
ggsave('results/plots/virus_1.png', width = 15, height = 6)

ggplot(viruses) + 
    geom_point(
        aes(x=contig_length,y=pident, size=max_score),
        alpha=0.1) + 
    facet_wrap(~max_score_group)
ggsave('results/plots/virus_max_score_group.png', width = 15, height = 6)

ggplot(viruses) + 
    geom_point(
        aes(x=contig_length,y=viral, color=checkv_quality, size=max_score),
        alpha=0.1) + 
    facet_wrap(~max_score_group)
ggsave('results/plots/viral.png', width = 15, height = 6)

ggplot(viruses) +
    geom_point(
        aes(x=contig_length,y=pident,color=checkv_quality,size=completeness),
        alpha=0.1) +
    facet_wrap(~family) +
    geom_vline(xintercept=150,colour='red',linetype='longdash') +
    geom_hline(yintercept=75,colour='red',linetype='longdash')
ggsave('results/plots/virus_2.png', width = 12, height = 6)

ggplot(viruses) +
    geom_point(
        aes(x=contig_length,y=pident,color=max_score_group,size=max_score),
        alpha=0.1) +
    facet_wrap(~family) +
    geom_vline(xintercept=150,colour='red',linetype='longdash') +
    geom_hline(yintercept=75,colour='red',linetype='longdash')
    # scale_x_continuous(limits = c(0, 15000), 
    #                  breaks = seq(0, 15000, length.out = 4))
ggsave('results/plots/max_score_group.png', width = 15, height = 8)




ggplot(viruses) +
    geom_point(
        aes(x=contig_length,y=pident,color=checkv_quality,size=viral_genes),
        alpha=0.1) +
    facet_wrap(~family) +
    geom_vline(xintercept=150,colour='red',linetype='longdash') +
    geom_hline(yintercept=75,colour='red',linetype='longdash')
ggsave('results/plots/viral_genes.png', width = 15, height = 6)


virusesFiltered = viruses %>% 
    filter(evalue<1e-10)

viruses = viruses %>% 
    mutate(filter=ifelse(evalue<1e-10,'pass','filter'))

ggplot(viruses) +
    geom_point(
        aes(x=contig_length,y=pident,color=filter),
        alpha=0.2) +
    facet_wrap(~family)
ggsave('results/plots/virus_3.png', width = 12, height = 6)


viruses = viruses %>% 
    mutate(filter=ifelse(contig_length>150 & pident>75,'pass','filter'))

ggplot(viruses) +
    geom_point(
        aes(x=contig_length,y=pident,color=filter),
        alpha=0.2) +
    facet_wrap(~family) +
    geom_vline(xintercept=150,colour='red',linetype='longdash') +    
    geom_hline(yintercept=75,colour='red',linetype='longdash')
ggsave('results/plots/virus_4.png', width = 15, height = 6)

# norovirus family: Caliciviridae
# enterovirus: Picornavirus
# rotavirus: Reoviridae
# pmmov: Virgaviridae
# crAssphage: Intestiviridae


# get pmmov counts
pmmovCounts = taxonCounts %>% 
    filter(taxonLevel=='family',taxonName=='Virgaviridae')

# plot
ggplot(pmmovCounts) +
    geom_bar(aes(x=sampleID,y=count),stat='identity') +
    coord_flip()
ggsave('results/plots/pmmov.png', width = 12, height = 6)



# get norovirus counts
noroCounts = taxonCounts %>% 
    filter(taxonLevel=='family',taxonName=='Caliciviridae')

# plot
ggplot(noroCountsCounts) +
    geom_bar(aes(x=sampleID,y=count),stat='identity') +
    coord_flip()
ggsave('results/plots/norovirus.png', width = 12, height = 6)


# get all viral family counts
viralCounts = taxonCounts %>% 
    filter(taxonLevel=='family',grepl('k_Viruses',taxonPath))

# plot
ggplot(viralCounts) +
    geom_bar(aes(x=sampleID,y=count,fill=taxonName),position='stack',stat='identity') +
    coord_flip()
ggsave('results/plots/all.viral.family.png', width = 12, height = 6)


# Answer for "Challenge: Filter your raw viral hits to only keep protein hits with an evalue < 1e-10"
virusesFiltered = viruses %>% 
    filter(alnType=='aa',evalue<1e-10)

# collect the filtered counts
viralFiltCounts = virusesFiltered %>% 
    group_by(sampleID,family) %>% 
    summarise(n = sum(percent))

ggplot(viralFiltCounts) +
    geom_bar(aes(x=sampleID,y=n,fill=family),position='fill',stat='identity') +
    coord_flip() +
    theme_bw()
ggsave('results/plots/all.viral.family.filtered.png', width = 12, height = 6)

#visualising groups
virgaCounts = virusesFiltered %>% 
    group_by(family,sampleID,sample_type) %>% 
    filter(family=='Virgaviridae') %>% 
    summarise(n = sum(percent))
ggplot(virgaCounts) +
    geom_jitter(aes(x=sample_type,y=n),width = 0.1) +
    theme_bw()
ggsave('results/plots/virus_5.png', width = 12, height = 6)


viralFamCounts = virusesFiltered %>% 
    group_by(family) %>% 
    summarise(n=sum(percent)) %>% 
    arrange(desc(n))
ggplot(viralFamCounts) +
    geom_bar(aes(x=family,y=n),stat='identity') +
    coord_flip()
ggsave('results/plots/virus_6.png', width =12, height =6)


## Part 6: Viral ecology
library(stringr)
library(tidyverse)
library(vegan)
data = read_csv('results/bigdata.disctict.csv')


meta = read_csv('results/sample_metadata.csv')
data = merge(data, meta, by = 'sampleID')
# write_csv(data, "results/bigdataFull.csv")


speciesCounts = data %>%
  group_by(sampleID,filter_type,collection_location,species) %>% 
  summarise(species_count = n()) %>% 
  # left_join(data, by = c("sampleID","filter_type","collection_location","species")) %>% 
  spread(key = species, value = species_count, fill=0 )


familyCounts = data %>%
  group_by(sampleID,filter_type,collection_location,family) %>% 
  summarise(family_count = n()) %>% 
  spread(key = family, value = family_count, fill=0 )

orderCounts = data %>%
  group_by(sampleID,filter_type,collection_location,order) %>% 
  summarise(order_count = n()) %>% 
  # left_join(data, by = c("sampleID","filter_type","collection_location","species")) %>% 
  spread(key = order, value = order_count, fill=0 )


speciesMtx = as.matrix(speciesCounts[,4:ncol(speciesCounts)])
familyMtx = as.matrix(familyCounts[,4:ncol(familyCounts)])
orderMtx = as.matrix(orderCounts[,4:ncol(orderCounts)])



#PERMANOVA
adonis2(speciesMtx ~ filter_type, data = speciesCounts, permutations =1000, method = "bray", by=NULL)
adonis2(familyMtx ~ filter_type, data = familyCounts, permutations = 1000, method="bray", by=NULL)
# Df SumOfSqs      R2     F Pr(>F)
# Model     4   1.7830 0.19327 1.138 0.2498
# Residual 19   7.4424 0.80673             
# Total    23   9.2254 1.00000       
adonis2(orderMtx ~ filter_type, data = orderCounts, permutations = 1000, method="bray", by=NULL)
# filter type
adonis2(familyMtx ~ collection_location, data = familyCounts, permutations = 1000, method="bray",by=NULL)
# adonis2(formula = familyMtx ~ collection_location, data = familyCounts, permutations = 1000, method = "bray", by = NULL)
# Df SumOfSqs      R2      F  Pr(>F)  
# Model     3   1.8843 0.20425 1.7111 0.01399 *
#   Residual 20   7.3412 0.79575                 
# Total    23   9.2254 1.00000    
adonis2(speciesMtx ~ collection_location, data = familyCounts, permutations = 1000, method="bray",by=NULL)

# position
familySourceSimper = with(familyCounts, simper(familyMtx, collection_location, permutations = 1000))
familyPositionSum = summary(familySourceSimper, ordered = T)
View(familyPositionSum$inner_outer)




#NMDS
 # calculate the NMDS from the count matrix
familyNmds = metaMDS(familyMtx, distance = 'bray', k = 5, try = 20, trymax = 10000)

# pull the dataframe for plotting
familyScores = as.data.frame(scores(familyNmds)$site)

# add the metada back in (speciesMtx contains the counts only)
familyScores$collection_location = familyCounts$collection_location
familyScores$filter_type = familyCounts$filter_type

# plot
ggplot(familyScores, aes(x=NMDS1,y=NMDS2, color=collection_location, fill=collection_location))+
  stat_ellipse(size=0.75)+
  geom_point(shape=21, color='black', size=3)+
  theme_bw()
ggsave("results/plots/nmds_on_location.png")

ggplot(familyScores, aes(x=NMDS1,y=NMDS2, color=filter_type, fill=filter_type))+
  stat_ellipse(size=0.75)+
  geom_point(shape=21, color='black', size=3)+
  theme_bw()
ggsave("results/plots/nmds_on_particle_size.png")

ggplot(speciesScores, aes(x=NMDS1,y=NMDS2, color=filter_type, fill=filter_type))+
  stat_ellipse(size=0.75)+
  geom_point(shape=21, color='black', size=3)+
  theme_bw()
ggsave("results/plots/nmds_on_particle_size.png")


# PCoA
speciesDist = vegdist(speciesMtx, method = 'bray')

View(speciesDist)

speciesDispersion = betadisper(speciesDist,group = speciesScores$filter_type)

permutest(speciesDispersion)

# pull the coordinates
dispVec = as.data.frame(speciesDispersion$vectors)

# add the sample metadata back in
dispVec$positionSource = speciesScores$positionSource

# reorder the samples - not necessary but for these samples it works better for the paired color palette
dispVec$positionSource = factor(
  dispVec$positionSource,
  levels=c('inner reef water',
           'outer reef water',
           'inner coral mucus', 
           'outer coral mucus'))

# plot!
ggplot(dispVec,aes(x=PCoA1,y=PCoA2,color=positionSource,fill=positionSource))+
  stat_ellipse(size=0.75, show.legend = F)+
  geom_point(shape=21,size=3,color='black')+
  scale_color_brewer(palette = 'Paired') +
  scale_fill_brewer(palette = 'Paired')+
  theme_bw()+
  guides(fill=guide_legend(title="Sample location & source"))