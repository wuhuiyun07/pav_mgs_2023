# https://hecatomb.readthedocs.io/en/latest/tutorialPt1/

# in R/Rstudio

# BiocManager::install(c("rstatix"))

library(tidyr)
library(dplyr)
library(ggplot2)
library(rstatix)
library(vegan)
getwd()
setwd("/Users/huiyunwu/Desktop/virus_particle/pav_mgs_2023/results/hecatomb")
data = read.csv('bigtable.tsv',header=T,sep='\t')
meta = read.csv('sample_metadata.csv', header=T)

taxonCounts = read.csv('taxonLevelCounts.tsv',header=T,sep='\t')

# save the merged tables to a new dataframe
# dataMeta = merge(data, meta, by='sampleID')

# If your metadata is incomplete and you want to keep all samples
dataMeta = merge(data, meta, by='sample_id', all=T)


# r$> colnames(dataMeta)
#  [1] "sampleID"            "seqID"               "count"               "percent"             "alnType"            
#  [6] "targetID"            "evalue"              "pident"              "fident"              "nident"             
# [11] "mismatches"          "qcov"                "tcov"                "qstart"              "qend"               
# [16] "qlen"                "tstart"              "tend"                "tlen"                "alnlen"             
# [21] "bits"                "targetName"          "taxMethod"           "kingdom"             "phylum"             
# [26] "class"               "order"               "family"              "genus"               "species"            
# [31] "baltimoreType"       "baltimoreGroup"      "host"                "sample_type"         "nucleotide_type"    
# [36] "collection_date"     "water_control"       "collection_location"


viruses = dataMeta %>% 
    filter(kingdom=="Viruses")

ggplot(viruses) + 
    geom_point(
        aes(x=alnlen,y=pident,color=alnType,size=count),
        alpha=0.1) + 
    facet_wrap(~family)
ggsave('plot/virus_1.png', width = 15, height = 6)

ggplot(viruses) +
    geom_point(
        aes(x=tlen,y=pident,color=alnType,size=count),
        alpha=0.1) +
    facet_wrap(~family) +
    geom_vline(xintercept=150,colour='red',linetype='longdash') +
    geom_hline(yintercept=75,colour='red',linetype='longdash')+
    scale_x_continuous(limits = c(0, 15000), 
                     breaks = seq(0, 15000, length.out = 4))
ggsave('plot/virus_2.png', width = 15, height = 6)


virusesFiltered = viruses %>% 
    filter(evalue<1e-10)

viruses = viruses %>% 
    mutate(filter=ifelse(evalue<1e-10,'pass','filter'))

ggplot(viruses) +
    geom_point(
        aes(x=alnlen,y=pident,color=filter),
        alpha=0.2) +
    facet_wrap(~family)
ggsave('plot/virus_3.png', width = 15, height = 6)


viruses = viruses %>% 
    mutate(filter=ifelse(alnlen>150 & pident>75,'pass','filter'))

ggplot(viruses) +
    geom_point(
        aes(x=alnlen,y=pident,color=filter),
        alpha=0.2) +
    facet_wrap(~family) +
    geom_vline(xintercept=150,colour='red',linetype='longdash') +
    geom_hline(yintercept=75,colour='red',linetype='longdash')
ggsave('plot/virus_4.png', width = 15, height = 6)

## Part 3: Visualising annotations

# get pmmov counts
pmmovCounts = taxonCounts %>% 
    filter(taxonLevel=='family',taxonName=='Virgaviridae')

# plot
ggplot(pmmovCounts) +
    geom_bar(aes(x=sampleID,y=count),stat='identity') +
    coord_flip()
ggsave('plot/pmmov.png', width = 12, height = 6)



# get enterovirus counts
enteroCounts = taxonCounts %>% 
    filter(taxonLevel=='family',taxonName=='Picornaviridae')

# plot
ggplot(enteroCounts) +
    geom_bar(aes(x=sampleID,y=count),stat='identity') +
    coord_flip()
ggsave('plot/enterovirus.png', width = 12, height = 6)


# get rotavirus counts
rotaCounts = taxonCounts %>% 
    filter(taxonLevel=='family',taxonName=='Reoviridae')

# plot
ggplot(rotaCounts) +
    geom_bar(aes(x=sampleID,y=count),stat='identity') +
    coord_flip()
ggsave('plot/rotavirus.png', width = 12, height = 6)


# get herpevirus counts
herpeCounts = taxonCounts %>% 
    filter(taxonLevel=='family',taxonName=='Hepeviridae')

# plot
ggplot(herpeCounts) +
    geom_bar(aes(x=sampleID,y=count),stat='identity') +
    coord_flip()
ggsave('plot/Hepeviridae.png', width = 12, height = 6)

# get all viral family counts
viralCounts = taxonCounts %>% 
    filter(taxonLevel=='family',grepl('k_Viruses',taxonPath))

# plot
ggplot(viralCounts) +
    geom_bar(aes(x=sampleID,y=count,fill=taxonName),position='stack',stat='identity') +
    coord_flip()
ggsave('plot/all.viral.family.png', width = 12, height = 6)


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
ggsave('plot/all.viral.family.filtered.png', width = 12, height = 6)

#visualising groups
virgaCounts = virusesFiltered %>% 
    group_by(family,sampleID,sample_type) %>% 
    filter(family=='Virgaviridae') %>% 
    summarise(n = sum(percent))
ggplot(virgaCounts) +
    geom_jitter(aes(x=sample_type,y=n),width = 0.1) +
    theme_bw()
ggsave('plot/virus_5.png', width = 12, height = 6)


viralFamCounts = virusesFiltered %>% 
    group_by(family) %>% 
    filter(family != 'unclassified Viruses family') %>%
    summarise(n=sum(percent)) %>% 
    arrange(desc(n))
ggplot(viralFamCounts) +
    geom_bar(aes(x=family,y=n),stat='identity') +
    coord_flip()
ggsave('plot/virus_6.png', width =12, height =6)

microCounts = virusesFiltered %>% 
    group_by(family,sampleID,sampleGroup) %>% 
    filter(family=='Microviridae') %>% 
    summarise(n = sum(percent))

ggplot(microCounts) +
    geom_jitter(aes(x=sampleGroup,y=n),width = 0.1) +
    theme_bw()

# collect counts
podoCounts = virusesFiltered %>% 
    group_by(family,sampleID,sampleGroup) %>% 
    filter(family=='Podoviridae') %>% 
    summarise(n = sum(percent))

# plot
ggplot(podoCounts) +
    geom_jitter(aes(x=sampleGroup,y=n),width = 0.1) +
    theme_bw()    


## Part4 Statistical tests

View(podoCounts)

t.test(
    podoCounts %>% 
        filter(sampleGroup=='A') %>% 
        pull(n),
    podoCounts %>% 
        filter(sampleGroup=='B') %>% 
        pull(n),
    alternative='two.sided',
    paired=F,
    var.equal=T)    

wilcox.test(
    podoCounts %>% 
        filter(sampleGroup=='A') %>% 
        pull(n),
    podoCounts %>% 
        filter(sampleGroup=='B') %>% 
        pull(n),
    alternative='t',
    paired=F)

# Dunn's test
# collect the family counts
viralFamCounts = virusesFiltered %>% 
    group_by(family) %>% 
    summarise(n=sum(percent)) %>% 
    arrange(desc(n))

# update factor levels to the sorted order
viralFamCounts$family = factor(viralFamCounts$family,levels=viralFamCounts$family)

ggplot(viralFamCounts) +
    geom_bar(aes(x=family,y=n),stat='identity') +
    coord_flip()


viralMajorFamCounts = viruses %>% 
    filter(family %in% c('Siphoviridae','Adenoviridae','Podoviridae','Microviridae')) %>% 
    group_by(sampleID,family,vaccine) %>% 
    summarise(n=sum(percent))

viralMajorFamCounts %>% 
    group_by(family) %>% 
    dunn_test(n ~ vaccine,p.adjust.method='holm') %>% 
    add_significance()

ggplot(viralMajorFamCounts) +
    geom_jitter(aes(x=vaccine,y=n)) +
    facet_wrap(~family)

# compare presence / absence
# apply a stringent filter
virusesStringent = viruses %>% 
    filter(evalue<1e-30,alnlen>150,pident>75,alnType=='aa')

# count hits for Myoviridae and score presence
# NOTE: the absent samples WONT be in the table yet, we have to add them in after
myovirPresAbs = virusesStringent %>% 
    filter(family=='Myoviridae') %>%
    group_by(sampleID) %>% 
    summarise(n=sum(percent)) %>%
    mutate(present=ifelse(n>0,1,0))

# merge in the metadata
# Using all=T will do an outer join and bring in the absent samples
myovirPresAbs = merge(myovirPresAbs,meta,by='sampleID',all=T)

# convert na's to zeros
myovirPresAbs[is.na(myovirPresAbs)] = 0

# matrix rows
mtxGroupA = c(
    myovirPresAbs %>% 
        filter(sampleGroup=='A',present==1) %>% 
        summarise(n=n()) %>% 
        pull(n),
    myovirPresAbs %>% 
        filter(sampleGroup=='A',present==0) %>% 
        summarise(n=n()) %>% 
        pull(n))
mtxGroupB = c(
    myovirPresAbs %>% 
        filter(sampleGroup=='B',present==1) %>% 
        summarise(n=n()) %>% 
        pull(n),
    myovirPresAbs %>% 
        filter(sampleGroup=='B',present==0) %>% 
        summarise(n=n()) %>% 
        pull(n))

# create the 2x2 matrix
myovirFishMtx = matrix(c(mtxGroupA,mtxGroupB),nrow = 2)

# this bit is not necessary, but lets add row and col names to illustrate the matrix layout
colnames(myovirFishMtx) = c('GroupA','GroupB')
row.names(myovirFishMtx) = c('present','absent')

# view
View(myovirFishMtx)

## contig  annotations
contigTable = read.csv('contigSeqTable.tsv.gz',header=T,sep='\t')

View(contigTable)

ggplot(contigTable %>% 
        filter(contigID=='contig_6')
    ) +
    geom_point(aes(x=start,y=genus,color=family))

ggplot(contigTable %>% 
        filter(contigID=='contig_6')
    ) +
    geom_point(aes(x=start,y=species,color=family))

ggplot(contigTable %>% 
           filter(contigID=='contig_56')
) +
    geom_point(aes(x=start,y=family,color=kingdom))    


ggplot(contigTable %>% 
           filter(contigID=='contig_56')
) +
    geom_point(aes(x=start,y=species,color=kingdom))   

## Part 6: Viral ecology
library(stringr)
library(tidyverse)
data = read.csv('bigtable.tsv', sep = '\t',header=T)
data <- data %>%
  mutate(sampleID = str_replace_all(sampleID, pattern = "_L001", replacement = "")) %>% 
  filter(tlen>1000)

meta = read_csv('sample_metadata.csv')
data = merge(data, meta, by = 'sampleID')

familyCounts = data %>%
  group_by(sampleID,filter_type,collection_location,family) %>%
  summarise(n=sum(percent)) %>%
  spread(family,n,fill=0)

View(familyCounts)


familyMtx = as.matrix(familyCounts[,4:ncol(familyCounts)])

View(familyMtx)

#PERMANOVA
adonis2(familyMtx ~ filter_type, data = familyCounts, permutations = 1000, method="bray", by=NULL)

# sample type
adonis2(familyMtx ~ collection_location, data = familyCounts, permutations = 1000, method="bray",by=NULL)


# position
familySourceSimper = with(familyCounts, simper(familyMtx, collection_location, permutations = 1000))
familyPositionSum = summary(familySourceSimper, ordered = T)
View(familyPositionSum$inner_outer)

# position x source
familyPositionSimper = with(familyCounts, simper(familyMtx, positionSource, permutations = 1000))
speciesPosSource = summary(speciesPositionSimper, ordered = T)
View(speciesPosSource$`outer coral mucus_outer reef water`)

# pull the data frame
outerReef = as.data.frame(speciesPosSource$`outer coral mucus_outer reef water`)

# filter for significant cummulative sums totalling 90% of variance
outerReef = outerReef %>% 
  filter(p<=0.05, cumsum < 0.9)

# plot
ggplot(outerReef) +
  geom_bar(aes(x=row.names(outerReef), y=ava, fill='a'), stat='identity') +
  geom_bar(aes(x=row.names(outerReef), y=-avb, fill='b'), stat='identity') +
  coord_flip() +
  theme_bw()

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

ggplot(familyScores, aes(x=NMDS1,y=NMDS2, color=filter_type, fill=filter_type))+
  stat_ellipse(size=0.75)+
  geom_point(shape=21, color='black', size=3)+
  theme_bw()


# PCoA
speciesDist = vegdist(speciesMtx, method = 'bray')

View(speciesDist)

speciesDispersion = betadisper(speciesDist,group = speciesScores$positionSource)

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


