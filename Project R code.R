knitr::opts_chunk$set(echo = TRUE)

knitr::opts_knit$set(root.dir = '~/Desktop/Microbiota_analysis_R/')

#Analyses of Phylogenetics and Evolution package. Required for tree calculations to be used with phyloseq
#install.packages("ape")
#install.packages("phangorn")
#install.packages("VennDiagram")

library(ape)
#This package will also help us more easily manipulate our data
library(dplyr)
#Graphing package used in phyloseq. To edit the default setting of a plot, you need to use functions in this package.
library(ggplot2)
#This package is used to calculate and plot Venn diagrams as well as heatmaps
library(gplots)
#Linear mixed-effects models like repeated measures analysis
library(lme4)
#used to read in mothur-formatted files
library(phangorn)
#The phyloseq package seeks to address issues with multiple microbiome analysis packages by providing a set of functions that internally manage the organizing, linking, storing, and analyzing of phylogenetic sequencing data. In general, this package is used for UniFrac analyses.
library(phyloseq)
#A package to create interactive web graphics of use in 3D plots
library(plotly)
#This package will help us more easily manipulate our data, which are matrices
library(tidyr)
#The vegan package provides tools for descriptive community ecology. It has most basic functions of diversity analysis, community ordination and dissimilarity analysis. In general, this package is used for Bray-Curtis and Jaccard analyses.
library(vegan)
#Pretty Venn disgrams
library(VennDiagram)
library(tibble)
#Loading Data
OTU= read.table("/Users/clapp5/Desktop/Microbiota_analysis_R/project/rarified-table.tsv", header=TRUE, sep="\t", stringsAsFactors = F)
#Taxonomy of each OTU
tax = read.table("/Users/clapp5/Desktop/Microbiota_analysis_R/project/taxonomy.tsv", header=TRUE, sep="\t")
#Metadata. Since we made this in Excel, not mothur, we can use the "row.names" modifier to automatically name the rows by the values in the first column (sample names)
meta = read.table("/Users/clapp5/Desktop/Microbiota_analysis_R/project/sample_info.tsv", header=TRUE, row.names=1, sep="\t")
evenness = read.table("/Users/clapp5/Desktop/Microbiota_analysis_R/project/evenness.tsv", header=TRUE, row.names=1, sep="\t")
faith_pd = read.table("/Users/clapp5/Desktop/Microbiota_analysis_R/project/faith_pd.tsv", header=TRUE, row.names=1, sep="\t")
observed_otus = read.table("/Users/clapp5/Desktop/Microbiota_analysis_R/project/observed_otus.tsv", header=TRUE, row.names=1, sep="\t")
shannon = read.table("/Users/clapp5/Desktop/Microbiota_analysis_R/project/shannon.tsv", header=TRUE, row.names=1, sep="\t")

#Transpose the OTU table
row.names(OTU) = OTU[,1]
OTU[,1] = NULL
OTU.clean = as.data.frame(t(as.matrix(OTU)))#OTU[,-1]

str(OTU.clean)

#Rename Rows by OTUS in Taxonomy table
tax2 = separate(tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#All the strings that need to be removed and replaced with NA
na_strings <- c(" s__", " g__", " f__", " o__", " c__")
#str(tax.clean)
tax3 = replace_with_na_all(tax2, condition = ~.x %in% na_strings)
tax3[] <- t(apply(tax3, 1, zoo::na.locf))
tax3 = tax3[,-1]
tax3_mat <- as.matrix(tax3)

row.names(tax3_mat) <- tax2[,1]
tax.clean <- as.data.frame(tax3_mat)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.
tax.clean = tax.clean[row.names(tax.clean) %in% colnames(OTU.clean),]

#Merge Diversity with Metadata table
alpha_diversity = merge(faith_pd, evenness, by.x = 0, by.y = 0)
alpha_diversity = merge(alpha_diversity, observed_otus, by.x = "Row.names", by.y = 0)
alpha_diversity = merge(alpha_diversity, shannon, by.x = "Row.names", by.y = 0)
meta = merge(meta, alpha_diversity, by.x = 0, by.y = "Row.names")
row.names(meta) = meta$Row.names
meta = meta[,-1]

#Order the Data
OTU.clean = OTU.clean[order(row.names(OTU.clean)),]
meta = meta[order(row.names(meta)),]

##Set seedWe will be running some processes that rely on the random number generater. To make your analysis reproducible, we set the random seed.

set.seed(8765)
#Create 2x2 plot environment 
par(mfrow = c(2, 2))

#Plots
hist(meta$shannon, main="Shannon diversity", xlab="", breaks=10)
hist(meta$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(meta$pielou_e, main="Evenness", xlab="", breaks=15)
hist(meta$observed_otus, main="Observed OTUs", xlab="", breaks=15)

#To test for normalcy statistically, we can run the Shapiro-Wilk test of normality.
shapiro.test(meta$shannon)
shapiro.test(meta$pielou_e)
shapiro.test(meta$faith_pd)
shapiro.test(meta$observed_otus)

#Run the ANOVA and save it as an object
aov.observed_otus.treatment = aov(observed_otus ~ treatment, data=meta)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.observed_otus.treatment)

aov.evenness.sex = aov(pielou_e ~ sex, data=meta)
summary(aov.evenness.sex)

aov.pielou_e.treatment = aov(pielou_e ~ treatment, data=meta)
summary(aov.pielou_e.treatment)

aov.shannon.treatment = aov(shannon ~ treatment, data=meta)
summary(aov.shannon.treatment)

aov.shannon.sex = aov(shannon ~ sex, data=meta)
summary(aov.shannon.sex)

aov.faith_pd.treatment = aov(faith_pd ~ treatment, data = meta)
summary(aov.faith_pd.treatment)

#Pairwise Comparison
TukeyHSD(aov.evenness.sex)
TukeyHSD(aov.shannon.sex)
TukeyHSD(aov.observed_otus.treatment)
TukeyHSD(aov.shannon.treatment)
TukeyHSD(aov.pielou_e.treatment)
TukeyHSD(aov.faith_pd.treatment)

#Plot Data
levels(meta$treatment)
#Re-order the groups because the default is 1yr-2w-8w
meta$treatment.ord = factor(meta$treatment, c("F_SHM", "F_OVX", "M_SHM", "M_OVX"))
levels(meta$treatment.ord)
#Return the plot area to 1x1
par(mfrow = c(1, 1))
#Plot
boxplot(pielou_e~ treatment, data=meta, ylab="pielou_e")

levels(meta$treatment)
par(mfrow = c(1, 1))
boxplot(pielou_e ~ treatment, data=meta, ylab="pielou_e")

pielou_e <- ggplot(meta, aes(treatment, pielou_e)) + 
  geom_boxplot(aes(color = treatment)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
pielou_e

#ANOVA and Tukey are tests based on the mean, but the boxplot plots the median. Its not wrong, its just not the best method. Unfortunately plotting the average and standard deviation is a little complicated.
observed_otus_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(treatment.ord) %>%   # the grouping variable
  summarise(mean_observed_otus = mean(observed_otus),  # calculates the mean of each group
            sd_observed_otus = sd(observed_otus), # calculates the standard deviation of each group
            n_observed_otus = n(),  # calculates the sample size per group
            se_observed_otus = sd(observed_otus/sqrt(n()))) # calculates the standard error of each group
#Make a barplot of means vs. treatment 
observed_otus_se <- ggplot(observed_otus_summary, aes(treatment.ord, mean_observed_otus, fill = treatment.ord)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_observed_otus - se_observed_otus, ymax = mean_observed_otus + se_observed_otus), width=0.2) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="observed_otus  ± s.e.", x = "") 
observed_otus_se

ggsave("output/evenness_se.png", evenness_se, height = 3, width = 3)

#We can test if the interaction of body site and days.since.experiment.start impacts diversity with a model that includes both of our variables. The `*` symbol is a shortcut for models. A*B is equivalent to A + B + A:B
aov.observed_otus.all = aov(observed_otus ~ treatment*sex, data=meta)
summary(aov.observed_otus.all)

aov.shannon.all = aov(shannon ~ treatment*sex, data=meta)
summary(aov.shannon.all)

aov.observed_otus.all = aov(observed_otus ~ sex*treatment, data=meta)
summary(aov.observed_otus.all)

#Now we will run the simplified model. This model statement only tests the main effects without the interaction: 
aov.observed_otus.all2 = aov(observed_otus ~ treatment+sex, data=meta)
summary(aov.observed_otus.all2)

aov.shannon.sex2 = aov(shannon ~ treatment+sex, data=meta)
summary(aov.shannon.sex2)

aov.evenness.all2 = aov(pielou_e ~ treatment+sex, data=meta)
summary(aov.evenness.all2)

#If all of our variables were categorical, we could run TukeyHSD like we did with body.site only.
TukeyHSD(aov.evenness.all2)
TukeyHSD(aov.faith.all)
TukeyHSD(aov.observed_otus.all)
TukeyHSD(aov.observed_otus.all2)
TukeyHSD(aov.observed_otus2)



BC.nmds = metaMDS(OTU.clean, distance="bray", k=2, trymax=1000)
#We see that we reached a convergent solution around 20 iterations and our stress is very low (<0.1), meaning that 2-axis are sufficient to view the data.
#if running line-by-line run this entire block at one time.
my_colors = c("blue", "green", "red", "black")
my_colors2 = c("blue", "green")
par(mfrow = c(1, 1))
#Create a blank plot for the nmds
plot(BC.nmds, type="n", main="Bray-Curtis") 
#Add the points colored by age
points(BC.nmds$points, display="sites", pch=20, col=my_colors[meta$treatment.ord])
#Add a legend
legend(-5, 2.5, legend=levels(meta$treatment.ord), col=my_colors, pch=20)

# Now with ggplot
BC.nmds$stress
# a lower stress means the spatial distances in your NMDS more accurately represent your calculated similarites. <0.2 is basically required, <0.1 is better but uncommon. If its bad, maybe you need to recalculate the mds line above and use transformed data 
#autotransform = TRUE

# I like to merge my NMDS coordinates in together with my metadata to make one big dataframe, I think this makes plotting easier later on

nmds <-as.data.frame(BC.nmds$points)
metanmds <- merge(meta, nmds, by.x =0, by.y = 0)
row.names(metanmds) <- metanmds[,1]
metanmds <- metanmds[,-1]
str(metanmds)
metanmds$treatment <- factor(metanmds$treatment)
# this generates a dataframe containing the group centroids

#The following line you may need to modify to indicate the column numbers that contain MDS1 and MDS2
NMDS.mean <- aggregate(metanmds[, 10:11], list(group=metanmds$sex), mean)
colnames(NMDS.mean) <- c('treatment', 'groupX', 'groupY')
# merging the group centroids with the rest of the NMDS data #
metanmds <- merge(metanmds, NMDS.mean, by.x = "sex", by.y="treatment")
str(metanmds)

ggplot(metanmds, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=treatment, shape=sex)) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(BC.nmds$stress, digits = 2))) +
  stat_ellipse(aes(color=treatment), level = 0.95) +
  theme(legend.title = element_blank()) 

ggsave("output/nmds_ellipses_all.png", height = 3, width = 4)

#A similar thing can be done for the Jaccard metric, which only takes into account presence/absence (*i.e.* richness).
J.nmds = metaMDS(OTU.clean, distance="jaccard", k=2, trymax=1000)

#if running line-by-line run this entire block at one time.
plot(J.nmds, type="n", main="Jaccard")
points(J.nmds, display="sites", pch=20, col=my_colors[meta$treatment.ord])
#Add a legend
legend(-5, 2.5, legend=levels(meta$treatment.ord), col=my_colors, pch=20)

#You see that the values are very different for Jaccard but the pattern of points is very similar to Bray-Curtis. This is because Jaccard is a transformation of Bray-Curtis with J = 2BC/(1+BC)
J.nmds$stress

# a lower stress means the spatial distances in your NMDS more accurately represent your calculated similarites. <0.2 is basically required, <0.1 is better but uncommon. If its bad, maybe you need to recalculate the mds line above and use transformed data 
#autotransform = TRUE

# I like to merge my NMDS coordinates in together with my metadata to make one big dataframe, I think this makes plotting easier later on

nmds <-as.data.frame(BC.nmds$points)
metanmds <- merge(meta, nmds, by.x = 0, by.y = 0)
row.names(metanmds) <- metanmds[,1]
metanmds <- metanmds[,-1]
str(metanmds)
metanmds$days.since.experiment.start <- factor(metanmds$treatment)



# this generates a dataframe containing the group centroids

#The following line you may need to modify to indicate the column numbers that contain MDS1 and MDS2
NMDS.mean <- aggregate(metanmds[,10:11], list(group=metanmds$treatment), mean)
colnames(NMDS.mean) <- c('design', 'groupX', 'groupY')

# merging the group centroids with the rest of the NMDS data #
metanmds <- merge(metanmds, NMDS.mean, by.x = "treatment", by.y="design")

str(metanmds)

#Plotting MDS
ggplot(metanmds, aes(x=MDS1, y=MDS2)) +
  geom_point(aes(color=treatment.ord, shape=sex)) +
  labs(x='NMDS 1', y= 'NMDS 2', caption = paste('Ordination stress: ', round(BC.nmds$stress, digits = 2))) +
  stat_ellipse(aes(color=treatment.ord), level = 0.95) +
  theme(legend.title = element_blank()) 

#3D plotting

#Calculate the Bray-Curtis nMDS for 3-axis
BC.nmds.3D = metaMDS(OTU.clean, distance="bray", k=3, trymax=1000)

#Extract x-y-z values for this nmds
BCxyz = scores(BC.nmds.3D, display="sites")
#This is a table that looks like 
BCxyz

#Plot the xyz coordinates and color by age
plot_ly(x=BCxyz[,1], y=BCxyz[,2], z=BCxyz[,3], type="scatter3d", mode="markers", color=meta$treatment.ord, colors=my_colors)

#Plotting 3D data in a 2D form  -> comming for journals
par(mfrow=c(1,2))
#Axis 1 and 2 (x and y)
plot(BCxyz[,1], BCxyz[,2], main="Bray-Curtis 1:2", pch=20, col=my_colors[meta$treatment.ord])
legend(-5.4, 0, legend=levels(meta$treatment.ord), col=my_colors, pch=20)
#Axis 1 and 3 (x and z)
plot(BCxyz[,1], BCxyz[,3], main="Bray-Curtis 1:3", pch=20, col=my_colors[meta$treatment.ord])

#qiime2R phylogenitic based metrics
library(qiime2R)
#metadata<-read.table("/Users/clapp5/Desktop/Microbiota_analysis_R/project/sample_info.tsv")
metadata = read.table("/Users/clapp5/Desktop/Microbiota_analysis_R/project/sample_info.tsv", header=TRUE, row.names=1, sep="\t")
wunifrac<-read_qza("/Users/clapp5/Desktop/Microbiota_analysis_R/project/weighted_unifrac_pcoa_results.qza")

vec = wunifrac$data$Vectors
vec %>%
  select(SampleID, PC1, PC2) %>%
  left_join(rownames_to_column(metadata), by = c('SampleID'='rowname')) %>%
  left_join(rownames_to_column(shannon), by = c('SampleID'='rowname')) %>%
  ggplot(aes(x=PC1, y=PC2, color=`treatment`, shape=`sex`)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(16,1), name="treatment") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="shannon") +
  scale_color_discrete(name="sex")
ggsave("output/wUF-PCoA.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave("output/wUF-PCoA_bigger.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches

shan = shannon$data$shannon.x
vec %>%
  select(SampleID, PC1, PC2) %>%
  left_join(rownames_to_column(metadata), by = c('SampleID'='rowname')) %>%
  left_join(rownames_to_column(shannon), by = c('SampleID'='rowname')) %>%
  ggplot(aes(x=PC1, y=PC2, color=`treatment`, shape=`sex`)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  stat_ellipse(aes(shape=`sex`), level = 0.95) +
  theme_q2r() +
  scale_shape_manual(values=c(16,1)) + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  #scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="Body Site")
ggsave("output/wUF-PCoA-ellipse.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave("output/wUF-PCoA-ellipse_bigger.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches


#if we want another axis
vec %>%
  select(SampleID, PC2, PC3) %>%
  left_join(rownames_to_column(metadata), by = c('SampleID'='rowname')) %>%
  #left_join(shannon) %>%
  ggplot(aes(x=PC2, y=PC3, color=`treatment`, shape=`sex`)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  stat_ellipse(aes(color=`sex`), level = 0.95) +
  theme_q2r() +
  scale_shape_manual(values=c(16,1)) + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  #scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="sex")
ggsave("output/wUF-PCoA-ellipse-2-3.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave("output/wUF-PCoA-ellipse-2-3_bigger.pdf", height=3, width=4, device="pdf") # save a PDF 3 inches by 4 inches



#For Bray Curtis
#if running line-by-line run this entire block at one time.

fit.BC = envfit(BC.nmds, meta) 
fit.BC

#plot(BC.nmds, type="n", main="Bray-Curtis")
#points(BC.nmds, pch=20, display="sites", col=my_colors[meta$treatment.ord])
#legend(-6, 2, legend=levels(meta$treatment.ord), col=my_colors, pch=20)
#Add fitted variables
#plot(fit.BC, col="black", p.max=0.01)

#plot(BC.nmds, type="n", main="Bray-Curtis")
#points(BC.nmds, pch=20, display="sites", col=my_colors[meta$treatment.ord])
#legend(-6, 2, legend=levels(meta$treatment.ord), col=my_colors, pch=20)
#Add fitted variables
#plot(fit.BC, col="black", p.max=0.05)
fit.BC.OTU = envfit(BC.nmds, OTU.clean[,1:10])
fit.BC.OTU
plot(BC.nmds, type="n", main="Bray-Curtis")
points(BC.nmds, pch=20, display="sites", col=my_colors[meta$treatment.ord])
legend(-6, 2, legend=levels(meta$treatment.ord), col=my_colors, pch=20)
#Add fitted variables
plot(fit.BC.OTU, col="black", p.max=0.05)

#Extract all OTUs within the genus Ruminococcus
tax.clean.complete = drop_na(tax.clean, Genus) 
OTU.Rumino = OTU.clean[,tax.clean.complete$Genus == "g__Ruminococcus"]
#Sum the abundances of the Ruminococcaceae OTUs into one variable (column)
OTU.Rumino$Rumino.sum = rowSums(OTU.Rumino)
#Statistically test Beta Diversity
#Calculate distance and save as a matrix
BC.dist=vegdist(OTU.clean, distance="bray")
#Run PERMANOVA on distances.
adonis(BC.dist ~ treatment.ord*sex, data = meta, permutations = 1000)
wUF.dist = UniFrac(physeq, weighted=TRUE, normalized=TRUE)

anosim(BC.dist, meta$treatment, permutations = 1000)

library(qiime2R)

physeq <- qza_to_phyloseq(
  features="/Users/clapp5/Desktop/Microbiota_analysis_R/project/rarefied_table.qza",
  tree="/Users/clapp5/Desktop/Microbiota_analysis_R/project/rooted-tree.qza",
  "/Users/clapp5/Desktop/Microbiota_analysis_R/project/taxonomy.qza",
  metadata = "/Users/clapp5/Desktop/Microbiota_analysis_R/project/sample_info.tsv"
)
plot_bar(physeq, fill="Phylum")
plot_bar(physeq, x="treatment",fill="Phylum")

plot_bar(physeq, x="treatment", fill="Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

#First get the OTU table from physeq
physeq_otu_table <- data.frame(otu_table(physeq))

#Assign as variables to be feed into phyloseq
OTU.physeq = otu_table(as.matrix(physeq_otu_table), taxa_are_rows=TRUE)

#our edited and formatted taxonomy table from the top of this script
tax.physeq = tax_table(as.matrix(tax.clean))    
meta.physeq = sample_data(meta)         

physeq_bar_plot = phyloseq(OTU.physeq, tax.physeq, meta.physeq)         
my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black"
)         
my_level <- c("Phylum", "Family", "Genus")
rm(taxa.summary)

ml ="Phylum"
for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(treatment, get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  names(taxa.summary)[2] <- ml
  
  physeq.taxa.average <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.average=mean(Abundance.average))
  names(physeq.taxa.average)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.average)
  
  abund_filter <- 0.01
  physeq_meta_filtered <- filter(physeq_meta, overall.average>abund_filter)
  str(physeq_meta_filtered)
  
  physeq_meta_filtered$body.site.ord = factor(physeq_meta_filtered$treatment, c("left palm", "right palm", "gut", "tongue"))
  
  # Plot 
  ggplot(physeq_meta_filtered, aes(x = treatment, y = Abundance.average, fill = get(ml))) + 
    #facet_grid(.~subject) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = my_colors) +
    # Remove x axis title
    #theme(axis.title.x = element_blank()) + 
    ylim(c(0,1)) +
    guides(fill = guide_legend(reverse = F, keywidth = .5, keyheight = .5, ncol = 1)) +
    theme(legend.text=element_text(size=8)) +
    #theme(legend.position="bottom") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(legend.title = element_blank()) +
    ylab("Relative Abundance") +
    ggtitle(paste0(ml, " (>", abund_filter * 100,"%) Composition of microbiome samples"))          
}         
#Check how many OTUs you have. If you have more than 1000, you may want to filter out rare OTUs, as shown below.
ntaxa(physeq)

# Set prunescale 
prunescale = 0.0001
seq_depth = sample_sums(physeq)[1]

tax.mean <- taxa_sums(physeq)/nsamples(physeq)

#Select the samples we want to compare
rf.samples = c("F_OVX", "F_SHM")

# Prune out rare OTUs by mean relative abundance set by prunescale
physeq.prune <- prune_taxa(tax.mean > prunescale*seq_depth, physeq)
physeq.prune2 <- prune_samples(sample_data(physeq.prune)$treatment == rf.samples[1] | sample_data(physeq.prune)$treatment == rf.samples[2], physeq.prune)
physeq.prune2

# Make a dataframe of training data with OTUs as column and samples as rows
predictors <- t(otu_table(physeq.prune2))
dim(predictors)

# Make one column for our outcome/response variable 
response <- as.factor(sample_data(physeq.prune2)$treatment)

# Combine them into 1 data frame
rf.data <- data.frame(response, predictors) 
str(rf.data)

set.seed(2)

physeq.classify <- randomForest(response~., data = rf.data, ntree = 100)
print(physeq.classify)

names(physeq.classify)

imp <- importance(physeq.classify)
imp <- data.frame(predictors = rownames(imp), imp)

# Order the predictor levels by importance
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)

#We need some editing to remove leading "X" in front of ASV names that should start with a numeral.
imp.sort$predictors <- str_remove(imp.sort$predictors, "[X]")

# Select the top 20 predictors
imp.20 <- imp.sort[1:20, ]

imp.20$predictors <-  factor(imp.20$predictors, levels = imp.20$predictors)

# ggplot
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle(paste("Most important ASVs between", rf.samples[1], "and ", rf.samples[2]))

# What are those OTUs?
otunames <- imp.20$predictors
r <- rownames(tax_table(physeq)) %in% otunames
rf.imp.asv <- as.data.frame(tax_table(physeq)[r, ])
imp.20_tax <- merge(imp.20, rf.imp.asv, by.x = "predictors", by.y = 0)

ggplot(imp.20_tax, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", aes(fill = Family)) +
  scale_fill_manual(values = my_colors) +
  coord_flip() +
  ggtitle(paste("Most important ASVs between", rf.samples[1], "and", rf.samples[2]))
ggsave(paste0("output/RandomForest_", rf.samples[1], "_", rf.samples[2], "_names.png"), height = 4, width = 7)
write.table(imp.20_tax, "output/rf_top20_taxa.tsv", sep="\t", row.names=FALSE, quote = FALSE)


library("DESeq2")
#To use DESeq, we need no zeros in our OTU table. So we will edit the table by multiplying by 2 and + 1
OTU.clean2 <- OTU.clean * 2 + 1
OTU.physeq = otu_table(as.matrix(OTU.clean2), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(tax.clean))
meta.physeq = sample_data(meta)
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)

diagdds = phyloseq_to_deseq2(physeq_deseq, ~ treatment)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

###Investigate test results table

my_contrast = c("treatment", "F_OVX", "F_SHM")
res = results(diagdds, contrast = my_contrast, cooksCutoff = FALSE)

alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq_deseq)[rownames(sigtab), ], "matrix"))

#Volcano Plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,15)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)}

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
DESeq_fig = ggplot(sigtab, aes(x=Genus, y = log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
DESeq_fig
ggsave(paste0("output/DESeq2-", my_contrast[2], "-", my_contrast[3], ".png"), DESeq_fig, height = 4, width = 10)

#Heat Maps
#Sort the OTUs by abundance and pick the top 20
top20OTU.names = names(sort(taxa_sums(physeq), TRUE)[1:20])
#Cut down the physeq data to only the top 10 Phyla
top20OTU = prune_taxa(top20OTU.names, physeq)
#We now see that we only have 20 taxa
top20OTU
#First, you can make a heatmap of OTU abundance across all samples
plot_heatmap(top20OTU)

plot_heatmap(top20OTU, sample.label="treatment", sample.order="treatment")
#We can label the OTU taxa
plot_heatmap(top20OTU, sample.label="treatment", sample.order="treatment", taxa.label="Genus")
#And group OTUs within the same Phyla
plot_heatmap(top20OTU, sample.label="treatment", sample.order="treatment", taxa.label="Genus", taxa.order="Phylum")

#alter color
plot_heatmap(top20OTU, sample.label="treatment", sample.order="treatment", taxa.label="Genus", taxa.order="Phylum", low="white", high="blue", na.value="grey")
#You can also have R automatically group your OTUs and samples by beta-diversity. This may yield the most easily interpreted heatmap but if you have a specific research question that is better addressed by your own ordering (like our age groups above), you should stick with that. We'll show Bray-Curtis as an example. Other options are
#bray
#jaccard
#wunifrac
#uwunifrac
plot_heatmap(top20OTU, "NMDS", "bray", title="Bray-Curtis")
#Bray-Curtis
heatmap.2(as.matrix(BC.dist))
#UniFrac-- not functioning yet
heatmap.2(as.matrix(wUF.dist))
#changing colors
rc <- rainbow(nrow(as.matrix(BC.dist)), start=0, end=0.9)
heatmap.2(as.matrix(BC.dist), col=rc)


#Venn diagrams
OTU.5017.2w = colnames(OTU.clean["5017.2w.F", apply(OTU.clean["5017.2w.F",], MARGIN=2, function(x) any(x >0))])

OTU.5017.8w = colnames(OTU.clean["5017.8w.F", apply(OTU.clean["5017.8w.F",], MARGIN=2, function(x) any(x >0))])

OTU.5017.1yr = colnames(OTU.clean["5017.1yr.F",apply(OTU.clean["5017.1yr.F",], MARGIN=2, function(x) any(x >0))])


#network
plot_net(prevotella, color="Species", type="taxa")

#co-occurance of OTUs
plot_net(physeq, color="treatment", distance="bray")

#publication figures
postscript("Fig1.ps", width = 7, height = 3, horizontal = FALSE, colormodel = "rgb", family = "ArialMT")

layout(matrix(c(1,2), 1, 2), widths=c(3,2), heights=c(1,1))

plot(BC.nmds, type="n", main="Bray-Curtis")
points(BC.nmds, display="sites", pch=20, col=c("blue", "green", "red")[meta$treatment])

boxplot(observed_otus ~ treatment.ord, data=meta, main="Diversity", ylab="Shannon's diversity", col=c("green", "red", "blue"))

dev.off()






