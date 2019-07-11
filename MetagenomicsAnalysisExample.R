

###SCRIPT FROM DR. HASLAM'S LECTURE 7/8/19

##CODE AND DATA PROVIDED BY DR. HASLAM

### Metagenomics Course Example Analysis ####


# Create a Folder on your desktop called 'MetagenomicsCourse'
# Put this R script, the sample key, and the AllSpeciesRaw.csv file in this folder
# make sure you can find this folder when you set working directory with setwd(".....") command below

# install and load required packages, with a little extra work to install Devtools

# ## Devtools installation:
#   # devtools is needed for the NBZIMM packages, which is required for negative binomial zero inflated modle
#   # We're following instructions from here: https://www.r-project.org/nosvn/pandoc/devtools.html
# # Windows:
#   # on windows, first install rtools by pasting this into browser: https://cran.r-project.org/bin/windows/Rtools/Rtools34.exe
# # Mac:
#   # may first need to install xcode: https://apps.apple.com/us/app/xcode/id497799835?mt=12
# # Linux: 
#   # may need to install a compiler but most versions have one preinstalled
# 
# # then do the following:
#   install.packages("devtools")
#   # Answer 'yes' to questions about downloading packages
#   # then do the following:
#   devtools::install_github("hadley/devtools")
# 
#   
# # now we can install the NBZIMM package
#   devtools::install_github("nyiuab/NBZIMM", build_opts = c("--no-resave-data", "--no-manual"), force = T)
  
  # install the other packages as needed:

 load.lib<-c("R.rsp", "reshape2", "plyr", "vegan", "ggplot2", "MASS", "FactoMineR", 
  "factoextra", "dunn.test", "NBZIMM")
 
 install.lib<-load.lib[!load.lib %in% installed.packages()]
 for(lib in install.lib) install.packages(lib,dependencies=TRUE)
 sapply(load.lib,require,character=TRUE)
 
# and get them into your working environment

library(reshape2)
library(plyr)
library(vegan)
library(ggplot2)
library(MASS)
library(FactoMineR)
library(factoextra)
library(dunn.test)
library(NBZIMM)
 
# Set working directory some useful paramater ## change working directory to your desktop folder

#setwd("C:/Users/dbhas/OneDrive/Desktop/MetagenomicsCourse")    #should noto be needed if you clone the repo

# generalized log2 function
glog2 <- function(x) ((asinh(x) - log(2))/log(2))


paramsBox<- function(x){theme(axis.text.x= element_text(size= 14, color="black")) +
    theme(axis.text.y = element_text(size= 18, color="black")) +
    theme(plot.title = element_text(size= 22, color="black")) +
    theme(axis.title.x = element_text(size=20), axis.title.y = element_text(size=22)) +
    theme(legend.title = element_text(size=18)) +
    theme(legend.text = element_text(size = 14)) +
    theme(legend.position="none") +
    theme(axis.text.x = element_text(size= 18, color="black")) +
    theme(axis.text.y = element_text(size= 16, angle=0, color="black")) +
    theme(strip.text.y = element_text(size = 14, angle = 270, colour = "black")) +
    theme(strip.text.x = element_text(size = 16, angle = 0, colour = "black"))
}

# Set some baseline colors
example.colors<-c("#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f", "#edc948", "#b07aa1",
                  "#ff9da7", "#9c755f", "#bab0ac")


# Noise removal function - the percent cutoff can be changed here, or when calling the function later
# most relevant to distance and PCoA calculation (e.g. set to 1% to get a few species, 0.01% to get more species)
noise.removal <- function(dataframe, percent=0.001, top=NULL){
  dataframe->Matrix
  AboveCutoff <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
  Matrix_1 <- Matrix[AboveCutoff,]
  print(percent)
  return(Matrix_1)
}


# Get the metadata. This table holds the information about patient (mouse) groups and other factors that might be 
# important for comparison between samples, such as age at sample collection, source of sample, clinical variables, etc., etc. 

MetaData<-read.csv("SimpleMetaData.csv", header = TRUE, stringsAsFactors = FALSE)
MetaData$Group<-factor(MetaData$Group, levels = c("Control", "Drug.1", "Drug.2"))
MetaData$SampleTime<-factor(MetaData$SampleTime, levels = c("Time.1", "Time.2"))


    # Next we want to get the taxonomy abundance data. Depending on the alignment step, the taxonomy count data
    # may come in different forms. From Kraken (and Bracken), the output is one file per sample
    # so we first need to get each of those files and merge them. We're goig to skip that step, but the code is below
    
    
    # This is how to generate the raw species count table from individual alignment files
    # First we use Kraken to align raw reads to a genome database
    # Optionally can then use Bracken to 'reassign' reads that were ambiguously assigned at higher taxonomic level
    # This little loop gets all of the files we need and combines them into one dataframe
    
    
    # setwd("~/Documents/Alignments/KrakenAlignments/")
    # 
    # # Find all the Kraken alignment files 
    # AllKrakenFiles<-list.files()
    # SpeciesFileList<-grep("_species_abundance.txt", AllKrakenFiles)
    # SpeciesFiles<-AllKrakenFiles[SpeciesFileList]
    # FileList<-gsub("_species_abundance.txt", "", SpeciesFiles)
    # 
    # # Get all the species count files
    # 
    # for(f in 1:length(SpeciesFiles)){
    #   fnr = SpeciesFiles[f]
    # 
    #   # assign is a neat little command that will call your data.frame whatever you want, then read in the data
    #   assign(fnr, read.csv(paste(fnr, sep=""), sep="\t", header = TRUE, stringsAsFactors = FALSE))
    # }
    # 
    # 
    # RawSpeciesCounts<-read.csv(paste(SpeciesFiles[1]), sep = "\t", header = TRUE)
    # names(RawSpeciesCounts)[1]<-"Species"
    # RawSpeciesCounts<-as.data.frame(RawSpeciesCounts[,c(1,6)])
    # names(RawSpeciesCounts)<-c("Species",paste(paste(FileList[1])))
    # 
    # # read in the rest of the files
    # for ( i in 2:length(SpeciesFiles)){
    #   x <- get( SpeciesFiles[i] )
    #   x<-as.data.frame(x[,c(1,6)])
    #   x<-subset(x, ! duplicated(x$Species))
    #   RawSpeciesCounts<-merge(x, RawSpeciesCounts, by = "Species", all = TRUE)
    #   
    # }
    # 
    # # add a column called Species
    # row.names(RawSpeciesCounts)<-RawSpeciesCounts$Species
    # # in case there are any NA in the table, replace them with 0
    # RawSpeciesCounts[is.na(RawSpeciesCounts)]<-0
    # AllSpeciesRaw<-t(RawSpeciesCounts)


  # So now we have a table called RawSpeciesCounts that has all the raw count data.
  # Since each sample may have different number of total reads, we can't just work with this data for most analysis
  # We have to normalize the counts in some way to make the total number of reads the same for every sample
  # however, for some analysis, that's not Groupessary and may lead to data loss, so in those cases (e.g. diversity calculations)
  # we willl use the raw count table



# However, rather than have a bunch of individual Kraken output files, Skip to already combined raw species count table:

AllSpeciesRaw<-read.csv("AllSpeciesRaw.csv", header = TRUE, stringsAsFactors = FALSE)
dim(AllSpeciesRaw)

# so there are 169 samples and 6209 species
# our database has human genome, so there will be human reads here
# we may want them for some projects, but not this one, so get rid of that column
# we also spike in Salniibacter as an internal standard so we ahve to get rid of those reads

HumanReads<-grep("Homo.sapiens", names(AllSpeciesRaw))
AllSpeciesRaw<-AllSpeciesRaw[, -HumanReads]
SaliniReads<-grep("Salinibacter.ruber", names(AllSpeciesRaw))
AllSpeciesRaw<-AllSpeciesRaw[, -SaliniReads]

# make the rownames the sample ID because we're going to have to get rid of Sample column at some point

SpeciesTable<-as.data.frame(AllSpeciesRaw)
row.names(SpeciesTable)<-SpeciesTable$Sample




# There are several decisions about data analysis from now on that vary between investigators and projects 
# Choices about filtering species based on abundance and getting rid of samples with low reads differ between investigtors
# also, normalization choice varies (many use rarefaction, which we're goint to use)
# statistical methods to compare overall microbiome structure and patient (mouse) groups also vary. 
  # In general non-parametric methods are used for microbiome count data
  # but depending on the investigator, the model used to simulate data distribution differs
  # recently, zero-inflated negative binomial regression has become common for microbiome data
  # but there's some evidence that other data distribution models are as good or better for certain data
  # e.g. Gaussian , Poisson, and other models
  # All of that being said, when compared side-to-side there are generally not huge differences in results
  # IMHO if there are real differences in the data, they will probably show up, regardless of how the count data distribution is modeled
  # one just needs to be consistent and clear about the choice of models and analytic methods


# So here we go: First filter species to be present in at least 10% of samples and account for > 0.01% of counts
  # This is to get rid of rare species that probably don't contribute much to overall difference 
  # and helps keep statistical analysis more robust, cuz we have to account for testing of many species
  # if there are a lot of low abundance species that likely don't contribute to differences in patient groups, 
  # they're going to hurt us during the 'multiple testing' correction. So get rid of them now
  # note though, this may not be appropriate if you think rare species may account for the disease/phenotype of interest

RawSpecies<-as.data.frame(SpeciesTable)
RawSpecies$Sample<-NULL

# Make sure we're only looking at samples for which we have metadata
SampleNamesPresent<-which(row.names(RawSpecies) %in% MetaData$Sample == TRUE)
RawSpecies<-RawSpecies[SampleNamesPresent,]

# Find  species which are present in > 10% samples
  # note this cutoff is easy to change, if we want to keep species present in at least 20%, etc of samples
  # 'floor' is a command that rounds down a number with decimal points. We want an integer here
  # TenPercentCutoff is just a calculation of the number of samples that constitute 10% of all samples
  # most of this code is re-useable, so that next time, if there are more species and more samples, we can use the same code


TenPercentCutoff<-floor(nrow(RawSpecies)/10)

NonZeroCounts<-list() # initiate an empty list. We can't call it in the loop below otherwise
# this loop goes column by column (species by species) to count the number of samples (rows) that have > 0 of that species
for (i in 1:ncol(RawSpecies)){
  NonZeroCounts[i]<-length(which(RawSpecies[,i] > 0))
}

# Which ones are greater than our cutoff of 10%, and keep only those columns (species)

TenPercentNotZero<-which(NonZeroCounts >= TenPercentCutoff)
SpeciesTableAboveCutoff<-RawSpecies[,TenPercentNotZero]

# Now there are only 1444 species left
# Get rid of samples with less than 1 million reads

LowSamples<-which(rowSums(SpeciesTableAboveCutoff) <= 1000000)
SpeciesCut<-SpeciesTableAboveCutoff[-LowSamples,]

# now noise removal, using the function we created up top:
# here we're choosing 0.005 cutoff which gets us ~ 500 species
SpeciesReduced<-as.data.frame(t(noise.removal(t(SpeciesCut), 0.001)))

# SpeciesNR is a dataframe with noise reduced Species (only species in > 10% of samples,
  # and only species accounting for at least 0.01% of overall read count. Also we're getting rid of samples
  # that have less than one million aligned reads which may be too stringent for some projects

# Now rarefying to 1 million reads, but can also rarefy to the number of reads in the sample with lowest number (minCount)
minCount<-min(rowSums(SpeciesReduced))
SpeciesNR<-data.frame(rrarefy(SpeciesReduced, 1000000))
SpeciesNR$Sample<-row.names(SpeciesNR)


# Ok now the data is normalized but for some analysis (e.g. diversity calculations and distance measures) we can use non-rarefied data
  # this may be better because with rarefication you can lose some information
  # in practice, I haven't found much difference when using rarefied or non-rarefied data for diversity or distance measures
  # but as a general principle, it's good to stick with non-rarefied when possible


# So we go back to the non-rarefied raw species table. We'll keep all species but still get rid of samples with low read counts
LowSamplesNonRarefied<-which(rowSums(RawSpecies) <= 1000000)
RawSpeciesNoLow<-RawSpecies[-LowSamplesNonRarefied,] # remember that RawSpecies is the dataframe before noise reduction and rarefaction
RawSpeciesNoLow$Sample<-row.names(RawSpeciesNoLow)

# add on the metadata, though in this example set I don't think we'll ever use this table
NonRarefiedSpecies<-merge(MetaData, RawSpeciesNoLow, by = "Sample", all.y = TRUE)

# Get some diversity metrics. Most people use Shannon Index but there are several options, most of which are here:
  # BTW, we keep adding and removing a column called "Sample" because we need that identifier to merge metadata
  # but it gets in the way for calculations, like the diversity calculation below
  # so we have to get rid of it by making it <-NULL

DiversitySpecies<-RawSpeciesNoLow
row.names(DiversitySpecies)<-DiversitySpecies$Sample
DiversitySpecies$Sample<-NULL

# Calculate sample diversity using 'diversity' function in vegan

H <- data.frame(diversity(DiversitySpecies))
simpson <- data.frame(diversity(DiversitySpecies, "simpson"))
shannon<-data.frame(diversity(DiversitySpecies, "shannon"))
invsimp <- data.frame(diversity(DiversitySpecies, "inv"))
alpha <- data.frame(fisher.alpha(DiversitySpecies))
## Species richness (S) and Pielou's evenness (J):
S <- data.frame(specnumber(DiversitySpecies))
J <- data.frame(H/log(S))
Diversity<-cbind(simpson, shannon, invsimp, alpha, S, J)
Diversity$Sample<-row.names(Diversity)
names(Diversity)<-c("Simpson", "Shannon",  "InvSimpson", "Alpha", "SpeciesNo", "Evenness", "Sample")
# reorder the columns
Diversity<-Diversity[,c(7,1,2,3,4,5,6)]

# this will let us look at how different diversity measures are related to each other
pairs(cbind(simpson, shannon, alpha, J, S), pch="+", col="blue")


# alrighty, add back the metadata to our new Diversity table so we can make some graphs
Diversity<-merge(MetaData, Diversity, by = "Sample", all.y  = TRUE)

# graph Shannon Index for each sample as boxplot, split by group

col = example.colors[1:3]

ShannonGraph<-ggplot(Diversity, 
                     aes(x=Group, y=as.numeric(Shannon ), fill=Group)) + geom_boxplot(lwd=1,aes(color=factor(Group),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Group))) + xlab(NULL) + 
  ylab("Shannon Diversity Index \n") + paramsBox() + scale_y_log10() + xlab(NULL) 
ShannonGraph
ggsave(filename = "ShannonGraph.pdf", plot = ShannonGraph, width = 6, 
       height = 8, limitsize = FALSE)
pairwise.wilcox.test(Diversity$Shannon, Diversity$Group)

col = example.colors[1:3]

# is there any difference based on sampling time?
ShannonGraph2<-ggplot(Diversity, 
                     aes(x=Group, y=as.numeric(Shannon ), fill=Group)) + geom_boxplot(lwd=1,aes(color=factor(Group),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Group))) + xlab(NULL) + facet_grid(facets = . ~ SampleTime) +
  ylab("Shannon Diversity Index \n") + paramsBox() + scale_y_log10() + xlab(NULL) 
ShannonGraph2
ggsave(filename = "ShannonGraph2.pdf", plot = ShannonGraph2, width = 10, 
       height = 8, limitsize = FALSE)


# Graph the number of species per sample by groups

SpeciesGraph<-ggplot(Diversity, 
                     aes(x=Group, y=as.numeric(SpeciesNo), fill=Group)) + geom_boxplot(lwd=1,aes(color=factor(Group),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Group))) + xlab(NULL) + 
  ylab("Number of Species \n") + paramsBox() + scale_y_log10() + xlab(NULL) 
SpeciesGraph
ggsave(filename = "SpeciesGraph.pdf", plot = SpeciesGraph, width = 6, 
       height = 8, limitsize = FALSE)

# check for statistical significance between patient groups
pairwise.wilcox.test(Diversity$SpeciesNo, Diversity$Group)

# How about that. We have significant difference in alpha-diversity between our  groups 
  # are there any differences in overall microbiome composition (probably, if alpha-diversity is different)
  # we'll use MRPP to figure that out and visualize by PCA 

# After that, we'll check for individual species differences
  # we can use several different ways to figure that out. Here we'll do a screening test using
  # a zero-inflated negative binomal model of the data which can account for other clinical factors, including time (repeated measures from the same patient)
  # this is called a zero-inflated negative binomial mixed model
  # from a practical standpoint, pairwise wilcoxon rank sum or Kruskal-Wallis are also often used and are more simple


## For the following analysis we're going to use the rarefied species count data
  # For PCoA plot we're going to also transform the data using the generalized log2 calculation function 
  # created at the beginning

SpeciesNR<-merge(MetaData, SpeciesNR, by = "Sample", all.y = TRUE)

# check to see how many samples in each group

table(SpeciesNR$Group)

# MRPP to check for overall differences in microbiome:

Group.mrpp<-mrpp(SpeciesNR[,5:ncol(SpeciesNR)], SpeciesNR$Group, distance = "bray")
Group.mrpp # this gives us a significance of 0.001


# compare MRPP using  transformed data:
SpeciesTransformed<-glog2(SpeciesNR[,5:ncol(SpeciesNR)])
GroupTransformed.mrpp<-mrpp(SpeciesTransformed, SpeciesNR$Group)
GroupTransformed.mrpp # also significant at 0.001


# Let's do PCA to look at overall composition distribution and differences between groups:

metadata <- SpeciesNR[,1:4]
cts <- as.matrix(SpeciesNR[, -(1:4)])
rownames(cts) <- metadata$Sample

cts_l2 <- glog2(cts)
grps <- SpeciesNR$Group

col = example.colors[1:3]
newPCA<-PCA(cts_l2, scale.unit = FALSE, ncp = 5, graph = FALSE)
pcaPlot<-fviz_pca_ind(newPCA,
                      geom.ind = "point", 
                      pointsize = 1.5,
                      point.alph = 0.1,
                      title = "Unsupervised Principal Coordinate Analysis",
                      subtitle = "Samples Colored by Patient Group",
                      col.ind = grps, 
                      addEllipses = TRUE, 
                      ellipse.alpha = 0.2,
                      ellipse.type = "confidence",
                      ellipse.level = 0.95,
                      legend.title = "Patient Group", 
                      legend.size = 11,
                      mean.point = FALSE,
                      palette = col,
                      axes.linetype = "blank"
)

GroupPCA<- pcaPlot + scale_fill_manual(values = col)
GroupPCA

# save the graph
pdf("GroupPCA.pdf") 
print(GroupPCA)
dev.off() 

# So overall microbiome composition is different between our patient groups.
# now we want to figure out if there are individual species that are significantly different between groups
# the main considerations here are: a) how do we model the data? (zero-inflated and overdispersed (variance > mean) binomial, zero-inflated gaussian, etc)
# b) do we want to use log transformed or non-transformed. A little trial and error is fine. Generally there are not huge differences.

# First let's use straightforward Wilcoxon rank sum test to compare each species abundance in the two groups
# we're using normalized (same number of counts per sample) but not transformed (log or scaled) data. Some people
# would probably use transformed data here as you remove the compositional nature of the data
# it's worth trying both ways and I defer to Nick's lectures on using count or transformed data for pairwise comparison

# Let's find out which species differ most between Control and Drug.2

GroupSpecies<-subset(SpeciesNR, SpeciesNR$Group %in% c("Control", "Drug.2"))
GroupSpecies<-GroupSpecies[, c(1,3,5:ncol(GroupSpecies))]


Wilcox<-pairwise.wilcox.test(GroupSpecies[,3], GroupSpecies[,2])
WilcoxTable<-as.data.frame(Wilcox$p.value)
names(WilcoxTable)<-names(GroupSpecies)[3]
row.names(WilcoxTable)<-"Wilcox"


for (i in 4:length(colnames(GroupSpecies))){
  x<-pairwise.wilcox.test(GroupSpecies[,i], GroupSpecies[,2])
  y<-as.data.frame(x$p.value)
  names(y)<-names(GroupSpecies)[i]
  # names(p)<-names(GroupSpecies)[i]
  WilcoxTable<-cbind(WilcoxTable, y)
}

# Correct for multiple testing (more species in the analysis means less likely to survive correction)
WilcoxTable<-as.data.frame(t(WilcoxTable))
names(WilcoxTable)<-c("Unadjusted_p")
# there are other options for p.adust function, including "bonferroni", "hommel", "holm", etc...
WilcoxTable$FDR<-p.adjust(WilcoxTable$Unadjusted_p, method = "fdr")
WilcoxTable$Species<-row.names(WilcoxTable)
WilcoxTable<-WilcoxTable[order(WilcoxTable$Unadjusted_p),]

# save this table
write.csv(WilcoxTable, file = "WilcoxTable.csv")

# we can find only species that have FDR < 0.05, which is actually very strict
SigSpecies<-subset(WilcoxTable, WilcoxTable$FDR <= 0.05)

# wow, there are a ton of differences
# let's look at a few of those species to see if they really look different between groups



# Graphs of individual species. We'll go back to SpeciesNR so we can look at all three groups


SmitisGraph<-ggplot(SpeciesNR, 
       aes(x=Group, y=as.numeric(Streptococcus.mitis ), fill=Group)) + geom_boxplot(lwd=1,aes(color=factor(Group),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Group))) + xlab(NULL) + 
  ylab("Streptococcus mitis\n") + paramsBox() + scale_y_log10() + xlab(NULL) + 
  theme(axis.title.y = element_text(face = "italic"))
SmitisGraph
ggsave(filename = "SmitisGraph.pdf", plot = SmitisGraph, width = 6, 
       height = 8, limitsize = FALSE)

pairwise.wilcox.test(SpeciesNR$Streptococcus.mitis, SpeciesNR$Group)


SaureusGraph<-ggplot(SpeciesNR, 
       aes(x=Group, y=as.numeric(Staphylococcus.aureus ), fill=Group)) + geom_boxplot(lwd=1,aes(color=factor(Group),fill = NA), outlier.size = 3)  +
  stat_summary(fun.y=mean, geom="point", shape=5, size=8) + scale_fill_manual(values = col) + scale_colour_manual(values = col ) +
  geom_point(size=4,aes(color = factor(Group))) + xlab(NULL) + 
  ylab("Staphylococcus aureus\n") + paramsBox() + scale_y_log10() + xlab(NULL) + 
  theme(axis.title.y = element_text(face = "italic"))
SaureusGraph
ggsave(filename = "SaureusGraph.pdf", plot = SaureusGraph, width = 6, 
       height = 8, limitsize = FALSE)
pairwise.wilcox.test(SpeciesNR$Staphylococcus.aureus, SpeciesNR$Group)



# this takes 10 to 20 mins to run so we'll skip it for class
# negative binomial, zero inflated mixed model
# 
# SampleTime = SpeciesNR$SampleTime
# Group = SpeciesNR$Group
# subject = SpeciesNR[, "PatientID"]; table(subject)
# 
# SpeciesData<-SpeciesNR[,5:ncol(SpeciesNR)]
# N<-rowSums(SpeciesData)
# non = nonzero(y = SpeciesData, total = N, plot = T)
# nonzero.p = non[[1]]
# 
# y = SpeciesData[, names(nonzero.p)[3]]
# 
# f = mms(y = SpeciesData, fixed = ~  SampleTime + Group + offset(log(N)),
#         random = ~ 1 | subject, correlation = corAR1(),
#         method = "zinb")
# 
# res = get.fixed(f, part="dist",vr.name="SampleTimeFinish",  sort.p=T)
# par(mfrow = c(1, 1), cex.axis = 1, mar = c(2, 10, 4, 4))
# plot.fixed(res, threshold=0.01, gap=300, main="Covariate: Patient Group",
#            cex.axis=1, cex.var=1)




