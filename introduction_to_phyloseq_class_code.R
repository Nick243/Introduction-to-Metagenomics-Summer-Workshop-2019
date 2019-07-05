


#######################################################

#Introduction to phyloseq

#Created by: Nicholas Ollberding
#On 7/10/19


#######################################################



#Loading required packages and phyloseq object
library(dada2); packageVersion("dada2")            #version: 1.12.1
library(phyloseq); packageVersion("phyloseq")      #version: 1.28.0
library(ggplot2); packageVersion("ggplot2")        #version: 3.1.1

ps <- readRDS("ps.rds")



#Accessing the sample information and sample metadata
nsamples(ps)
sample_names(ps)

sample_variables(ps)
head(sample_data(ps))
sample_data(ps)$When
table(sample_data(ps)$When)
median(sample_data(ps)$Day)

metadata <- data.frame(sample_data(ps))
head(metadata)



#Examining the number of reads for each sample
sample_sums(ps)
sort(sample_sums(ps))
hist(sample_sums(ps), main="Histogram: Read Counts", xlab="Total Reads", 
     border="blue", col="green", las=1, breaks=12)

metadata$total_reads <- sample_sums(ps)



#Examining the OTU table
ntaxa(ps)
head(taxa_names(ps))
head(taxa_sums(ps))
(asv_tab <- data.frame(otu_table(ps)[1:5, 1:5]))



#Examining the taxonomy 
rank_names(ps)
head(tax_table(ps))
head(tax_table(ps)[, 2])
table(tax_table(ps)[, 2])
(tax_tab <- data.frame(tax_table(ps)[50:55, ]))



#Examining the reference sequences
head(refseq(ps))
dada2::nwhamming(as.vector(refseq(ps)[1]), as.vector(refseq(ps)[2]))
(ref_tab <- data.frame(head(refseq(ps))))



#Agglomerating and subsetting taxa
(ps_phylum <- tax_glom(ps, "Phylum"))
taxa_names(ps_phylum)
taxa_names(ps_phylum) <- tax_table(ps_phylum)[, 2]
taxa_names(ps_phylum)
otu_table(ps_phylum)[1:5, c(1:3, 5, 7)]



#Subsetting taxa
(ps_bacteroides <- subset_taxa(ps, Genus == "Bacteroides"))
tax_table(ps_bacteroides)
prune_taxa(taxa_sums(ps) > 100, ps) 
filter_taxa(ps, function(x) sum(x > 10) > (0.1*length(x)), TRUE)   



#Subsetting samples and tranforming counts
ps_early <- subset_samples(ps, When == "Early")
(ps_early = prune_taxa(taxa_sums(ps_early) > 0, ps_early))
sample_data(ps_early)$When
sort(sample_sums(ps))
(ps_reads_GT_5k = prune_samples(sample_sums(ps) > 5000, ps))
sort(sample_sums(ps_reads_GT_5k))

ps_relabund <- transform_sample_counts(ps, function(x) x / sum(x))
otu_table(ps_relabund)[1:5, 1:5]

(ps_rare <- rarefy_even_depth(ps, sample.size = 4000, rngseed = 123, replace = FALSE))
sample_sums(ps_rare)



##EXAMPLE ANANLYSIS and GRAPHICAL CAPABILITES 
#Alpha-diversity
head(estimate_richness(ps))
(p <- plot_richness(ps, x = "When", color = "When", measures = c("Observed", "Shannon")))
p + labs(x = "", y = "Alpha Diversity Measure\n") + 
  theme_bw()



#Beta-diversity ordination
ps_rare_bray <- ordinate(ps_rare, "NMDS", "bray")
plot_ordination(ps_rare, ps_rare_bray, type="samples", color="When") + geom_point(size = 3) 



#Bar plots of composition
plot_bar(ps, fill="Phylum")

plot_bar(ps_relabund, fill="Phylum") + 
  geom_bar(aes(color = Phylum, fill = Phylum), stat="identity", position="stack") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(panel.background = element_blank())



#Heatmaps
(ps_fam <- tax_glom(ps, "Family"))
(ps_fam_rare <- rarefy_even_depth(ps_fam, sample.size = 4000, rngseed = 123, replace = FALSE))
plot_heatmap(ps_fam_rare, sample.label = "When", taxa.label = "Family")






