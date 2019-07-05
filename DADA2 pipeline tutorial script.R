

#############################################################

#DADA2 Pipeline Tutorial (version 1.12)

#notes: - tutorial and data found at: https://benjjneb.github.io/dada2/tutorial.html
#       - instructions to install DADA2 provided at: https://benjjneb.github.io/dada2/dada-installation.html
#       - running DADA2 requires R version 3.6.0 and Bioconductor version 3.9
#       - pipeline requires DADA2, phyloseq and Biostring packages to be installed (and all dependencies)
#       - example fastq files and reference databases provided in MiSeq_SOP and tax folders



#Tutorial and data accessed on: 06/05/19


##############################################################




#Loading required packages
library(dada2); packageVersion("dada2")            #version: 1.12.1
library(phyloseq); packageVersion("phyloseq")      #version: 1.28.0


#Set working directory
setwd("dada2 tutorial")       #change this to point to the folder where you have saved the seq data and tax databases


#Listing path to sequence files
path <- "MiSeq_SOP"           #folder with the fastq files
list.files(path)


#Getting matched lists of the forward and reverse fastq files and extracting sample names
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)                 #extracting string prior to the first underscore

head(fnFs)                   #sanity check that the files are in order and the names are as expected
head(fnRs)
sample.names


#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


#Filter and trim forward and reverse reads
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))    #Place filtered files in new filtered subdirectory
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names               #add names
names(filtRs) <- sample.names

head(filtFs)                                #sanity check

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)                    #On Windows set multithread=FALSE (may need to throughout)
head(out)


#Learn error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
errF$err_out


#Sample inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[1]]


#Merge F and R reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])


#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


#Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))


#Filter chimeric reads
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim) 
sum(seqtab.nochim)/sum(seqtab)    #ensuring retained majority of reads


#Tracking reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)



#Assigning  taxonomy
#For more information and examples see: https://benjjneb.github.io/dada2/assign.html
taxa <- assignTaxonomy(seqtab.nochim, "tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)   #assign taxonomy using RDP classifier and Silva ref database
taxa <- addSpecies(taxa, "tax/silva_species_assignment_v132.fa.gz")                            #updating with genus/species binomials

taxa.print <- taxa                     #Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)



#Evaluate accuracy on mock sample
unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")



#Create phyloseq object
samples.out <- rownames(seqtab.nochim)                            #cleaning metadata for import to phyloseq
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
samdf

(ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa)))

(ps <- prune_samples(sample_names(ps) != "Mock", ps)) # Remove mock sample


dna <- Biostrings::DNAStringSet(taxa_names(ps))       #convert ASV names to ASV1, ASV2, ASV...
names(dna) <- taxa_names(ps)
(ps <- merge_phyloseq(ps, dna))
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

saveRDS(ps, "ps.rds")








