#----Prepartion----#

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.11")
library(dada2)

path <- "~/R/ComputationalGenomics/biocomputingclass/MiSeq_SOP" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#----Inspect read quality profiles----#

plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

#----Filter and trim----#

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Place filtered files in filtered/ subdirectory
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#----Learn the error rate----#

#learn the error rate for forward reads(16*41 matrix)
# 41 is cresponding to the Quality value
errF <- learnErrors(filtFs, multithread=TRUE)

#learn the error rate for reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

#plot the results
plotErrors(errF, nominalQ=TRUE)

#----Sample inference----#

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

#----Merge paired reads----#

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#----Construct ASV sequence table----#

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#----Remove chimeras----#

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#----Assign taxonomy by naive Bayesian classifier----#

taxa <- assignTaxonomy(seqtab.nochim, "~/R/ComputationalGenomics/biocomputingclass/silva_nr99_v138_wSpecies_train_set.fa.gz", multithread=TRUE)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#----Evaluate accuracy----#

unqs.mock <- seqtab.nochim["Mock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")

mock.ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match.ref <- sum(sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

#----How to calculate Shannon diversity----#

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender  <- substr(subject,1,1)
subject <- substr(subject,2,999)
day     <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf   <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out
samdf$col <- rep('red', nrow(samdf))
samdf$col[samdf$When == 'Early'] <- 'navy'

#input should be samples x taxonomy
shannonDiversity <- function(df) 
{
  p <- t( apply(df, 1, function(x){x/sum(x)}) )
  H <- apply(p , 1, function(x){x <- x[x > 0]; -sum( log(x) * x )})
  return(H)
}

SimpsonDiversity <- function(df) 
{
  p <- t( apply(df, 1, function(x){x/sum(x)}) )
  H <- apply(p , 1, function(x){x <- x[x > 0];1-sum( (x ^ 2))})
  return(H)
}
shD <- shannonDiversity(seqtab.nochim)
siD <- SimpsonDiversity(seqtab.nochim)
par(mfrow = c(1,2))
plot(samdf$Day,shD, col = samdf$col, pch = 19 ,las = 1, 
     xlab = 'Days', ylab = 'Shannon Diversity')
plot(samdf$Day,siD, col = samdf$col, pch = 19 ,las = 1, 
     xlab = 'Days', ylab = 'Simpson Diversity')

ps.prop <- t(apply(seqtab.nochim,1,function(x){x/sum(x)}))

##install package vegan
##install.packages('vegan')
library(vegan)

NMDS <- metaMDS(ps.prop[-20,],     # remove mock line
                distance = 'bray') # The number of reduced dimensions

plot(NMDS$points[,1], NMDS$points[,2], col = samdf$col,
     pch = 19, las = 1, xlab = 'MDS1', ylab = 'MDS2')

dim(taxa)
dim(seqtab.nochim)

idTop20   <- order(colSums(seqtab.nochim), decreasing = T)[1:20]
top20     <- ps.prop[-nrow(seqtab.nochim),idTop20] #remove Mock
top20Taxa <- data.frame( taxa[idTop20, ] )

#encode color for top20Taxa$Family
library("RColorBrewer")
colNeed <- data.frame(1:9, brewer.pal(9, "Set1"))
top20Taxa$category <- as.numeric(as.factor(top20Taxa$Family))
top20Taxa$category[is.na(top20Taxa$category)] <- 7

top20Taxa$color    <- colNeed[match(top20Taxa$category, colNeed[,1]),2]
samdfUse <- samdf[-nrow(samdf),]
top20prop <- top20[, order(top20Taxa$Family, decreasing = T)]
top20Taxa <- top20Taxa[order(top20Taxa$Family,decreasing = T), ]

#pdf('~/Desktop/aa.pdf', height = 8,width = 10)
par(mfrow = c(1,3))
barplot( t(top20prop[samdfUse$When == 'Early',]), col = top20Taxa$color,
         las = 2)
barplot( t(top20prop[samdfUse$When == 'Late',]), col = top20Taxa$color,
         las = 2)

##for legend
legendMatrix <- unique(cbind(top20Taxa$Family,
                             top20Taxa$color))
par(mar = c(0,0,0,0))
plot(NA,NA, xlab = '', ylab = '',
     xaxt = 'n', yaxt = 'n',xlim = c(0,1), ylim = c(0,1),
     bty = 'n')
legend('center', legend = legendMatrix[,1],
       col = legendMatrix[,2], pch = 15,
       bty = 'n')