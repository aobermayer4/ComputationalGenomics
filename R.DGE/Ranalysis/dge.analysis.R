# Differential gene-expression analysis using DEseq2
## authors: Chengqi Wang, J. Oberstaller

## set source (installs packages if required, loads CW's go_enrichment function for later)
source("Rsource/source.R")

### Load count data from featureCounts

input <- "Rdata/vehicle_drug_feature_counts.txt"
df <- read.table(input,
                      header = T,
                      sep = '\t',
                      row.names = 1)

df <- df[,6:9]


### Check the distribution of counts per sample

par(mfrow = c(2,2), mar = c(4,4,1,1))
print( apply(df, 2, function(x) quantile(as.numeric(x))) )
hist.plots <- apply(df, 2, function(x) {hist(log2(x), breaks = 100)})

## Low-expressed genes can make for noisy data. We'll filter out the lowest-expressed genes (bottom 25%), 
## which corresponds to less than ~16 (or 2^4) reads.

### Keep genes with expression bigger than 2^4

idExp <- apply(df, 1, function(x) any(x > 16))
dfExp <- df[idExp, ]
print(dim(dfExp))


## now we have filtered data ready for normalization and expression-analysis with DESeq2.

### Build DEseq class

library('DESeq2')
table2 <- data.frame(name = colnames(dfExp),
                     condition = c('control','control', 'treatment', 'treatment'))
dds <- DESeqDataSetFromMatrix(dfExp, 
                              colData=table2, design= ~ condition)
dds <- DESeq(dds)

#normalized reads count
norCounts <- counts(dds, normalized=TRUE)

#DEseq results
res <- results(dds)


####Getting differentially expressed genes
# take a look at our results object ("res")
res
# res is a dataframe object, so we can check out metadata for what the columns mean
mcols(res, use.names=TRUE)

## extract the genes with adjusted P-value < 0.01
resSig <- res[ which(res$padj < 0.01), ]
dim(resSig)


# sort by log2FoldChange to get the significant genes with the strongest down-regulation
head( resSig[ order( resSig$log2FoldChange ), ] )

# and strongest up-regulation
tail( resSig[ order( resSig$log2FoldChange ), ] )

resSig.sorted.df <- as.data.frame(resSig[ order( resSig$log2FoldChange ), ])
# write the full list to a sorted, tab-delimited table for your records. The file will be created in your Routput folder.
write.table(resSig.sorted.df,
            file = "Routput/resSig.sorted.tab.txt",
            sep = "\t",
            row.names = TRUE,
            quote = FALSE)


##heatmap plot
##plot heatmap of normalized read-counts
library("RColorBrewer")
library("gplots")

sigNorData <- norCounts[which(res$padj < 0.01),]
hmcol <-  colorRampPalette(brewer.pal(9, "GnBu"))(100)
## expression heatmap
heatMapDF_nor <- t( apply(sigNorData, 1, function(x){(x-mean(x))/sd(x)}) )

colnames(heatMapDF_nor) <- c('control1','control2',
                             'treat1'  , 'treat2')

## comment out lines 83 & 85 to view your plot in your "plots" window, or 
## leave as-is to save your plot to .pdf
pdf('Rfigs/heatmap.pdf', height = 10, width = 10)
heatmap.2(heatMapDF_nor, col = hmcol, trace="none", margin=c(10, 10),labRow = F)
dev.off()

### Volcano plots

# Generate a data frame with a column assigning color to significant deferentially expressed genes
res_plot      <- data.frame( res )
res_plot$col  <- 'gray40'

# setting cutoffs for significantly up (red) and down (blue) genes; everything else is gray
res_plot$col[res_plot$log2FoldChange > 1 & res_plot$padj < 0.01] <- 'red'
res_plot$col[res_plot$log2FoldChange < -1 & res_plot$padj < 0.01] <- 'cornflowerblue'


## comment out lines 99 & 105 to view your plot in your "plots" window, or leave as-is to save your plot to .pdf
pdf('Rfigs/volcano.pdf', height = 10, width = 10)
plot(res_plot$log2FoldChange,
     -log10(res_plot$padj),
     col = res_plot$col, pch = 19, xlab = 'log2(fold change)',
     ylab = '-log10(p-adj)'
)
dev.off()
# Functional enrichment analysis

## load GO gaf file (GO annotation)
geneGO <- read.delim2('Rdata/PlasmoDB-48_Pfalciparum3D7_GO.gaf.txt', header = F, sep = '\t')
geneGO <- geneGO[,c(2,5)]


## load ontology
library(GO.db)
go_db <- Term(GOTERM)
go_On <- Ontology(GOTERM)
go_inf<- data.frame(ontology = go_On,
                    term     = go_db)

## now we'll run CW's go_enrichment function perform a Fisher test to find any GO terms 
## that appear in our differentially expressed genes more often than would be expected by 
## chance (see "RSource/source.R" for details)
goEnrichment <- go_enrichment(rownames(resSig)[resSig$log2FoldChange > 1],
                              rownames(df),
                              geneGO)

idMatch <- match(rownames(goEnrichment), rownames(go_inf))
goEnrichment <- data.frame(goEnrichment,
                           go_inf[idMatch, ])

# list the enriched GO terms
goEnrichment$term[goEnrichment$padj < 0.1]



### bar plot for go enrichment

goSig   <- goEnrichment[goEnrichment$padj < 0.01,]
goSig$expection <- goSig$X.BackgroundWithGOterm/(goSig$X.BackgroundWithGOterm +
                                                   goSig$X.BackgroundWithoutGOterm) * (goSig$X.QueryWithGOterm + 
                                                                                         goSig$X.QueryWithoutGOterm)

goSig <- goSig[order(goSig$X.QueryWithGOterm), ]
goSigDraw <- t( goSig[,c(1,9)] )
colnames(goSigDraw) <- goSig$term


## comment out lines 147 & 150 to view your plot in your "plots" window, or leave as-is to save your plot to .pdf
pdf('Rfigs/GO_barplot.pdf', height = 10, width = 10)
par(mar = c(4,20,1,1))
barplot(goSigDraw, horiz = T, las = 1)
dev.off()

## upload the three plot .pdfs you generated (saved in your "Rfigs" directory; heatmap.pdf, volcano.pdf, and GO_barplot.pdf) to Canvas to complete your assignment.
