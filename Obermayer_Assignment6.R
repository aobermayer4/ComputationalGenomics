#Get File
df <- read.delim("~/R/ComputationalGenomics/ComputationalGenomics/vehicle_drug_feature_counts.txt", comment.char="#")
df <- vehicle_drug_feature_counts

#Load Packages
library(ggplot2)
library(gplots)
library(reshape2)
library(gridExtra)

#----------Part 1----------#

#Generate FPKM data frame
df_fpkm <- apply(df[,7:10], 2,
                 function(x){
                   x/df$Length * 10^9/sum(x)
                 })

df_fpkm <- data.frame(df_fpkm) #generate data frame
colnames(df_fpkm) <- colnames(df[,7:10]) #assign column names
head(df_fpkm)

#Convert FPKM data frame values into log2 scale
df_fpkm_log2 <- log2(df_fpkm + 0.01) #use small number, b/c without number it takes log to infinity
head(df_fpkm_log2)

#----------Part 2----------#

#Generate scatter plot for each pair of biological replicates

scatterControl <- ggplot(df_fpkm_log2, 
                         aes(x = vehicle_rep1.bam,
                             y = vehicle_rep2.bam)) +
  geom_point(shape = 21, color = 'black', bg = rgb(235, 207, 196, max = 255)) +
  ggtitle("log2 FPKM Control") +
  xlab("log2(Control1_FPKM)") +
  ylab("log2(Control2_FPKM)")
scatterControl

scatterTreatment <- ggplot(df_fpkm_log2, 
                         aes(x = drug_rep1.bam,
                             y = drug_rep2.bam)) +
  geom_point(shape = 21, color = 'black', bg = rgb(172, 130, 149, max = 255)) +
  ggtitle("log2 FPKM Treatment") +
  xlab("log2(Treatment1_FPKM)") +
  ylab("log2(Treatment2_FPKM)")
scatterTreatment

ScatterPlots <- grid.arrange(scatterControl, scatterTreatment, nrow = 1)

#----------Part 3----------#

#Generate box plots according to FPKM values of each assay
par(mar = c(6,4,1,1)) #set margins
boxplot <- boxplot(df_fpkm_log2$vehicle_rep1.bam,
        df_fpkm_log2$vehicle_rep2.bam,
        df_fpkm_log2$drug_rep1.bam,
        df_fpkm_log2$drug_rep2.bam,
        main = "log2 FPKM Boxplot",
        names = c('vehicle_rep1', 'vehicle_rep2',
                  'drug_rep1', 'drug_rep2'),
        ylab = 'log2(FPKM)',
        las = 2, #axis label orientations
        col = c(rgb(235, 207, 196, max = 255),
                rgb(235, 207, 196, max = 255),
                rgb(172, 130, 149, max = 255),
                rgb(172, 130, 149, max = 255)),
        medcol = rep('gray30', 5))

#----------Part 4----------#

#Generating histograms according to FPKM values of each assay

histFPKMc1 <- ggplot(df_fpkm_log2,
       aes(x = vehicle_rep1.bam)) + 
  geom_histogram(binwidth = 1, 
                 color = 'gray',
                 fill = rgb(235, 207, 196, max = 255)) +
  scale_x_continuous(name = 'vehicle_rep1') +
  scale_y_continuous(name = 'Frequency')
histFPKMc2 <- ggplot(df_fpkm_log2,
       aes(x = vehicle_rep2.bam)) + 
  geom_histogram(binwidth = 1, 
                 color = 'gray',
                 fill = rgb(235, 207, 196, max = 255)) +
  scale_x_continuous(name = 'vehicle_rep1') +
  scale_y_continuous(name = 'Frequency')
histFPKMt1 <- ggplot(df_fpkm_log2,
       aes(x = drug_rep1.bam)) + 
  geom_histogram(binwidth = 1, 
                 color = 'gray',
                 fill = rgb(172, 130, 149, max = 255)) +
  scale_x_continuous(name = 'drug_rep1') +
  scale_y_continuous(name = 'Frequency')
histFPKMt2 <- ggplot(df_fpkm_log2,
       aes(x = drug_rep2.bam)) + 
  geom_histogram(binwidth = 1, 
                 color = 'gray',
                 fill = rgb(172, 130, 149, max = 255)) +
  scale_x_continuous(name = 'drug_rep2') +
  scale_y_continuous(name = 'Frequency')

HistPlots <- grid.arrange(histFPKMc1, histFPKMc2, histFPKMt1, histFPKMt2, 
                          nrow = 2,
                          top = "log2 FPKM Histograms: Control vs. Treatment")

#----------Part 5----------#

corFPKM <- cor(df_fpkm)
colnames(corFPKM) <- c('c1','c2','t1','t2')
rownames(corFPKM) <- c('c1','c2','t1','t2')
colPlate <- colorRampPalette(c(rgb(235, 207, 196, max = 255),
                               rgb(244, 238, 225, max = 255),
                               rgb(77, 36, 61, max = 255)))(100) #returns 100 colors between the palette given
heatmapFPKM <- heatmap.2(corFPKM,
          Colv = F, #Determines reordering of col, if true columns treated as rows
          dendrogram = 'row', #draws the tree on row/col but defaults to both
          trace = 'n', #trace line on the color key, defaults to cyan, n makes that false
          denscol = 'black', #density display color
          srtCol = 0, #angel of row/col labels
          col = colPlate, #color palette described above
          key.title = "Pearson Correllation",
          main = "FPKM Heatmap: Control vs. Treatment")

