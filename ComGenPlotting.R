#Load packages
library(ggplot2)
library(gplots)
library(reshape2)

#Load data
data("airquality")

#---------Boxplot---------#

#Base Boxplot
boxplot(airquality$Ozone ~ airquality$Month,
        names = c('May','June','July','Aug','Sept'),
        xlab = 'Month', las = 2)

#Create boxplot, no outliers
boxplot(airquality$Ozone ~ airquality$Month,
        names = c('May','June','July','Aug','Sept'),
        xlab = 'Month', 
        las = 2, #0-3 orientation of labels on asix
        outline = F, #Takes away outlier dots outside of box plot
        ylab = 'Ozone',
        col = c(rgb(187, 181, 164, max = 255),
                rgb(235, 207, 196, max = 255),
                rgb(244, 238, 225, max = 255),
                rgb(172, 130, 149, max = 255),
                rgb(77, 36, 61, max = 255)),
        medcol = rep('gray30', 5))



#ggplot boxplot

aq <- airquality

#Change month column from numbers to names of month
aq$Month <- factor(airquality$Month,
                   labels = c('May','June','July','Aug','Sept'))

pp <- ggplot(aq, aes(x = Month, y = Ozone)) + 
  geom_boxplot(fill = 'pink') +
  scale_y_continuous(name = 'Mean Ozone in Parts per Billion',
                     breaks = seq(0,175,25))
pp



#-------Histogram---------#

#Base R histogram
par(mar = c(4,3,3,1)) #Setting graphical parameters, helpful if you need more space
hist(aq$Ozone,
     breaks = 50,
     col = 'cornflowerblue',
     las = 1,
     xlab = 'Ozone',
     main = 'Frequency Histogram if Mean Ozone')



#ggplot version

ggplot(aq, aes(x = Ozone)) + 
  geom_histogram()

#Specific Months, only for August
ggplot(aq[aq$Month == 8, ], aes(x = Ozone)) + 
  geom_histogram()

#ggplot histogram, make it fancy
barfill <- 'gold'
barlines <- 'gray'
ggplot(aq, aes(x = Ozone)) + 
  geom_histogram(binwidth = 5, #width of bars
                 color = barlines,
                 fill = barfill) +
  scale_x_continuous(name = 'Mean Ozone in Parts per Billion',
                     breaks = seq(0,175,25),
                     limits = c(0,175)) +
  ggtitle('Frequency Histogram of Mean Ozone')


#-----------ScatterPlot------------#

#get data
data("mtcars")


#Base R
par(mar = c(4,4,1,1))
plot(mtcars$mpg, mtcars$wt,
     xlab = 'mpg',
     ylab = 'wt',
     las = 1,
     col = 'cornflowerblue',
     pch = 19, #point type
     cex = 1.5, #size of points
     cex.axis = 1.5, #size of axis labels
     cex.lab = 1.5) #size of labels
abline(lm(mtcars$wt ~ mtcars$mpg),
       lwd = 3, #line width
       col = 'forestgreen')

#ggplot

ggplot(mtcars, aes(x = mpg, y = wt)) +
  geom_point(size = 2,
             shape = 23) +
  geom_text(label = rownames(mtcars))

ggplot(mtcars, aes(x = mpg, y = wt, size = cyl)) +
  geom_point() +
  geom_smooth(method = lm, #linear regression
              se = T) #T adds gray background to line



#trouble with with
#ggplot(mtcars, aes(x = mpg, y = wt)) +
  geom_point(aes(size = qsec))


#---------HeatMap-----------#
  
#base R
heatmap(as.matrix(mtcars)) #requires to be matrix

#gplots
heatmap.2(as.matrix(mtcars),
          trace = 'none', #takes vertical lines away
          adjCol = c(1,1),
          offsetRow = 0,
          offsetCol = 0,
          srtRow = 0,
          srtCol = 90)

heatmap.2(as.matrix(mtcars),
          trace = 'none',
          scale = 'col') #normalizing data

#ggplot requires 2 dimensional long list 
#create correlation data frame and transpose it
cor_df <- cor(t(mtcars))

#melt converts data frame into long format
melt_corDF <- melt(cor_df)

#checking dimensions of data frame
dim(cor_df)
dim(melt_corDF)

ggplot(data = melt_corDF,
       aes(x = Var1, y = Var2, 
           fill = value)) + #fill based off value
       scale_fill_gradient2(low = 'navy', 
                            high = 'red', 
                            mid = 'white', 
                            midpoint = 0.925, 
                            limit = c(0.85,1)) +
  geom_tile(color = 'white') #color=white surrounds the tiles with white outline


#using class materials
#already loaded data

#normalize data
df_fpkm <- apply(df[,6:9], 2,
                 function(x){
                   x/df$Length * 10^9/sum(x)
                 })

df_fpkm <- data.frame(df_fpkm)
colnames(df_fpkm) <- colnames(df[,6:9])
head(df_fpkm)

#making scatter plot with FPKM control data
plot(log2(df_fpkm$vehicle_rep1.bam), log2(df_fpkm$vehicle_rep2.bam),
     las = 1,
     pch = 21, #19 w/ col for plain circles, 21 with bg to make black circles with background color
     bg = 'red',
     xlab = 'log2(control1_FPKM)',
     ylab = 'log2(control2_FPKM)')

#making boxplot

#need to convert FPKM value into log2 scale
df_fpkm_log2 <- log2(df_fpkm + 0.01) #use small number because without number it takes log to infinity
head(df_fpkm_log2)

#make boxplot
par(mar = c(6,4,1,1)) #set margins
boxplot(df_fpkm_log2$vehicle_rep1.bam,
        df_fpkm_log2$vehicle_rep2.bam,
        df_fpkm_log2$drug_rep1.bam,
        df_fpkm_log2$drug_rep2.bam,
        names = c('control 1', 'control2',
                  'treatment1', 'treatment2'),
        ylab = 'log2(FPKM + 0.01',
        las = 2, #axis label orientations
        col = c('red','red','navy','navy'))


#make histograms
par(mfrow = c(2,2))
hist(df_fpkm_log2$vehicle_rep1.bam)
hist(df_fpkm_log2$vehicle_rep2.bam)
hist(df_fpkm_log2$drug_rep1.bam)
hist(df_fpkm_log2$drug_rep2.bam)

hist_own <- function(val, xlabel, mainTitle){
  hist(val, col = 'cornflowerblue',
       xlab = xlabel,
       main = mainTitle,
       las = 1,
       xlim = c(-5,15)) #all x-axis have the same limit
}
par(mfrow = c(2,2))
par(mar = c(4,4,1,1))
hist_own(df_fpkm_log2$vehicle_rep1.bam,
         'log2(FPKM + 0.01)',
         'control1')
hist_own(df_fpkm_log2$vehicle_rep2.bam,
         'log2(FPKM + 0.01)',
         'control2')
hist_own(df_fpkm_log2$drug_rep1.bam,
         'log2(FPKM + 0.01)',
         'treatment1')
hist_own(df_fpkm_log2$drug_rep2.bam,
         'log2(FPKM + 0.01)',
         'treatment2')


corDF <- cor(df_fpkm)
colnames(corDF) <- c('c1','c2','t1','t2')
rownames(corDF) <- c('c1','c2','t1','t2')
colPlate <- colorRampPalette(c('red','white','navy'))(100) #adding 100 seems to give it smoother gradient
heatmap.2(corDF,
          Colv = F, #Determines reordering of col, if true columns treated as rows
          dendrogram = 'row', #draws the tree on row/col but defaults to both
          trace = 'n', #trace line on the color key, defaults to cyan, n makes that false
          denscol = 'yellow', #density display color
          srtCol = 0, #angel of row/col labels
          col = colPlate)












