# 031420 hands-on differential expression analysis: start-up and functions for analysis
# Author: J. Oberstaller

##### Start-up environment ####
# install necessary general packages/dependencies if they aren't already installed
packages <- c("BiocManager","gplots","knitr","RColorBrewer")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos='https://cloud.r-project.org') 
}

# install our primary analysis-packages from bioconductor if they aren't already installed
bioconductor.packages <- c("GO.db", "DESeq2")
if (length(setdiff(bioconductor.packages, rownames(installed.packages()))) > 0){
  BiocManager::install(bioconductor.packages,
                       ask = FALSE,
                       quiet = TRUE,
                       verbose = FALSE)
}

# clean up leftover variables from above
rm("packages", "bioconductor.packages")


## build function for go enrichment
go_enrichment <- function(geneID_query,
                          geneID_background,
                          geneAnnotation) 
{
  geneQueryInf      <- geneAnnotation[geneAnnotation[,1] %in% geneID_query,]
  geneBackgroundInf <- geneAnnotation[geneAnnotation[,1] %in% geneID_background,]
  queryGOterm <- unique(geneQueryInf[,2])
  
  goResult <- c()
  aa <- sapply(queryGOterm, function(id) 
  {
    numQuery        <- length( which(geneQueryInf[,2] == id) )
    numBackground   <- length( which(geneBackgroundInf[,2] == id))
    
    numQuery_no     <- length(geneID_query) - numQuery
    numBackground_no<- length(geneID_background) - numBackground
    
    #print(c(numBackground, numBackground_no))
    fishTest <- fisher.test(rbind( c(numQuery, numQuery_no),
                                   c(numBackground, numBackground_no) ),
                            alternative = 'greater')
    infReturn <- c(numQuery,
                   numQuery_no,
                   numBackground,
                   numBackground_no, fishTest$p.value)
    goResult <<- rbind(goResult, infReturn)
  })
  rownames(goResult) <- queryGOterm
  colnames(goResult) <- c('#QueryWithGOterm',
                          '#QueryWithoutGOterm',
                          '#BackgroundWithGOterm',
                          '#BackgroundWithoutGOterm', 'pvalue')
  goResult <- data.frame(goResult)
  goResult$padj <- p.adjust(goResult$pvalue, method = 'fdr')
  return(goResult)
}
