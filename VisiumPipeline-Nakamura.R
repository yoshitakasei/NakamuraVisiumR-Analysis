# 10X Visium Spatial Transcriptomics Analysis Script
# Nakamura Lab
# v1.0.0

#### Required Libraries to run the script ####
library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
library(matrixStats)
library(numbers)
library(ggrepel)
library(ggpubr)
library(gridExtra) 



#### Pre-processing steps ####
# Load the sample information file that describes the experimental design
VisiumSample <- read_excel('SampleInfo.xlsx')

# Key functions needed for the compareGroups
totalExpressionAuto <- function(VisiumSample){
  #From sample info file, load the gene expression data for each pad
  csv <- list()
  
  for(x in 2:length(VisiumSample[3,])){
    if(!is.na(VisiumSample[3,x])){
      csv[[x-1]] <- read_csv(as.character(VisiumSample[3,x]))  
    }
  }
  
  
  for(x in 1:length(csv)){
    csv[[x]] <- csv[[x]] %>%
      select(FeatureID, FeatureName, contains('Average'))
  }
  
  for(x in 1:length(csv)){
    colnames(csv[[x]]) <- paste(colnames(csv[[x]]), x, sep = "_")
  }
  for(x in 1:length(csv)){
    colnames(csv[[x]])[1] <- 'FeatureID'
    colnames(csv[[x]])[2] <- 'FeatureName'
  }
  
  total <- data.frame(matrix(ncol=2,nrow=0))
  colnames(total) <- c("FeatureID","FeatureName")
  
  for(x in 1:length(csv)){
    total <- merge(total,csv[[x]],by=c("FeatureID","FeatureName"),all=TRUE)
  }
  
  return(total) 
}
regionalExpression <- function(totalExpression, VisiumSample){
  regionalValues <- list()
  for(x in 2:length(VisiumSample[4,])){
    if(!is.na(VisiumSample[4,x])){
      regionalValues[[x-1]] <- totalExpression %>%
        select(FeatureID, FeatureName, contains(as.character(VisiumSample[4,x])))
    }
  }
  return(regionalValues)
}
getRegions <- function(VisiumSample){
  regions <- list()
  for(x in 2:length(VisiumSample[4,])){
    if(!is.na(VisiumSample[4,x])){
      regions[x-1] <- VisiumSample[4,x]
    }
  }
  return(regions)
}
getGroups <- function(VisiumSample){
  group <- list()
  for(x in 2:length(VisiumSample[5,])){
    if(!is.na(VisiumSample[5,x])){
      group[x-1] <- VisiumSample[5,x]
    }
  }
  return(group)
}
getGroupRange <-function(group, regionalExpression){
  groupRange <- grep(group, colnames(regionalExpression))
  return(groupRange)
}
getRegionHitstt <- function(region,threshold,ranget,rangec){
  
  region <- region %>%
    #adds column to df for t-test p-value (between columns in ranget and rangec)
    mutate(PValue = apply(region, 1, function(x) tryCatch({
      (t.test(as.numeric(x[rangec]),as.numeric(x[ranget]),var.equal = TRUE)$p.value)
    }, error = function(e) {NA}))) %>%
    #adds column to df for -log(p-value)
    mutate(neglog10PVal = -log10(PValue)) %>%
    #adds column to df for boolean indicating if p-value is <0.05
    mutate(pHit = ifelse((neglog10PVal>-log10(0.05)),TRUE,FALSE)) %>%
    #adds column to df for fold change
    mutate(fold= log2(rowMeans(.[, ranget])/rowMeans(.[, rangec]))) %>%
    #adds column to df for signal to noise ratio (SNR)
    mutate(SNR= apply (region,1, function (x) tryCatch({
      as.numeric((mean(as.numeric(x[ranget])) - mean(as.numeric(x[rangec]))) / (sd(as.numeric(x[rangec])) + sd(as.numeric(x[ranget]))))
    }, error = function(e) {NA}) ) ) %>%
    #adds column to df for fold score (fold change * p-value) 
    mutate(foldscore= fold * neglog10PVal) %>%
    #adds column to df for SNR score (SNR * p-value)
    mutate(SNRscore= SNR * neglog10PVal) %>%
    #arranges df such that p-value is from lowest to highest
    arrange(desc(neglog10PVal)) 
  
  #new df "pHits" filters df "region" for p-value hits
  pHits <- region %>%
    filter(pHit==TRUE) %>%
    select(FeatureID, FeatureName,pHit,foldscore, SNRscore)
  
  #
  SNRneg <- pHits %>%
    filter(SNRscore<0 & foldscore!=-Inf) %>%
    select(FeatureID, FeatureName,SNRscore) %>%
    arrange(SNRscore)
  
  #
  SNRpos <- pHits %>%
    filter(SNRscore>0 & foldscore!=Inf) %>%
    select(FeatureID, FeatureName,SNRscore) %>%
    arrange(desc(SNRscore))
  
  #
  SNRnegthreshold <- SNRneg$SNRscore[nrow(SNRneg)*threshold]
  SNRposthreshold <- SNRpos$SNRscore[nrow(SNRpos)*threshold]
  
  #
  foldneg <- pHits %>%
    filter(foldscore<0 & foldscore!=-Inf) %>%
    select(FeatureID, FeatureName,foldscore) %>%
    arrange(foldscore)
  
  #
  foldpos <- pHits %>%
    filter(foldscore>0 & foldscore!=Inf) %>%
    select(FeatureID, FeatureName,foldscore) %>%
    arrange(desc(foldscore))
  
  #
  foldnegthreshold <- foldneg$foldscore[nrow(foldneg)*threshold]
  foldposthreshold <- foldpos$foldscore[nrow(foldpos)*threshold]
  
  #
  region <- region %>%
    mutate(SNRHit = ifelse((pHit==TRUE)&((SNRscore>SNRposthreshold)|(SNRscore<SNRnegthreshold)),TRUE,FALSE)) %>%
    mutate(foldHit = ifelse((pHit==TRUE)&(foldscore!=-Inf)&(foldscore!=Inf)&((foldscore>foldposthreshold)|(foldscore<foldnegthreshold)),TRUE,FALSE)) %>%
    mutate(overlapHit = ifelse((SNRHit==TRUE)&(foldHit==TRUE),TRUE,FALSE))
  
  #
  foldHits <- region %>%
    filter(foldHit==TRUE) %>%
    select(FeatureID, FeatureName,foldscore,foldHit)
  
  #
  SNRHits <- region %>%
    filter(SNRHit==TRUE) %>%
    select(FeatureID, FeatureName,SNRscore,SNRHit)
  
  #
  overlapHits <- region %>%
    filter(overlapHit==TRUE) %>%
    select(FeatureID, FeatureName, neglog10PVal, foldscore, SNRscore)
  
  #
  return(list(region,foldHits,SNRHits,overlapHits,foldnegthreshold,foldposthreshold,SNRnegthreshold,SNRposthreshold,SNRpos,SNRneg,foldneg,foldpos))
  
}

# The main function to identify differentially expressed genes between 2 groups within a region of interest
compareGroups <- function(VisiumSample, region, groupExp, groupCtrl, threshold){
  #all data
  total <- totalExpressionAuto(VisiumSample)
  
  #all data sorted into regions
  regionValues <- regionalExpression(total, VisiumSample)
  
  #values from target region
  regionIDs <- getRegions(VisiumSample)
  regionIndex <- match(region, regionIDs)
  targetRegionValues <- regionValues[[regionIndex]] %>%
    select(FeatureID, FeatureName, contains(as.character(regionIDs[[regionIndex]])))
  
  #values from target groups
  groupIDs <- getGroups(VisiumSample)
  groupExpRange <- getGroupRange(as.character(groupExp), targetRegionValues)
  groupCtrlRange <- getGroupRange(as.character(groupCtrl), targetRegionValues)
  
  hitGenes <- getRegionHitstt(targetRegionValues, threshold, groupExpRange, groupCtrlRange)
  
  return(hitGenes)
}

# Function to list the hits identified in the comparedGroups function
getHitList <- function(comparedGroups){
  hitList <- as.data.frame(comparedGroups[4]) %>%
    select(FeatureName)
  return(hitList)
}

# Function to graph volcano plots with highlighted hits
plotregionvolcano <- function(regionmaster, highlight){
  regionmaster <- as.data.frame(regionmaster[1])
  genelist <- as.list(geneSample$Genes)
  
  if(highlight=="hits"){
    regionplotdf <- regionmaster %>%
      select(FeatureID, FeatureName, neglog10PVal, fold, overlapHit) %>%
      filter(fold!=Inf&fold!=-Inf) %>%
      mutate(plotname=ifelse(overlapHit==TRUE,FeatureName,""))
    
    
    return_item <- ggplot(regionplotdf, aes(fold,neglog10PVal, colour = overlapHit,label=plotname, size=overlapHit)) +
      geom_point() +
      geom_text_repel(size=3, max.overlaps = Inf) +
      scale_size_manual(values=c(1,2))
    
    return(return_item)
  }
  else if(highlight=="goi"){
    regionplotdf <- regionmaster %>%
      select(FeatureID, FeatureName, neglog10PVal, fold, overlapHit) %>%
      filter(fold!=Inf&fold!=-Inf) %>%
      mutate(plotname=ifelse(FeatureName %in% genelist, FeatureName,"")) %>%
      mutate(goi=ifelse(FeatureName %in% genelist, TRUE, FALSE))
    
    regionplotdf <- regionplotdf[order(regionplotdf$goi),]
    
    
    return_item <- ggplot(regionplotdf, aes(fold,neglog10PVal, colour = goi,label=plotname, size=goi)) +
      geom_point() +
      geom_text_repel(size=4, max.overlaps = Inf) +
      scale_size_manual(values=c(1,2)) 
    
    return(return_item)
  }
}

# Function to format gene expression data frame for plotting with ggplot
formatGeneExpression <- function(VisiumSample, region, gene, comparedGroups, maxN){
  #get values for all genes
  geneValues <- comparedGroups[[1]]
  
  #identify groups of columns for each group
  groupIDs <- getGroups(VisiumSample)
  groupRange <- list()
  for(x in 1:length(groupIDs)){
    groupRange[[x]] <- getGroupRange(as.character(groupIDs[[x]]), geneValues)
  }
  
  #get average gene expression for a single gene in given region
  geneRow <- filter(geneValues, FeatureName==gene)
  geneAvgs <- data.frame()
  if(nrow(geneRow) > 0){
    for(x in 1:length(groupIDs)){
      for(y in 1:maxN){
        if(!is.na(groupRange[[x]][y])){
          geneAvgs[x, y] <- geneRow[groupRange[[x]][y]]
        }else{
          geneAvgs[x, y] <- NA
        }
      }
    }
  }else{
    for(x in 1:length(groupIDs)){
      for(y in 1:maxN){
        geneAvgs[x, y] <- NA
      }
    }
  }
  
  #make data frame for graphing
  genePlotDim <- length(groupIDs)*maxN
  z <- 1
  expression <- data.frame()
  for(x in 1:length(groupIDs)){
    for(y in 1:maxN){
      expression[z,1] <- geneAvgs[x,y]
      expression[z,2] <- as.character(groupIDs[[x]])
      expression[z,3] <- as.character(region)
      z <- z+1
    }
  }
  colnames(expression) <- c("expression", "group", "region")
  
  return(expression)
  
}

# Functions to plot a series of bar graphs
plotBarGene <- function(df, gene, region){
  df_summary <- df %>%
    group_by(group) %>%
    summarise(
      sd = sd(expression, na.rm = TRUE),
      expression = mean(expression, na.rm = TRUE)
    )
  
  plot <- ggplot(df, aes(group, expression)) +
    ggtitle(gene, subtitle = region) + 
    geom_col(data = df_summary, fill = NA, color = "black") +
    geom_jitter(width = 0.2, height = 0, color = "black") + 
    geom_errorbar( aes(ymin = expression-sd, ymax = expression+sd), 
                   data = df_summary, width = 0.2)
  
  return(plot)
  
}
plotBarGeneArray <- function(VisiumSample, region, geneList, comparedGroups, maxN, saveFileName){
  plotlist <- list()
  count <- 1
  while(count < nrow(geneList)+1){
    gene_summary <- formatGeneExpression(VisiumSample, region, as.character(geneList[count,1]), comparedGroups, maxN)
    nacount <- colSums(is.na(gene_summary))
    if(nacount[1]!=(ncol(comparedGroups[[1]])-1)){
      plotlist[[count]] <- plotBarGene(gene_summary,as.character(geneList[count,1]), as.character(region))
    }
    count = count+1
  }
  multi.page <- ggarrange(plotlist = plotlist, nrow = 3, ncol = 3)
  ggexport(multi.page, filename = saveFileName)
}

# Load file if there are specific genes you want to plot
geneSample <- read_excel('genesOfInterest.xlsx')





#### Analysis of the data ####

#### Example analysis
# Compare gene expression values between Group1 and Group2 in Region1 with no thresholding
masterListComparison <- compareGroups(VisiumSample, "Region1", "Group1", "Group2", 1)

# Extract a list of Differentially Expressed Genes (DEGs) from the comparison above
hitDEGs <- getHitList(masterListComparison)

# Generate a volcano plot of the hit genes
plotregionvolcano(masterListComparison, "hits")

# Generate a PDF of bar graphs for hit gene expression levels with a maximum replicate number of 4
plotBarGeneArray(VisiumSample, "Region1", hitDEGs, masterListComparison, 4, "Outputs/testHITS.pdf")

# Generate a volcano plot of other genes of interest from geneSample
plotregionvolcano(masterListComparison, "goi")

# Generate a PDF of bar graphs for hit gene expression levels with a maximum replicate number of 4
plotBarGeneArray(VisiumSample, "Region1", geneSample, masterListComparison, 4, "Outputs/testGOI.pdf")

# Save a CSV file of the masterListComparison
write_csv(masterListComparison[[1]],"Outputs/testMasterListComparison.csv")
