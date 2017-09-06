library("DESeq2")
library("ggplot2")
library("pcaExplorer")
library("VennDiagram")
library("reshape")
library(VennDiagram)
library(devtools)
library(ggbiplot)
install_github("ggbiplot", "vqv")


# pathway; map to ‘htseq_count’ file.
setwd("/Users/langyiying/Dropbox_Work/Dropbox/Msc_Lectures/RP2/htseq_count")
directory <- "/Users/langyiying/Dropbox_Work/Dropbox/Msc_Lectures/RP2/htseq_count"

### initiation
files <- list.files(directory)
type <- c("vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","vvd","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt","wt")
temp <- c(25,25,25,25,25,25,25,25,25,25,28,28,28,28,28,28,28,28,28,28,25,25,25,25,25,25,25,25,25,25,28,28,28,28,28,28,28,28,28,28)
time <- c(0,0,4,4,8,8,12,12,16,16,0,0,4,4,8,8,12,12,16,16,0,0,4,4,8,8,12,12,16,16,0,0,4,4,8,8,12,12,16,16)
rep <- c("rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02","rep01","rep02")
luminance <- c("light","light","dark","dark","dark","dark","dark","dark","dark","dark","light","light","dark","dark","dark","dark","dark","dark","dark","dark","light","light","dark","dark","dark","dark","dark","dark","dark","dark","light","light","dark","dark","dark","dark","dark","dark","dark","dark")
temp <- factor(temp)
time <- factor(time)
#colData <- cbind(type,temp,time)
#rownames(colData) <- files
sampleTable <- data.frame(sampleName=files,fileName=files,type,temp,time,rep,luminance)


### All
DEdataset_all <- DESeqDataSetFromHTSeqCount(sampleTable,directory,~type)
DESeq_all <- DESeq(DEdataset_all)
res_all <- results(DESeq_all)
GeneNames <- rownames(res_all)


##This second way for creating count table is for the PCA package plotting. Dataframe row names and the like probably need changing for different datasets... (HINT)
#for this dataset, find ‘MergedHTSeq.txt’ file. It is htseq_count merged.
CountTable <- read.table("/Users/Neurospora/Desktop/Joseph_Things/MergedHTSeq.txt", header= TRUE)  #Made with custom python script for all 40 files. header important.
CountTable <- head(CountTable,-5)
row.names(CountTable) <- GeneNames
colData <- data.frame(row.names = c("VVD25DD0R1", "VVD25DD0R2", "VVD25DD4R1", "VVD25DD4R2", "VVD25DD8R1", "VVD25DD8R2", "VVD25DD12R1", "VVD25DD12R2", "VVD25DD16R1", "VVD25DD16R2",
                                    "VVD28DD0R1", "VVD28DD0R2", "VVD28DD4R1", "VVD28DD4R2", "VVD28DD8R1", "VVD28DD8R2", "VVD28DD12R1", "VVD28DD12R2", "VVD28DD16R1", "VVD28DD16R2",
                                    "WT25DD0R1", "WT25DD0R2", "WT25DD4R1", "WT25DD4R2", "WT25DD8R1", "WT25DD8R2", "WT25DD12R1", "WT25DD12R2", "WT25DD16R1", "WT25DD16R2",
                                    "WT28DD0R1", "WT28DD0R2", "WT28DD4R1", "WT28DD4R2", "WT28DD8R1", "WT28DD8R2", "WT28DD12R1", "WT28DD12R2", "WT28DD16R1", "WT28DD16R2"), 
                                    genotypes = rep(c("VVD", "WT"),  each = 20), temp  = rep(c("25", "28"),2,  each = 10), repetition = rep(c("R1", "R2"),20), 
                                    times = rep(c("DD0", "DD4", "DD8", "DD12", "DD16"), 2, each = 2) )
dds <- DESeqDataSetFromMatrix(countData = CountTable, colData = colData, design = ~ genotypes) #maybe add time/reps?

#This stops an error code of trying to get log10(0).
CountTable[CountTable == 0] <- 0.1

#Creates certain files. Always run.
dda <- DESeq(dds)
res <- results(dda)
res
plotMA(dda, ylim=c(-2,2), main = "DESeq2")
rld <- rlog(dda)

#PLot PCAs with program. in the line show all options (rep, geno, time & temp). Change to suit your needs.
plotPCA(rld, intgroup = c("repetition", "times", "genotypes", "temp"))

#This is for ggplot custom modification of plotPCA plots. makes it prettier; look! axis labels are custom, so run program first to find out true values.
StrainTemp <- rep(c("VVD25", "VVD28", "WT25", "WT28"), each = 10)
TimePoints <- rep(c("0-DD", "4-DD", "8-DD", "DD-12", "DD-16"), each = 2, 4)
TimePoints
#HOLY CODE, BATMAN! this customises plot
plottedpca <- plotPCA(rld, ylim(-50, 50), intgroup = c("genotypes", "times", "temp"), returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(plottedpca, aes(PC1, PC2, color= StrainTemp, shape=TimePoints)) +
  geom_point(size=3) +
  theme_light()+
  xlab(paste0("PC1: 60 % variance")) +
  ylab(paste0("PC2: 10 % variance"))


#Started re-creating PCA plot with x-y plot. 
GenoTemps <- rep(c('VVD25', "VVD28", "WT25", "WT28"), each = 10)
TimePoints <- rep(c("DD00", "DD04", "DD08", "DD12", "DD16"), each = 2, 4)
xy <- as.numeric(c(-18.67729, -17.65599, -13.98881, -16.33067, -7.831102, -7.647084, -4.082487, 2.221523, 16.47063, 9.0403, -11.16064, -13.44571, -8.982194, -11.27937, 1.129003, -3.481472, 8.96641, 8.360086, 24.7726, 11.16207, -13.59813, -16.78496, -7.411196, -9.15474, -6.14304, -2.44799, 4.207044, 9.798175, 16.32092, 13.93885, -5.920845, -12.8372, -1.346649, -4.698874, 9.273313, 4.103207, 16.35892, 16.61436, 20.60195, 21.56708))    # the eruption durations 
yx <- as.numeric(c(-7.626281, -7.641205, -1.492138, 3.809963, -1.528365, 6.892218, -1.6507, 3.378833, 0.261571, 5.313718, -9.251099, -7.707719, -3.543133, -0.9247928, 0.1811116, 0.1368201, -5.290745, -1.453386, -12.80915, -1.320349, -2.969427, 3.728188, -1.056797, 6.921788, 1.108382, 7.673831, 1.706793, 7.291645, 2.121573, 6.331316, -4.387491, 5.826837, 1.205871, 7.174366, 1.298804, 6.505103, -1.589629, 1.912893, -6.788624, -1.750595))         # the waiting interval 
qplot(xy, yx, colour = GenoTemps, shape = TimePoints,xlab="PC1 60% variance", ylab="PC2 10% variance")



#### this plots all Linear Regressions and saves them to Joseph_`Things on desktop
#Run from next line, expand image to reolution wanted then run 'dev.copy2pdf' line. THEN close image.
x11()
par(mfrow=c(4,5))

#This whole thing does Linear regression for each condition and it's replicate.
CountTable[CountTable == 0] <- 0.1 #replaces '0' for '0.1' in df. This is to allow Log10 to work.
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39)){
  rep1 = log10(CountTable[,i])
  rep2 = log10(CountTable[,i+1])
  StatForTitle <- lm(rep2 ~rep1)
  title <- paste(StatForTitle$coefficients[1], StatForTitle$coefficients[2], sep = ' ')
  plot(rep1, rep2, xlab = colnames(CountTable[i]),ylab = colnames(CountTable[i +1]), main = "Log Linear Regression of two repetitions")
  legend("topleft", bty="n", legend=paste("Linear Regression Line is ", signif(StatForTitle$coefficients[1],4), " + ", signif(StatForTitle$coefficients[2],4), "x", sep = ''))
  legend("bottomright", bty="n", legend=paste("R2 is ", signif(summary(StatForTitle)$r.squared, 4)))
  abline(lm(rep2~rep1))
}
dev.copy2pdf(file = "/Users/Neurospora/Desktop/Joseph_Things/LRModels.pdf")

#####

#### same as above, but title is now Line Eq + R2 value. Prints better on 4*5
x11()
par(mfrow=c(4,5))
forcsv <- c()

#This whole thing does Linear regression for each condition and it's replicate.
CountTable[CountTable == 0] <- 0.1 #replaces '0' for '0.1' in df. This is to allow Log10 to work.
for (i in c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39)){
  rep1 = log10(CountTable[,i])
  rep2 = log10(CountTable[,i+1])
  StatForTitle <- lm(rep2 ~rep1)
  title <- paste(StatForTitle$coefficients[1], StatForTitle$coefficients[2], sep = ' ', cex = 2)
  plot(rep1, rep2,  xlab = colnames(CountTable[i]),ylab = colnames(CountTable[i +1]), main = paste("Line = ", signif(StatForTitle$coefficients[1],4), " + ", signif(StatForTitle$coefficients[2],4), "x, R2 = ", signif(summary(StatForTitle)$r.squared, 4),sep = '', cex = 2))
  abline(lm(rep2~rep1))
  addtocsv <- paste(unique(strsplit(colnames(CountTable[i]), 'R')[[1]])[1], signif(summary(StatForTitle)$r.squared, 4), sep = " ")
  forcsv <- c(forcsv, addtocsv)
}

forcsv1 <- c("Sample and it's Replicate")
forcsv2 <- c("R-Squared Value for similarity")
for ( l in forcsv){
  data <- strsplit(l,' ')
  forcsv1 <- c(forcsv1, toString(data[[1]][1]))
  forcsv2 <- c(forcsv2, toString(data[[1]][2]))
}
poyyy <- data.frame(forcsv1, forcsv2)
write.csv(x = poyyy, file = "/Users/Neurospora/Desktop/Joseph_Things/TableOfRSquaredOfReplicates.csv")
dev.copy2pdf(file = "/Users/Neurospora/Desktop/Joseph_Things/LRModels.pdf")

#####

#remnants from original file.
plotCounts(DESeq_all,"NCU16528",intgroup=c("type","temp","time")) #hypothetical protein
plotCounts(DESeq_all,"NCU03967",intgroup=c("type","temp","luminance"),normalized=TRUE) #vvd
VVDcounts<-plotCounts(DESeq_all,"NCU08626",intgroup=c("type","temp", "time"),normalized=TRUE,returnData = TRUE) #vvd
gg<-plotCounts(DESeq_all,"NCU16528",intgroup=c("type","temp","luminance"),normalized=TRUE)
gg <- plotCounts(DESeq_all,"NCU08626",intgroup=c("type","temp","luminance", "time"),normalized=TRUE,returnData = TRUE) #hypothetical protein

#this tests (by placing in diff NCU numbers) count data for genes upreg. at diff temps
sourcefile <- read.table("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_temp25to28_DD16", header = TRUE, row.names = 1)
for (i in rownames(sourcefile)){
  testgene <- plotCounts(DESeq_all, i,intgroup=c("type","temp", "time"),normalized=TRUE,returnData = TRUE)
  counts <- testgene$count
  barplot(counts, main = paste("Counts for test gene", i, sep = ''), xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)  
  }


#Input location of DESeq output file to be analysed. Greps name, then uses 'barmaker' to create barplots for expression
tempgenes <- read.csv("/Users/Neurospora/PycharmProjects/HistogramMaker/UpstreamFound_vvd25_0_4_names.txt", header = FALSE, sep = "\"")
barmaker(tempgenes$V1) 

for (i in secondfilelist){ #This checks DE genes and prints barplots for total counts and gene presence per timepoint
  singletonbars(i)
}

#Test here any Accession number and it'll plot you some pretty plots.
singletonbars("NCU03967")
templist <- c("NCU11424", "NCU00582", "NCU08626")
barmaker(templist)


#calculation for if genes are up or down reg....
PercentAbove0("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_DD0_2528.txt")
PercentAbove0("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_DD4_2528.txt")
PercentAbove0("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_DD8_2528.txt")
PercentAbove0("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_DD12_2528.txt")
PercentAbove0("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_DD16_2528.txt")

#find genes shared in DE b/w 25 and 28C with WT and VVD. Then, GO later.
A1 <- namegrabber(("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD0_28.txt"))
B1 <- namegrabber(("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD4_28.txt"))
C1 <- namegrabber(("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD8_28.txt"))
D1 <- namegrabber(("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD12_28.txt"))
E1 <- namegrabber(("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD16_28.txt"))

A2 <- namegrabber(("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD0_25.txt"))
B2 <- namegrabber(("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD4_25.txt"))
C2 <- namegrabber(("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD8_25.txt"))
D2 <- namegrabber(("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD12_25.txt"))
E2 <- namegrabber(("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD16_25.txt"))

AA <- setdiff(A1, A2) #it's LL. ignore.
AA
BB <- setdiff(B1, B2)
BB
CC <- setdiff(C1, C2)
CC
DD <- setdiff(D1, C2)
DD
EE <- setdiff(E1, E2) #only VVD here. ignore it.
EE

BBCC <- intersect(BB,CC)
BBCC
BBCCDD <- intersect(BBCC,DD)
BBCCDD

#Making github CSVs

  for (times in c(0,4,8,12,16)){
    MakingCSVs(paste("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_VVD_vs_WT_temp25_DD",times,".txt", sep = '' ),paste("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_VVD_vs_WT_temp28_DD",times,".txt", sep = '' ), paste("VVD_DD", times, "_25C", sep = ''), paste("WT_DD", times,"_28C", sep = ''), paste("/Users/Neurospora/Desktop/Joseph_Things/FilesForGitHub/","VVD_vs_WT_DD", times, "for_25 vs_28.csv", sep=''))
  }


MakingCSVs <- function(file1, file2, name1, name2, outputnameanddirectory){
forfirst <- (read.csv(file1 , sep = "\""))$baseMean
forsecond <- (read.csv(file2, sep = "\""))$baseMean
firstfile <- read.csv(file1 , sep = " ")
secondfile <-read.csv(file2 , sep = " ")
both <- intersect(forfirst, forsecond)
OnlyForFirst <- setdiff(forfirst, forsecond)
OnlyForSecond <- setdiff(forsecond, forfirst)
Listings <- paste("GeneName", name1, name2, sep = ' ')
for (u in both){
  if (as.double(firstfile[u, "log2FoldChange"]) > 0 && as.double(secondfile[u, "log2FoldChange"]) > 0){
    Listings <- c(Listings, paste(u, "Upregulated", "Upregulated",sep = ' '))
  }
  else if (as.double(firstfile[u, "log2FoldChange"]) < 0 && as.double(secondfile[u, "log2FoldChange"]) > 0){
    Listings <- c(Listings, paste(u, "Downregulated", "Upregulated",sep = ' '))
  }
  else if (as.double(firstfile[u, "log2FoldChange"]) > 0 && as.double(secondfile[u, "log2FoldChange"]) < 0){
    Listings <- c(Listings, paste(u, "Upregulated", "Downregulated",sep = ' '))
  }
  else {
    Listings <- c(Listings, paste(u, "Downregulated", "Downregulated",sep = ' '))
  }
}
for (u in OnlyForFirst){
  if (as.double(firstfile[u, "log2FoldChange"]) > 0){
    Listings <- c(Listings, paste(u, "Upregulated", "---",sep = ' '))
  }
  else{
    Listings <- c(Listings, paste(u, "Downregulated", "---",sep = ' '))
  }
}
for (u in OnlyForSecond){
  if (as.double(secondfile[u, "log2FoldChange"]) > 0){
    Listings <- c(Listings, paste(u, "---", "Upregulated",sep = ' '))
  }
  else {
    Listings <- c(Listings, paste(u, "---", "Downregulated",sep = ' '))
  }
}
Listings1 <- c()
Listings2 <- c()
Listings3 <- c()
for (i in Listings){
Listings1 <- c(Listings1, toString(strsplit(i, ' ')[[1]][1]))
Listings2 <- c(Listings2, toString(strsplit(i, ' ')[[1]][2]))
Listings3 <- c(Listings3, toString(strsplit(i, ' ')[[1]][3]))
}
piyyy <- data.frame(Listings1, Listings2, Listings3)
write.csv(x = piyyy, file = outputnameanddirectory)
}


#setdiif(x,y) finds only-X.

barmaker("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_wt_DD16_28.txt")


#Run this for finding similarity, make bars, etc.
firstbit <- "/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_TESTYvvd.txt"
secondbit<- "/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_TESTYwt.txt"
firstfilelist <- namegrabber(firstbit)
secondfilelist <- namegrabber(secondbit)
firsttitle <- strsplit(toString(sapply(strsplit(firstbit,"/"), `[`, 7)), ".txt")
secondtitle <- strsplit(toString(sapply(strsplit(secondbit,"/"), `[`, 7)), ".txt")
similaritycheck(as.character(firstfilelist), as.character(secondfilelist)) #two-way Venn diagram creation
potate <- tablemaker(firstfilelist, secondfilelist) #Don't do these two next lines unless you want output csv
write.csv(file = "/Users/Neurospora/Desktop/Joseph_Things/OutputA5.csv", x = potate)



#More tests!
listofgenes <- c('NCU08907', 'NCU07787', 'NCU00399', 'NCU07817', 'NCU05395', 'NCU08909', 'NCU07569')
barmaker(listofgenes)
for (i in listofgenes){
  singletonbars("NCU05259")
}


#This creates PCA... Old code.
PCADataPrep(firstfilelist)


#plots all expression values on a single plot. ugly & didnt help, so I didnt adjust it.
MultipleLine(secondfilelist, 50000)



MultipleLine <- function(GeneNames, lowerend){
  #Initialise blank vectors of 40 sets of 0.

  cat("Final Countdown\n") #Its funny.
  HighlyExpressedGenes <- NULL
  numba <- NROW(GeneNames)  #Number of genes involved in the RNASeq data\
  MAXIMAL <- 10
  for (l in GeneNames){ #Genenames are the names of the genes, extracted previously. Loop through, extract count data, add to main plot.
    wholeplotcount <- plotCounts(DESeq_all, l,intgroup=c("type","temp", "time"),normalized=TRUE,returnData = TRUE)
    singlecountdata <- wholeplotcount$count 
    for (u in (1:40)){
      if (as.double(singlecountdata[u]) >= as.double(MAXIMAL)){
        MAXIMAL <- singlecountdata[u]
      }
      
    }
  }
  
  plot(1, type="n", main = "DE Gene counts plotted together", xlab="20 Conditions", ylab="Counts", xlim = c(0,20), ylim = c(0, MAXIMAL))
  AllNCUs <- ''
  for (l in GeneNames){ #Genenames are the names of the genes, extracted previously. Loop through, extract count data, add to main plot.
    numba <- numba -1
    cat(numba)
    cat("...") #This little bit is just so I know progress of loop... though it slows it all down alot...
    wholeplotcount <- plotCounts(DESeq_all, l,intgroup=c("type","temp", "time"),normalized=TRUE,returnData = TRUE)
    singlecountdata <- wholeplotcount$count 
    SC1 <- (singlecountdata[1] + singlecountdata[2])/2
    SC2 <- (singlecountdata[3] + singlecountdata[4])/2
    SC3 <- (singlecountdata[5] + singlecountdata[6])/2
    SC4 <- (singlecountdata[7] + singlecountdata[8])/2
    SC5 <- (singlecountdata[9] + singlecountdata[10])/2
    SC6 <- (singlecountdata[11] + singlecountdata[12])/2
    SC7 <- (singlecountdata[13] + singlecountdata[14])/2
    SC8 <- (singlecountdata[15] + singlecountdata[16])/2
    SC9 <- (singlecountdata[17] + singlecountdata[18])/2
    SC10 <- (singlecountdata[19] + singlecountdata[20])/2
    SC11 <- (singlecountdata[21] + singlecountdata[22])/2
    SC12 <- (singlecountdata[23] + singlecountdata[24])/2
    SC13 <- (singlecountdata[25] + singlecountdata[26])/2
    SC14 <- (singlecountdata[27] + singlecountdata[28])/2
    SC15 <- (singlecountdata[29] + singlecountdata[30])/2
    SC16 <- (singlecountdata[31] + singlecountdata[32])/2
    SC17 <- (singlecountdata[33] + singlecountdata[34])/2
    SC18 <- (singlecountdata[35] + singlecountdata[36])/2
    SC19 <- (singlecountdata[37] + singlecountdata[38])/2
    SC20 <- (singlecountdata[39] + singlecountdata[40])/2
    yt <- c(SC1, SC2, SC3, SC4, SC5, SC6, SC7, SC8, SC9, SC10, SC11, SC12, SC13, SC14, SC15, SC16, SC17, SC18, SC19, SC20 )
    if(as.double(SC5) > as.double(lowerend)){
      lines(yt,col="black",lty=1,lwd=1, pch = 20)
      AllNCUs <- c(AllNCUs, l)
    }
  }
  cat("Genes on graph are: ", as.character(AllNCUs))
  }
  



###THIS IS PCA STUFF; input gene list as vector of strings.... all useless now.
RunPCA <- function(genelisty){
testy <- PCADataPrep(genelisty)
#bit of summation data.

sumtesty1 <- rep(0, 40)
for (u in c(1:40)){
  sumtesty1[u] <- sum(testy[u+1])
}
#1 = VVDR1, 2 = VVDR2, 3= WTR1, 4= WTR2
DD0.25 <- c(sumtesty1[c(1,2,21,22)])
DD4.25 <- c(sumtesty1[c(3,4,23,24)])
DD8.25 <- c(sumtesty1[c(5,6,25,26)])
DD12.25 <- c(sumtesty1[c(7,8,27,28)])
DD16.25 <- c(sumtesty1[c(9,10,29,30)])
DD0.28 <- c(sumtesty1[c(11,12,31,32)])
DD4.28 <- c(sumtesty1[c(13,14,33,34)])
DD8.28 <- c(sumtesty1[c(15,16,35,36)])
DD12.28 <- c(sumtesty1[c(17,18,37,38)])
DD16.28 <- c(sumtesty1[c(19,20,39,40)])
Names = c("VVDR1", "VVDR2", "WTR1", "WTR2")
sumtesty <- data.frame(Names = Names,DD0.25 = DD0.25, DD4.25 = DD4.25, DD8.25 = DD8.25, DD12.25 = DD12.25, DD16.25 = DD16.25, DD0.28 = DD0.28, DD4.28 = DD4.28, DD8.28 = DD8.28, DD12.28 = DD12.28, DD16.28 = DD16.28)

log.ir <- log(testy[, 2:41])   #this is for all columns, all genes, not summed, 40 conditions
#log.ir <- log(testy[, c(3,4,23,24)]) #same, but only DD0 and 25C.
#log.ir <- log(testy[, c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)])
log.ir2 <- log(sumtesty[, 2:11])
# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
ir.pca <- prcomp(log.ir,
                 center = TRUE,
                 scale. = TRUE) 
# advisable, but default is FALSE. 
ir.pca2 <- prcomp(log.ir2,
                 center = TRUE,
                 scale. = TRUE) 
print(ir.pca)
print(ir.pca2)
g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
  ellipse = TRUE)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal', labels = testy$GeneName,
               legend.position = 'top')
print(g)

g6 <- ggbiplot(ir.pca2, obs.scale = 1, var.scale = 1, labels = sumtesty$Names, 
              ellipse = TRUE)
g6 <- g6 + scale_color_discrete(name = '')
g6 <- g6 + theme(legend.direction = 'horizontal', 
               legend.position = 'top')
print(g6)
}
###This is PCA STUFF


PCADataPrep <- function(genelist){
  vvd25DD0R1 <- rep(0, length(genelist))
  vvd25DD0R2 <- rep(0, length(genelist))
  vvd25DD4R1 <- rep(0, length(genelist))
  vvd25DD4R2 <- rep(0, length(genelist))
  vvd25DD8R1 <- rep(0, length(genelist))
  vvd25DD8R2 <- rep(0, length(genelist))
  vvd25DD12R1 <- rep(0, length(genelist))
  vvd25DD12R2 <- rep(0, length(genelist))
  vvd25DD16R1 <- rep(0, length(genelist))
  vvd25DD16R2 <- rep(0, length(genelist))
  vvd28DD0R1 <- rep(0, length(genelist))
  vvd28DD0R2 <- rep(0, length(genelist))
  vvd28DD4R1 <- rep(0, length(genelist))
  vvd28DD4R2 <- rep(0, length(genelist))
  vvd28DD8R1 <- rep(0, length(genelist))
  vvd28DD8R2 <- rep(0, length(genelist))
  vvd28DD12R1 <- rep(0, length(genelist))
  vvd28DD12R2 <- rep(0, length(genelist))
  vvd28DD16R1 <- rep(0, length(genelist))
  vvd28DD16R2 <- rep(0, length(genelist))
  wt25DD0R1 <- rep(0, length(genelist))
  wt25DD0R2 <- rep(0, length(genelist))
  wt25DD4R1 <- rep(0, length(genelist))
  wt25DD4R2 <- rep(0, length(genelist))
  wt25DD8R1 <- rep(0, length(genelist))
  wt25DD8R2 <- rep(0, length(genelist))
  wt25DD12R1 <- rep(0, length(genelist))
  wt25DD12R2 <- rep(0, length(genelist))
  wt25DD16R1 <- rep(0, length(genelist))
  wt25DD16R2 <- rep(0, length(genelist))
  wt28DD0R1 <- rep(0, length(genelist))
  wt28DD0R2 <- rep(0, length(genelist))
  wt28DD4R1 <- rep(0, length(genelist))
  wt28DD4R2 <- rep(0, length(genelist))
  wt28DD8R1 <- rep(0, length(genelist))
  wt28DD8R2 <- rep(0, length(genelist))
  wt28DD12R1 <- rep(0, length(genelist))
  wt28DD12R2 <- rep(0, length(genelist))
  wt28DD16R1 <- rep(0, length(genelist))
  wt28DD16R2 <- rep(0, length(genelist))

  for (i in c(1: length(genelist))){
    tempcounts <- plotCounts(DESeq_all, as.character(genelist[i]), intgroup=c("type","temp", "time"),normalized=TRUE,returnData = TRUE)
    CountsOnly <- tempcounts$count
    vvd25DD0R1[i] <- CountsOnly[1]
    vvd25DD0R2[i] <- CountsOnly[2]
    vvd25DD4R1[i] <- CountsOnly[3]
    vvd25DD4R2[i] <- CountsOnly[4]
    vvd25DD8R1[i] <- CountsOnly[5]
    vvd25DD8R2[i] <- CountsOnly[6]
    vvd25DD12R1[i] <- CountsOnly[7]
    vvd25DD12R2[i] <- CountsOnly[8]
    vvd25DD16R1[i] <- CountsOnly[9]
    vvd25DD16R2[i] <- CountsOnly[10]
    vvd28DD0R1[i] <- CountsOnly[11]
    vvd28DD0R2[i] <- CountsOnly[12]
    vvd28DD4R1[i] <- CountsOnly[13]
    vvd28DD4R2[i] <- CountsOnly[14]
    vvd28DD8R1[i] <- CountsOnly[15]
    vvd28DD8R2[i] <- CountsOnly[16]
    vvd28DD12R1[i] <- CountsOnly[17]
    vvd28DD12R2[i] <- CountsOnly[18]
    vvd28DD16R1[i] <- CountsOnly[19]
    vvd28DD16R2[i] <- CountsOnly[20]
    wt25DD0R1[i] <- CountsOnly[21]
    wt25DD0R2[i] <- CountsOnly[22]
    wt25DD4R1[i] <- CountsOnly[23]
    wt25DD4R2[i] <- CountsOnly[24]
    wt25DD8R1[i] <- CountsOnly[25]
    wt25DD8R2[i] <- CountsOnly[26]
    wt25DD12R1[i] <- CountsOnly[27]
    wt25DD12R2[i] <- CountsOnly[28]
    wt25DD16R1[i] <- CountsOnly[29]
    wt25DD16R2[i] <- CountsOnly[30]
    wt28DD0R1[i] <- CountsOnly[31]
    wt28DD0R2[i] <- CountsOnly[32]
    wt28DD4R1[i] <- CountsOnly[33]
    wt28DD4R2[i] <- CountsOnly[34]
    wt28DD8R1[i] <- CountsOnly[35]
    wt28DD8R2[i] <- CountsOnly[36]
    wt28DD12R1[i] <- CountsOnly[37]
    wt28DD12R2[i] <- CountsOnly[28]
    wt28DD16R1[i] <- CountsOnly[39]
    wt28DD16R2[i] <- CountsOnly[40]
  }
  PCATable <- data.frame(GeneName = genelist, vvd25DD0R1 = vvd25DD0R1, vvd25DD0R2 = vvd25DD0R2, vvd25DD4R1 = vvd25DD4R1, vvd25DD4R2 = vvd25DD4R2, vvd25DD8R1 = vvd25DD8R1, vvd25DD8R2 = vvd25DD8R2, vvd25DD12R1 = vvd25DD12R1, vvd25DD12R2 = vvd25DD12R2,
                         vvd25DD16R1 = vvd25DD16R1, vvd25DD16R2 = vvd25DD16R2, vvd28DD0R1 = vvd28DD0R1, vvd28DD0R2 = vvd28DD0R2, vvd28DD4R1 = vvd28DD4R1, vvd28DD4R2 = vvd28DD4R2, vvd28DD8R1 = vvd28DD8R1, vvd28DD8R2 = vvd28DD8R2, 
                         vvd28DD12R1 = vvd28DD12R1, vvd28DD12R2 = vvd28DD12R2, vvd28DD16R1 = vvd28DD16R1, vvd28DD16R2 = vvd28DD16R2,  wt25DD0R1 = wt25DD0R1, wt25DD0R2 = wt25DD0R2, wt25DD4R1 = wt25DD4R1, wt25DD4R2 = wt25DD4R2,
                         wt25DD8R1 = wt25DD8R1, wt25DD8R2 = wt25DD8R2, wt25DD12R1 = wt25DD12R1, wt25DD12R2 = wt25DD12R2, wt25DD16R1 = wt25DD16R1, wt25DD16R2 = wt25DD16R2, wt28DD0R1 = wt28DD0R1, wt28DD0R2 = wt28DD0R2, wt28DD4R1 = wt28DD4R1, wt28DD4R2 = wt28DD4R2,
                         wt28DD8R1 = wt28DD8R1, wt28DD8R2 = wt28DD8R2, wt28DD12R1 = wt28DD12R1, wt28DD12R2 = wt28DD12R2, wt28DD16R1 = wt28DD16R1, wt28DD16R2 = wt28DD16R2)
return(PCATable)
}


# This will only check for DD0, as that is main focus of my write-up. But can be modified for any stat test (singlecountdata[x]...)
StatTester <- function(A){
  
  wholeplotcount <- plotCounts(DESeq_all, A ,intgroup=c("type","temp", "time"),normalized=TRUE,returnData = TRUE)
  singlecountdata <- wholeplotcount$count
  x <- c(singlecountdata[1],singlecountdata[2], singlecountdata[11], singlecountdata[12])
  y <- c(singlecountdata[21],singlecountdata[22], singlecountdata[31], singlecountdata[32])
  FTest <- var.test(x, y, ratio = 1,
           alternative = c("two.sided", "less", "greater"),
           conf.level = 0.95)
  FTest$p.value
  if ((FTest$p.value) > 0.5){
    return(t.test(x, y, var.equal=TRUE))
  }
  else {
    return(t.test(x, y, var.equal=FALSE))
  }

}

##grab gene names of DEseq file with this, creates list
namegrabber <- function(A){
  filesource <- read.csv(A, sep = "\"")
  namelist <- filesource$baseMean
  return(namelist)
}

#Number of genes, as percentage, in file with expression above 0 (=> higher in WT). Also prints gene total.
PercentAbove0 <- function(fileA){
  FileAReader <- read.csv(fileA, sep = ' ')
  Over0 <- sum(FileAReader$log2FoldChange >= 0)
  Percentage <- (Over0/NROW(FileAReader))*100
  TotalGenes <- NROW(FileAReader)
  return(paste("Percentage upregulated:", Percentage," total genes: ", TotalGenes))
}

#Make table of intersections
tablemaker <- function(A,B){
  intersection <- intersect(A,B)
  OnlyA <- setdiff(A,B)
  OnlyB <- setdiff(B,A)
  totalvector <- c(intersection, OnlyA, OnlyB)
  vector1 = rep(0, (length(intersection)) + (length(OnlyA)) + (length(OnlyB)))
  vector2 = rep(0, (length(intersection)) + (length(OnlyA)) + (length(OnlyB)))
  vector3 = rep(0, (length(intersection)) + (length(OnlyA)) + (length(OnlyB)))
  
  for (i in c(1:(length(intersection)))){
    vector1[i] = totalvector[i]
    vector2[i] = totalvector[i]
    vector3[i] = totalvector[i]
  }
  for (a in c((length(intersection)+1):(length(intersection) + length(OnlyA)))){
    vector1[a] = totalvector[a]
    vector2[a] = 'NULL'
    vector3[a] = 'NULL'
  }
  for (b in c((length(intersection) + length(OnlyA) +1):length(totalvector))){
    vector1[b] = 'NULL'
    vector2[b] = 'NULL'
    vector3[b] = totalvector[b]
  }
  return(data.frame(Vector1 = vector1, Vector2 = vector2, Vector3 = vector3))
}


##Similarity check of two condition's gene DE. creates Venn diagram. it's pretty.
similaritycheck <- function(A,B){
  intersection <- intersect(A,B)
  OnlyA <- setdiff(A,B)
  OnlyB <- setdiff(B,A)
  grid.newpage()
  draw.pairwise.venn(margin = 0.08, area1 = length(intersection)+ length(OnlyA), area2 = length(intersection)+ length(OnlyB), cross.area = length(intersection), category = c(firsttitle, secondtitle), fill = c("blue", "red"))
  return(cat("Intersect is: ",toString(length(intersection)), "\nIntersect dataset is: ", intersection,  "\nFirst dataset-only are: ", toString(length(OnlyA)) , "\nFirst dataset uniques are: ", toString(OnlyA), "\nSecond dataset-only are: ", toString(length(OnlyB)), "\nSecond dataset uniques are: ", toString(OnlyB)))
}

#main script; creates tables of gene expression/count numebers. Count= RNA count. Gene numbers= 'if gene expressed at this time-point'
barmaker <- function(GeneNames){
#Initialise blank vectors of 40 sets of 0.
countdata <- rep(0, 40)
countdata2 <- rep(0, 40)
countdata3 <- rep(0, 40)
countdata4 <- rep(0, 40)
countdata5 <- rep(0, 40)
countdata6 <- rep(0, 40)
countdata7 <- rep(0, 40)
countdata8 <- rep(0, 40)
countdata10 <- rep(0, 40)
countdata12 <- rep(0, 40)
countdata13 <- rep(0, 40)
countdata14 <- rep(0, 40)
countdata15 <- rep(0, 40)
countdata16 <- rep(0, 40)
countdata17 <- rep(0, 40)
countdata18 <- rep(0, 40)
cat("Final Countdown") #Its funny.
HighlyExpressedGenes <- NULL
numba <- NROW(GeneNames)  #Number of genes involved in the RNASeq data
for (l in GeneNames){ #Genenames are the names of the genes, extracted previously. Loop through, extract count data, add to main plot.
  numba <- numba -1
  cat(numba)
  cat("...") #This little bit is just so I know progress of loop... though it slows it all down alot...
  wholeplotcount <- plotCounts(DESeq_all, l,intgroup=c("type","temp", "time"),normalized=TRUE,returnData = TRUE)
  singlecountdata <- wholeplotcount$count 
  for (i in c(1:40)){
        countdata[i] <- countdata[i] + singlecountdata[i]       #Total RNA count for genes in each exp. point
    if (singlecountdata[i] > 10) {
          countdata10[i] <- countdata10[i] + 1                    #Total Gene count in each exp. point (>10 count = expressed.)
        }
    if (singlecountdata[i] > 0 & singlecountdata[i] < 10){
        countdata2[i] <- countdata2[i] + singlecountdata[i] 
        countdata12[i] <- countdata12[i] + 1
    }
    else if (singlecountdata[i] >= 10 & singlecountdata[i] < 100){
      countdata3[i] <- countdata3[i] + singlecountdata[i]   #find three thresholds. count >0, >5, >10 . see then if peaks differ. else, higher diff.
      countdata13[i] <- countdata13[i] + 1
      }
    else if (singlecountdata[i] >=100 & singlecountdata[i] < 5000){
        countdata4[i] <- countdata4[i] + singlecountdata[i] 
        countdata14[i] <- countdata14[i] + 1
    }
    else if (singlecountdata[i] >= 5000 & singlecountdata[i] < 10000){
      countdata5[i] <- countdata5[i] + singlecountdata[i]   
      countdata15[i] <- countdata15[i] + 1
    }
    else if (singlecountdata[i] >= 10000 & singlecountdata[i] < 60000){
      countdata6[i] <- countdata6[i] + singlecountdata[i]   
      countdata16[i] <- countdata16[i] + 1
      HighLength <- length(HighlyExpressedGenes)
      HighlyExpressedGenes[HighLength+1] <- l
            }
    else if (singlecountdata[i] >= 60000 & singlecountdata[i] <100000){
      countdata7[i] <- countdata7[i] + singlecountdata[i]  
      countdata17[i] <- countdata17[i] + 1
      }
    else if (singlecountdata[i] >=100000){
      countdata8[i] <- countdata8[i] + singlecountdata[i]   
      countdata18[i] <- countdata18[i] + 1
    }
        
        }
}


#This part is for creating barplots (histograms, but with counts already done!) from count data for a particular gene.
##Could also maybe do for diff. expression..?
barnames <- NULL
colours <- c(rep("red", 10), rep("blue", 10), rep("mistyrose", 10), rep("lavender",10))
for (jake in c(1:40)){
  barnames[jake] <- paste(toString(VVDcounts[jake,2]),toString(VVDcounts[jake,3]),toString(VVDcounts[jake,4]),  collapse = '', sep = '')
}


#No paste loop thing due to... reasons. This will create bar plots.
barplot(countdata, main = "Counts for All Expresed Genes", xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)
barplot(countdata2, main = "Counts for All Expressed Genes, each with counts of 0 to 10 ", xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)
barplot(countdata3, main = "Counts for All Expressed Genes, each with counts of 10 to 100", xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)
barplot(countdata4, main = "Counts for All Expressed Genes, each with counts of 100 to 5000", xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)
barplot(countdata5, main = "Counts for All Expressed Genes, each with counts of 5000 to 10,000", xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)
barplot(countdata6, main = "Counts for All Expressed Genes, each with counts of 10,000 to 60,000", xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)
barplot(countdata7, main = "Counts for All Expressed Genes, each with counts of 60,000 t0 100,000", xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)
barplot(countdata8, main = "Counts for All Expressed Genes, each with counts of > 100,000", xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)
barplot(countdata10, main = "Number Of Expressed Genes", xlab = "Experiment Sample", ylab = "Number of Genes", col = colours, names.arg = barnames)
barplot(countdata12, main = "Number Of Expressed Genes, each with counts of 0 to 10", xlab = "Experiment Sample", ylab = "Number of Genes", col = colours, names.arg = barnames)
barplot(countdata13, main = "Number Of Expressed Genes, each with counts of 10 to 100", xlab = "Experiment Sample", ylab = "Number of Genes", col = colours, names.arg = barnames)
barplot(countdata14, main = "Number Of Expressed Genes, each with counts of 100 to 5000", xlab = "Experiment Sample", ylab = "Number of Genes", col = colours, names.arg = barnames)
barplot(countdata15, main = "Number Of Expressed Genes, each with counts of 5000 to 10,000", xlab = "Experiment Sample", ylab = "Number of Genes", col = colours, names.arg = barnames)
barplot(countdata16, main = "Number Of Expressed Genes, each with counts of 10,000 to 60,000", xlab = "Experiment Sample", ylab = "Number of Genes", col = colours, names.arg = barnames)
barplot(countdata17, main = "Number Of Expressed Genes, each with counts of 60,000 to 100,000", xlab = "Experiment Sample", ylab = "Number of Genes", col = colours, names.arg = barnames)
barplot(countdata18, main = "Number Of Expressed Genes, each with counts of > 100,000", xlab = "Experiment Sample", ylab = "Number of Genes", col = colours, names.arg = barnames)
}

########Next produces two bar graphs for each gene in DEseq file; first, gene count accross organism, other is '1' if gene at timepoint.
singletonbars <- function(A){
  countdata = rep(0:40)
  countdata10 = rep(0, 40)
  wholeplotcount <- plotCounts(DESeq_all, A,intgroup=c("type","temp", "time"),normalized=TRUE,returnData = TRUE)
  singlecountdata <- wholeplotcount$count
  for (i in c(1:40)){
  if (singlecountdata[i] > 10) {
    countdata10[i] <- 1                   #Total Gene count in each exp. point (>10 count = expressed.)
  }}
  
  barplot(singlecountdata, main = paste("RNA count at each timepoint for", A), xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)
  barplot(countdata10, main = paste("Gene expression at timepoints for", A), xlab = "Experiment Sample", ylab = "Counts", col = colours, names.arg = barnames)
  
  }
######


####    Try sort this to create PCA plot
sampleTable2 <- sampleTable
sampleTable2$countdata <- countdata
analysis.train <- sampleTable2[, c(4,5,8)]
analysis.species <- sampleTable2[, 1]
analysis.pca <- prcomp( analysis.train, center = TRUE, scale. = TRUE)
biplot(analysis.pca, scale = 0)


pcs <- prcomp(sampleTable2, center = T, scale = T)
prinComp <- cbind(sampleTable2, pcs$rotation)
plot(prinComp[, c("sampleName", "fileName", "type", "temp", "time", "rep", "luminance","countdata", "PC1", "PC2", "PC3")], pch = 19, cex = 0.8)
Group <- sampleTable2$type
ggplot(princomp, aes(x = temp, y = time, colour = Group)) + geom_point()
####

### Condition: compare vvd DD0 25 &28C;n then same for wt
sample_25_light <- subset(sampleTable,  !sampleTable$time == 4 & !sampleTable$time == 8 & !sampleTable$time == 12 & !sampleTable$time == 16 & sampleTable$type == 'wt')
DEdataset_25_light <- DESeqDataSetFromHTSeqCount(sample_25_light,directory, ~temp)
DESeq_25_light <- DESeq(DEdataset_25_light)
res_25_light <- results(DESeq_25_light)
summary(res_25_light)
DEgenes_25_light <- subset(res_25_light,res_25_light$padj<0.01)
DEgenes_25_light<-subset(DEgenes_25_light,abs(DEgenes_25_light$log2FoldChange)>log2(1.5))
DEgenes_25_light <- DEgenes_25_light[order(DEgenes_25_light$log2FoldChange,decreasing=TRUE),]
#DEgenes_25_light_up <- subset(DEgenes_25_light,DEgenes_25_light$log2FoldChange>0)
#DEgenes_25_light_down <- subset(DEgenes_25_light,DEgenes_25_light$log2FoldChange<0)
#write.table(DEgenes_25_light_up,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_25_light_up.txt")
#write.table(DEgenes_25_light_down,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_25_light_down.txt")
write.table(DEgenes_25_light,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_DD0_2528.txt")


### stuff
sample_25_light <- subset(sampleTable, sampleTable$type == "vvd")
DEdataset_25_light <- DESeqDataSetFromHTSeqCount(sample_25_light,directory, ~temp)
DESeq_25_light <- DESeq(DEdataset_25_light)
res_25_light <- results(DESeq_25_light)
summary(res_25_light)
DEgenes_25_light <- subset(res_25_light,res_25_light$padj<0.01)
DEgenes_25_light<-subset(DEgenes_25_light,abs(DEgenes_25_light$log2FoldChange)>log2(1.5))
DEgenes_25_light <- DEgenes_25_light[order(DEgenes_25_light$log2FoldChange,decreasing=TRUE),]
#DEgenes_25_light_up <- subset(DEgenes_25_light,DEgenes_25_light$log2FoldChange>0)
#DEgenes_25_light_down <- subset(DEgenes_25_light,DEgenes_25_light$log2FoldChange<0)
#write.table(DEgenes_25_light_up,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_25_light_up.txt")
#write.table(DEgenes_25_light_down,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_25_light_down.txt")
write.table(DEgenes_25_light,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_TESTYwt.txt")

#I'm so sorry to vandalise your beautiful code.


#wt, comparing timepoints at diff temperatures (25 vs 28)
testah<- DESeq(subset(sampleTable, sampleTable$type == "wt" & sampleTable$time == "0" & sampleTable$temp == "28" & sampleTable$rep == "rep01"))
#compares each strain at each timepoint for different temperatures.
for (typez in c("wt", "vvd")){
for (times in c(0,4,8,12,16)){
  sample_temp_2528 <- subset(sampleTable, sampleTable$type == typez & sampleTable$time== times)
  DEdataset_temp_2528 <- DESeqDataSetFromHTSeqCount(sample_temp_2528, directory, ~temp)
  DESeq_temp_2528 <- DESeq(DEdataset_temp_2528)
  res_temp_2528 <- results(DESeq_temp_2528)
  DEgenes_temp_2528 <- subset(res_temp_2528, res_temp_2528$padj <=0.01)
  DEgenes_temp_2528 <- subset(DEgenes_temp_2528, abs(DEgenes_temp_2528$log2FoldChange) > log2(1.5))
  DEgenes_temp_2528 <- DEgenes_temp_2528[order(DEgenes_temp_2528$log2FoldChange, decreasing = TRUE), ]
  write.table(DEgenes_temp_2528, paste("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_", typez,"_temp25to28_DD", times, sep = '' ))
  
}}


  sample_temp_2528 <- subset(sampleTable, sampleTable$temp == 25)
  DEdataset_temp_2528 <- DESeqDataSetFromHTSeqCount(sample_temp_2528, directory, ~type)
  DESeq_temp_2528 <- DESeq(DEdataset_temp_2528)
  res_temp_2528 <- results(DESeq_temp_2528)
  DEgenes_temp_2528 <- subset(res_temp_2528, res_temp_2528$padj <=0.01)
  DEgenes_temp_2528 <- subset(DEgenes_temp_2528, abs(DEgenes_temp_2528$log2FoldChange) > log2(1.5))
  DEgenes_temp_2528 <- DEgenes_temp_2528[order(DEgenes_temp_2528$log2FoldChange, decreasing = TRUE), ]
  write.table(DEgenes_temp_2528, paste("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvdTOwt25", sep = '' ))
  


#comparing DD0 and DD4, for both temps. So then to compare rise vs VVD and WT (both sides, then see which up and which down reg.)
for (temps in c(25,28)){
for (types in c("vvd", "wt")){

  sample_temp_2528 <- subset(sampleTable, sampleTable$temp== temps & sampleTable$type== types & !sampleTable$time== 0 & !sampleTable$time== 4 & !sampleTable$time== 8 )
  DEdataset_temp_2528 <- DESeqDataSetFromHTSeqCount(sample_temp_2528, directory, ~time)
  DESeq_temp_2528 <- DESeq(DEdataset_temp_2528)
  res_temp_2528 <- results(DESeq_temp_2528)
  DEgenes_temp_2528 <- subset(res_temp_2528, res_temp_2528$padj <=0.01)
#  DEgenes_temp_2528 <- subset(DEgenes_temp_2528, abs(DEgenes_temp_2528$log2FoldChange) > log2(1.5))  ##This part is removed for now; see upreg AND downreg genes. Python will then sort.
  DEgenes_temp_2528 <- DEgenes_temp_2528[order(DEgenes_temp_2528$log2FoldChange, decreasing = TRUE), ]
  write.table(DEgenes_temp_2528, paste("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/TESTEED", types, temps,"DD12-DD16", sep = "-" ))
}}

#VVD vs WT for 25 & 28 degrees, all timepoints
for (tempz in c (25, 28)){
  for (times in c(0,4,8,12,16)){
  sample_temp_2528 <- subset(sampleTable, sampleTable$temp == tempz & sampleTable$time== times)
  DEdataset_temp_2528 <- DESeqDataSetFromHTSeqCount(sample_temp_2528, directory, ~type)
  DESeq_temp_2528 <- DESeq(DEdataset_temp_2528)
  res_temp_2528 <- results(DESeq_temp_2528)
  DEgenes_temp_2528 <- subset(res_temp_2528, res_temp_2528$padj <=0.01)
  DEgenes_temp_2528 <- subset(DEgenes_temp_2528, abs(DEgenes_temp_2528$log2FoldChange) > log2(1.5))
  DEgenes_temp_2528 <- DEgenes_temp_2528[order(DEgenes_temp_2528$log2FoldChange, decreasing = TRUE), ]
  write.table(DEgenes_temp_2528, paste("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_VVD_vs_WT_temp",tempz ,"_DD",times,".txt", sep = '' ))
  
}}
  
  
#VVD vs WT for 28 degrees, multiple timepoints
for (times in c(0,4,8,12,16)){
  sample_temp_2528 <- subset(sampleTable, sampleTable$temp == '28' & sampleTable$time==times)
  DEdataset_temp_2528 <- DESeqDataSetFromHTSeqCount(sample_temp_2528, directory, ~type)
  DESeq_temp_2528 <- DESeq(DEdataset_temp_2528)
  res_temp_2528 <- results(DESeq_temp_2528)
  DEgenes_temp_2528 <- subset(res_temp_2528, res_temp_2528$padj <=0.01)
  DEgenes_temp_2528 <- subset(DEgenes_temp_2528, abs(DEgenes_temp_2528$log2FoldChange) > log2(1.5))
  DEgenes_temp_2528 <- DEgenes_temp_2528[order(DEgenes_temp_2528$log2FoldChange, decreasing = TRUE), ]
  write.table(DEgenes_temp_2528, paste("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_VVD_vs_WT_temp28_DD",times ))
  
}
#WT DD4 vs WT DD12 at both temperatures... also loops for VVD. LOOPCEPTION
for (strain in c("wt", "vvd")){
for (timevalues in c(25, 28)){
  sample_temp_2528 <- subset(sampleTable, sampleTable$type == strain & sampleTable$temp == timevalues & !sampleTable$time== 0 & !sampleTable$time== 8 & !sampleTable$time== 12)
  DEdataset_temp_2528 <- DESeqDataSetFromHTSeqCount(sample_temp_2528, directory, ~time)
  DESeq_temp_2528 <- DESeq(DEdataset_temp_2528)
  res_temp_2528 <- results(DESeq_temp_2528)
  DEgenes_temp_2528 <- subset(res_temp_2528, res_temp_2528$padj <=0.01)
  DEgenes_temp_2528 <- subset(DEgenes_temp_2528, abs(DEgenes_temp_2528$log2FoldChange) > log2(1.5))
  DEgenes_temp_2528 <- DEgenes_temp_2528[order(DEgenes_temp_2528$log2FoldChange, decreasing = TRUE), ]
  write.table(DEgenes_temp_2528, paste("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_", strain, "_DD4vs16_Temp", timevalues, sep = ""))
}}

  
  
  vennmaker1 <- read.csv("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_DD4vs12_Temp25", sep = " ")
  vennmaker2 <- read.csv("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_DD4vs12_Temp28", sep = " ") 
  vennmaker3 <- read.csv("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_DD4vs12_Temp25", sep = " ")
  vennmaker4 <- read.csv("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_DD4vs12_Temp28", sep = " ")
  vennmaker5 <- read.csv("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_DD4vs16_Temp25", sep = " ")
  vennmaker6 <- read.csv("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_DD4vs16_Temp28", sep = " ") 
  vennmaker7 <- read.csv("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_DD4vs16_Temp25", sep = " ")
  vennmaker8 <- read.csv("/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_vvd_DD4vs16_Temp28", sep = " ")
  #wtDD412_25 == WT, DD 4 to 12, temp 25
  wtDD412_25 <- rownames(vennmaker1)
  wtDD412_28 <- rownames(vennmaker2) 
  vvdDD412_25 <- rownames(vennmaker3)
  vvdDD412_28 <- rownames(vennmaker4) 
  wtDD416_25 <- rownames(vennmaker5)
  wtDD416_28 <- rownames(vennmaker6) 
  vvdDD416_25 <- rownames(vennmaker7)
  vvdDD416_28 <- rownames(vennmaker8)
  
  
venn.diagram(
  x = list(
    wt25_4_12 = wtDD412_25,
    wt28_4_12 = wtDD412_28,
    vvd25_4_16 = vvdDD416_25,
    vvd28_4_16 = vvdDD416_28
  ),
  filename = "/Users/langyy/Desktop/Differential_expression/DEgenes_tables/VennDiagram2.png",
  col = "black",
  lty = "dotted",
  lwd = 3,
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1,
  cat.fontfamily = "serif"
);



  # back to  older code.
### Condition: dark, 25; vvd vs wt
sample_25_dark <- subset(sampleTable,sampleTable$temp==25 & sampleTable$luminance=="dark")
DEdataset_25_dark <- DESeqDataSetFromHTSeqCount(sample_25_dark,directory,~type)
DESeq_25_dark <- DESeq(DEdataset_25_dark)
res_25_dark <- results(DESeq_25_dark)
DEgenes_25_dark <- subset(res_25_dark,res_25_dark$padj<=0.01)
DEgenes_25_dark <- subset(DEgenes_25_dark,abs(DEgenes_25_dark$log2FoldChange)>log2(1.5))
DEgenes_25_dark <- DEgenes_25_dark[order(DEgenes_25_dark$log2FoldChange,decreasing=TRUE),]
write.table(DEgenes_25_dark,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_25_dark.txt")
plotCounts(DESeq_25_dark,"NCU09968",intgroup=c("type","time"),main="NCU09968 counts at 25??C,dark")

### Condition: light, 28; vvd vs wt
sample_28_light <- subset(sampleTable,sampleTable$temp==28 & sampleTable$luminance=="light")
DEdataset_28_light <- DESeqDataSetFromHTSeqCount(sample_28_light,directory, ~type)
DESeq_28_light <- DESeq(DEdataset_28_light)
res_28_light <- results(DESeq_28_light)
DEgenes_28_light <- subset(res_28_light,res_28_light$padj<0.01)
DEgenes_28_light <- subset(DEgenes_28_light,abs(DEgenes_28_light$log2FoldChange)>log2(1.5))
DEgenes_28_light <- DEgenes_28_light[order(DEgenes_28_light$log2FoldChange,decreasing=TRUE),]
write.table(DEgenes_28_light,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_28_light.txt")
#plotCounts(DESeq_28_light,"NCU09968",intgroup=c("type","time"),main="NCU09968 counts at 28??C,light")


### Condition: dark, 28; vvd vs wt
sample_28_dark <- subset(sampleTable,sampleTable$temp==28 & sampleTable$luminance=="dark")
DEdataset_28_dark <- DESeqDataSetFromHTSeqCount(sample_28_dark,directory,~type)
DESeq_28_dark <- DESeq(DEdataset_28_dark)
res_28_dark <- results(DESeq_28_dark)
DEgenes_28_dark <- subset(res_28_dark,res_28_dark$padj<=0.01)
DEgenes_28_dark <- subset(DEgenes_28_dark,abs(DEgenes_28_dark$log2FoldChange)>log2(1.5))
DEgenes_28_dark <- DEgenes_28_dark[order(DEgenes_28_dark$log2FoldChange,decreasing=TRUE),]
write.table(DEgenes_28_dark,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_28_dark.txt")
plotCounts(DESeq_28_dark,"NCU09968",intgroup=c("type","time"),main="NCU09968 counts at 28??C,dark")
write.table(DEgenes_28_dark,"/Users/langyy/Desktop/Differential_expression/DEgenes_28_dark.txt")

### wt vs vvd, 28, timecourse
sample_28 <- subset(sampleTable,sampleTable$temp==28)
DEdataset_28 <- DESeqDataSetFromHTSeqCount(sample_28,directory,~type)
DESeq_28 <- DESeq(DEdataset_28)
res_28 <- results(DESeq_28)
plotCounts(DESeq_28,"NCU03639",intgroup = c("type","time"),main="NCU03639 counts at 28??C")
plotCounts(DESeq_28,"NCU07825",intgroup = c("type","time"),main="NCU07825 counts at 28??C")
plotCounts(DESeq_28,"NCU08770",intgroup = c("type","time"),main="NCU08770 counts at 28??C")

### wt vs vvd, 25, timecourse
sample_25 <- subset(sampleTable,sampleTable$temp==25)
DEdataset_25 <- DESeqDataSetFromHTSeqCount(sample_25,directory,~type)
DESeq_25 <- DESeq(DEdataset_25)
res_25 <- results(DESeq_25)
plotCounts(DESeq_25,"NCU08770",intgroup = c("type","time"),main="NCU08770 counts at 25??C")

### wt vs vvd, dark
sample_dark <- subset(sampleTable,sampleTable$luminance=="dark")
DEdataset_dark <- DESeqDataSetFromHTSeqCount(sample_dark,directory, ~type)
DESeq_dark <- DESeq(DEdataset_dark)
res_dark <- results(DESeq_dark)
plotCounts(DESeq_dark,"NCU09969",intgroup = c("type","time"),main="NCU09969 counts")

### wt 25 DD0 vs wt 28 DD0
sample_wt_0 <- subset(sampleTable,sampleTable$type=="wt" & sampleTable$time==0)
DEdataset_wt_0 <- DESeqDataSetFromHTSeqCount(sample_wt_0,directory,~temp)
DESeq_wt_0 <- DESeq(DEdataset_wt_0)
res_wt_0 <- results(DESeq_wt_0)
summary(res_wt_0)
DEgenes_wt_0 <- subset(res_wt_0,res_wt_0$padj<0.01)
DEgenes_wt_0 <- subset(DEgenes_wt_0,abs(DEgenes_wt_0$log2FoldChange)>log2(1.5))
DEgenes_wt_0 <- DEgenes_wt_0[order(DEgenes_wt_0$log2FoldChange,decreasing=TRUE),]
write.table(DEgenes_wt_0,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_0.txt")
wt_0_up <- subset(DEgenes_wt_0,DEgenes_wt_0$log2FoldChange>0)
wt_0_down <- subset(DEgenes_wt_0,DEgenes_wt_0$log2FoldChange<0)

### wt 25 DD4 vs wt 28 DD4
sample_wt_4 <- subset(sampleTable,sampleTable$type=="wt" & sampleTable$time==4)
DEdataset_wt_4 <- DESeqDataSetFromHTSeqCount(sample_wt_4,directory,~temp)
DESeq_wt_4 <- DESeq(DEdataset_wt_4)
res_wt_4 <- results(DESeq_wt_4)
summary(res_wt_4)
DEgenes_wt_4 <- subset(res_wt_4,res_wt_4$padj<0.01)
DEgenes_wt_4 <- subset(DEgenes_wt_4,abs(DEgenes_wt_4$log2FoldChange)>log2(1.5))
DEgenes_wt_4 <- DEgenes_wt_4[order(DEgenes_wt_4$log2FoldChange,decreasing=TRUE),]
write.table(DEgenes_wt_4,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_4.txt")
wt_4_up <- subset(DEgenes_wt_4,DEgenes_wt_4$log2FoldChange>0)
wt_4_down <- subset(DEgenes_wt_4,DEgenes_wt_4$log2FoldChange<0)

### wt 25 DD8 vs wt 28 DD8
sample_wt_8 <- subset(sampleTable,sampleTable$type=="wt" & sampleTable$time==8)
DEdataset_wt_8 <- DESeqDataSetFromHTSeqCount(sample_wt_8,directory, ~ temp)
DESeq_wt_8 <- DESeq(DEdataset_wt_8)
res_wt_8 <- results(DESeq_wt_8)
summary(res_wt_8)
DEgenes_wt_8 <- subset(res_wt_8,res_wt_8$padj<0.01)
DEgenes_wt_8 <- subset(DEgenes_wt_8,abs(DEgenes_wt_8$log2FoldChange)>log2(1.5))
DEgenes_wt_8 <- DEgenes_wt_8[order(DEgenes_wt_8$log2FoldChange,decreasing=TRUE),]
write.table(DEgenes_wt_8,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_8.txt")
wt_8_up <- subset(DEgenes_wt_8,DEgenes_wt_8$log2FoldChange>0)
wt_8_down <- subset(DEgenes_wt_8,DEgenes_wt_8$log2FoldChange<0)

### wt 25 DD12 vs wt 28 DD12
sample_wt_12 <- subset(sampleTable,sampleTable$type=="wt" & sampleTable$time==12)
DEdataset_wt_12 <- DESeqDataSetFromHTSeqCount(sample_wt_12,directory,~temp)
DESeq_wt_12 <- DESeq(DEdataset_wt_12)
res_wt_12 <- results(DESeq_wt_12)
summary(res_wt_12)
DEgenes_wt_12 <- subset(res_wt_12,res_wt_12$padj<0.01)
DEgenes_wt_12 <- subset(DEgenes_wt_12,abs(DEgenes_wt_12$log2FoldChange)>log2(1.5))
DEgenes_wt_12 <- DEgenes_wt_12[order(DEgenes_wt_12$log2FoldChange,decreasing=TRUE),]
write.table(DEgenes_wt_12,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_12.txt")
wt_12_up <- subset(DEgenes_wt_12,DEgenes_wt_12$log2FoldChange>0)
wt_12_down <- subset(DEgenes_wt_12,DEgenes_wt_12$log2FoldChange<0)

### wt 25 DD16 vs wt 28 DD16
sample_wt_16 <- subset(sampleTable,sampleTable$type=="wt" & sampleTable$time==16)
DEdataset_wt_16 <- DESeqDataSetFromHTSeqCount(sample_wt_16,directory,~temp)
DESeq_wt_16 <- DESeq(DEdataset_wt_16)
res_wt_16 <- results(DESeq_wt_16)
summary(res_wt_16)
DEgenes_wt_16 <- subset(res_wt_16,res_wt_16$padj<0.01)
DEgenes_wt_16 <- subset(DEgenes_wt_16,abs(DEgenes_wt_16$log2FoldChange)>log2(1.5))
DEgenes_wt_16 <- DEgenes_wt_16[order(DEgenes_wt_16$log2FoldChange,decreasing=TRUE),]
write.table(DEgenes_wt_16,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_wt_16.txt")
wt_16_up <- subset(DEgenes_wt_16,DEgenes_wt_16$log2FoldChange>0)
wt_16_down <- subset(DEgenes_wt_16,DEgenes_wt_16$log2FoldChange<0)

### up and down regulated genes count
count_up <- c(nrow(wt_0_up),nrow(wt_4_up),nrow(wt_8_up),nrow(wt_12_up),nrow(wt_16_up))
count_down <- c(nrow(wt_0_down),nrow(wt_4_down),nrow(wt_8_down),nrow(wt_12_down),nrow(wt_16_down))
count_up_down <- cbind(count_up,count_down)
row.names(count_up_down) <- c("0h","4h","8h","12h","16h")
count_up_down <- melt(count_up_down)
colnames(count_up_down) <- c("time","regulation","count")
ggplot(count_up_down,aes(x=time,y=count,fill=regulation)) + 
  geom_bar(stat="identity",position=position_dodge()) +
  geom_text(aes(label=count),position=position_dodge(0.9),vjust=-0.5)+
  scale_fill_discrete(labels=c("down","up"),guide = guide_legend(reverse=TRUE))

### wt 25 vs wt 28
sample_wt <- subset(sampleTable,sampleTable$type=="wt")
DEdataset_wt <- DESeqDataSetFromHTSeqCount(sample_wt,directory,~temp)
DESeq_wt <- DESeq(DEdataset_wt)
res_wt <- results(DESeq_wt)
DEgenes_wt <- subset(res_wt,res_wt$padj<0.01)
DEgenes_wt <- subset(DEgenes_wt,abs(DEgenes_wt$log2FoldChange)>log2(1.5))
DEgenes_wt <- DEgenes_wt[order(DEgenes_wt$log2FoldChange,decreasing = TRUE),]

### wt 25 vs wt 28 time course
sample_wt_time <- subset(sampleTable,sampleTable$type=="wt")
DEdataset_wt_time <- DESeqDataSetFromHTSeqCount(sample_wt_time,directory,~temp+time)
DESeq_wt_time <- DESeq(DEdataset_wt_time)
res_wt_time <- results(DESeq_wt_time)

### pass data to ggplot2
d <- plotCounts(DEdataset_wt_time,"NCU00701",intgroup = c("time","temp","rep"),returnData = TRUE)
d <- plotCounts(DEdataset_wt_time,"NCU02287",intgroup = c("time","temp","rep"),returnData = TRUE)
d <- plotCounts(DEdataset_wt_time,"NCU02370",intgroup = c("time","temp","rep"),returnData = TRUE)
d <- plotCounts(DEdataset_wt_time,"NCU03415",intgroup = c("time","temp","rep"),returnData = TRUE)
d <- plotCounts(DEdataset_wt_time,"NCU03771",intgroup = c("time","temp","rep"),returnData = TRUE)
d <- plotCounts(DEdataset_wt_time,"NCU04298",intgroup = c("time","temp","rep"),returnData = TRUE)
d <- plotCounts(DEdataset_wt_time,"NCU04603",intgroup = c("time","temp","rep"),returnData = TRUE)
d <- plotCounts(DEdataset_wt_time,"NCU05185",intgroup = c("time","temp","rep"),returnData = TRUE)
d <- plotCounts(DEdataset_wt_time,"NCU05209",intgroup = c("time","temp","rep"),returnData = TRUE)
ggplot(d, aes(x=time, y=count,color=temp,group=interaction(temp,rep),shape=rep))+
  geom_line()+geom_point()+
  geom_point(position=position_jitter(w=0.1,h=0)) 
  

  

### generate gene list
genelist <- c(rownames(DEgenes_wt_0),rownames(DEgenes_wt_4),rownames(DEgenes_wt_8),rownames(DEgenes_wt_12),rownames(DEgenes_wt_16))  
genelist <- genelist[order(genelist)]
write(genelist,"/Users/langyiying/Dropbox_Work/Dropbox/Msc_Lectures/RP2/genelist.txt")

### Venn diagram between: 1)Condition: light, 25; vvd vs wt 2)Condition: dark, 25; vvd vs wt
#                         3)Condition: light, 25; vvd vs wt 4)Condition: dark, 25; vvd vs wt
venn.plot <- venn.diagram(x=list("light 25"=rownames(DEgenes_25_light),
                                 "dark 25"=rownames(DEgenes_25_dark),
                                 "light 28"=rownames(DEgenes_28_light),
                                 "dark 28"=rownames(DEgenes_28_dark)),
                                 filename = "/Users/langyy/Desktop/Differential_expression/plot/venn_vvd_vs_wt.tiff",
                                 col="transparent",fill=c("cornflowerblue","green","yellow","darkorchid1"),
                                 cex=1.5)

#plotCounts(DESeq_25_light,"NCU09968",intgroup=c("type","luminance"))

# intersection of the comparison of vvd vs wt: 25_light, 25_dark, 28_light, 28_dark 
group_e_h <- Reduce(intersect,list(rownames(DEgenes_25_light),rownames(DEgenes_28_light),rownames(DEgenes_28_dark)))
group_e_h_f <- Reduce(intersect,list(rownames(DEgenes_28_light),rownames(DEgenes_28_dark)))
group_h_i <- Reduce(intersect,list(rownames(DEgenes_25_dark),rownames(DEgenes_28_dark)))
group_e_h_d <- Reduce(intersect,list(rownames(DEgenes_25_light),rownames(DEgenes_28_dark)))
group_b_e_h <- Reduce(intersect,list(rownames(DEgenes_25_light),rownames(DEgenes_28_light)))
group_b_d_e_h <- union(group_e_h_d,group_b)
group_b_e_f_h <- union(group_b, group_e_h_f)
group_d_e_f_h <- union(group_e_h_f,group_d)
group_d_e_f_h_i <- union(group_d_e_f_h,group_i)

group_a <- setdiff(rownames(DEgenes_25_light),group_b_d_e_h)
group_b <- setdiff(group_b_e_h,group_e_h)
group_c <- setdiff(rownames(DEgenes_28_light),group_b_e_f_h)
group_d <- setdiff(group_e_h_d,group_e_h)
group_e <- setdiff(group_e_h,group_h)
group_f <- setdiff(group_e_h_f,group_e_h)
group_g <- setdiff(rownames(DEgenes_28_dark),group_d_e_f_h_i)
group_h <- Reduce(intersect,list(rownames(DEgenes_25_light),rownames(DEgenes_25_dark),rownames(DEgenes_28_light),rownames(DEgenes_28_dark)))
group_i <- setdiff(group_h_i,group_h)

### NCU03639
plotCounts(DESeq_28_light,"NCU03639",intgroup=c("type","time"),main="NCU03639 counts at 28??C,light")
plotCounts(DESeq_28_dark,"NCU03639",intgroup=c("type","time"),main="NCU03639 counts at 28??C,dark")







### compare time course of temperature 25 and 28 in vvd knocked-out type
temp_in_vvd <- DESeqDataSetFromHTSeqCount(sampleTable[sampleTable$type=="vvd",],directory,~ temp + time + temp:time)
temp_in_vvd_TC <- DESeq(temp_in_vvd) 
res_temp_in_vvd <- results(temp_in_vvd_TC)
ordered_res_temp_vvd <- res_temp_in_vvd[order(res_temp_in_vvd$padj),]
a<-subset(ordered_res_temp_vvd,ordered_res_temp_vvd$padj<0.05)
data <- plotCounts(temp_in_vvd_TC,"NCU08721",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU08721 in vvd") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()
data <- plotCounts(temp_in_vvd_TC,"NCU09506",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU09506 in vvd") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()
data <- plotCounts(temp_in_vvd_TC,"NCU00987",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU00987 in vvd") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()
data <- plotCounts(temp_in_vvd_TC,"NCU09505",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU09505 in vvd") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()
data <- plotCounts(temp_in_vvd_TC,"NCU08179",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU08179 in vvd") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()
data <- plotCounts(temp_in_vvd_TC,"NCU09773",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU09773 in vvd") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()


### compare time course of temperature 25 and 28 in wild-type type
temp_in_wt <- DESeqDataSetFromHTSeqCount(sampleTable[sampleTable$type=="wt",],directory,~ temp + time + temp:time)
temp_in_wt_TC <- DESeq(temp_in_wt)
res_temp_in_wt <- results(temp_in_wt_TC)
ordered_res_temp_wt <- res_temp_in_wt[order(res_temp_in_wt$padj),]
b<-subset(ordered_res_temp_wt,ordered_res_temp_wt$padj<0.01)

head(res_temp_in_wt[order(res_temp_in_wt$padj),],200)
data <- plotCounts(temp_in_wt_TC,"NCU08721",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU08721 in wt") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()
data <- plotCounts(temp_in_wt_TC,"NCU09506",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU09506 in wt") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()
data <- plotCounts(temp_in_wt_TC,"NCU00987",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU00987 in wt") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()
data <- plotCounts(temp_in_wt_TC,"NCU09505",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU09505 in wt") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()
data <- plotCounts(temp_in_wt_TC,"NCU08179",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU08179 in wt") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()
data <- plotCounts(temp_in_wt_TC,"NCU09773",intgroup=c("temp","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=temp,group=temp)) + ggtitle("NCU09773 in wt") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()

### intersection
intersect(rownames(a),rownames(b))

### Test the effect of different types in different time point at 25
types_at_25 <- DESeqDataSetFromHTSeqCount(sampleTable[sampleTable$temp==25,],directory,~ type + time + type:time)
types_at_25_TC <- DESeq(types_at_25)
res_types_at_25 <- results(types_at_25_TC)
ordered_res_types_25 <- res_types_at_25[order(res_types_at_25$padj),]
c<-subset(ordered_res_types_25,ordered_res_types_25$padj<0.01)

data <- plotCounts(types_at_25_TC,"NCU02031",intgroup=c("type","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=type,group=type)) + ggtitle("NCU02031 at 25") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()




### Test the effect of different types in different time point at 28
types_at_28 <- DESeqDataSetFromHTSeqCount(sampleTable[sampleTable$temp==28,],directory,~ type + time + type:time)
types_at_28_TC <- DESeq(types_at_28)
res_types_at_28 <- results(types_at_28_TC)
ordered_res_types_28 <- res_types_at_28[order(res_types_at_28$padj),]
d<-subset(ordered_res_types_28,ordered_res_types_28$padj<0.01)

data <- plotCounts(types_at_28_TC,"NCU02031",intgroup=c("type","time"),returnData = TRUE)
ggplot(data, aes(x=time,y=count,color=type,group=type)) + ggtitle("NCU02031 at 28") + geom_point()+ stat_smooth(se=FALSE,method="loess")+  scale_y_log10()


### intersection
intersect(rownames(c),rownames(d))


















design(temp_in_vvd)<-formula(~ temp + time + temp:time)
test <- DESeq(test)
restest <- results(test)
data3 <- plotCounts(test, which.min(restest$padj), 
                    intgroup=c("temp","time"), returnData=TRUE)
ggplot(data3, aes(x=time, y=count, color=temp, group=temp)) 
+ 
  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()


restest <- results(test,contrast = c("temp","25","28"))

dds <- DESeq(ddshtseq)
res <-results(dds)
res_temp <- results(dds,contrast = c("temp",25,28))

dds_temp <- DESeqDataSetFromHTSeqCount(sampleTable,directory,~ temp )
dds_temp1 <- DESeq(dds_temp)
res_temp1 <- results(dds_temp1)

summary(res)
sum(res$padj<0.1,na.rm=TRUE)
resordered <- res[order(res$padj),]
head(resordered,150)
plotCounts(dds, gene=which(row.names(res)=="NCU03967"), intgroup="type")
plotMA(res, main="DESeq2", ylim=c(-2,2))
pcaplot(dds,intgroup = c("type"))



### heatmap
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("type")])
pheatmap(assay(dds), cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)




