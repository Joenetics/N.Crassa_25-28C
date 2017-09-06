library("DESeq2")
library("ggplot2")
library("pcaExplorer")
library("VennDiagram")
library("reshape")

# pathway… MODIFY THIS TO ‘htseq_count’ folder/directory!
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

plotCounts(DESeq_all,"NCU09968",intgroup=c("type","temp","time")) #hypothetical protein
plotCounts(DESeq_all,"NCU03967",intgroup=c("type","temp","luminance"),normalized=TRUE) #vvd
d<-plotCounts(DESeq_all,"NCU03967",intgroup=c("type","temp"),normalized=TRUE,returnData = TRUE) #vvd
plotCounts(DESeq_all,"NCU16528",intgroup=c("type","temp","luminance"),normalized=TRUE)
plotCounts(DESeq_all,"NCU08626",intgroup=c("type","temp","luminance")) #hypothetical protein


### Condition: light, 25; vvd vs wt
sample_25_light <- subset(sampleTable,sampleTable$temp==25 & sampleTable$luminance=="light")
DEdataset_25_light <- DESeqDataSetFromHTSeqCount(sample_25_light,directory, ~type)
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
write.table(DEgenes_25_light,"/Users/langyy/Desktop/Differential_expression/DEgenes_tables/DEgenes_25_light.txt")

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




