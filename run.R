###Code for WGBS project######
options(scipen = 20)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(plyr)
library(DelayedMatrixStats)
library(tibble)
library(valr)
library(ggbeeswarm)
library(ggpubr)
library(VennDiagram)
library(dmrseq)
library(bsseq)
library(gridExtra)
library(readr)
library(vipor)
library(dplyr)
library(formattable)
library(stringr)
library(reshape2)
library(BSgenome)
library(MutationalPatterns)
################################################meta_info################################################
EAC_Tumor_sampleList=c(paste0("EAC_", c(1:4,6)), paste0("GEJ_", 1:7))
ESCC_Tumor_sampleList=paste0("ESCC_", c(1:17,19:22))
EAC_Nonmalignant_sampleList=paste0("GEJ_Nonmalignant_", 1:7)
ESCC_Nonmalignant_sampleList=paste0("ESCC_Nonmalignant_", c(1:3, "53F", "54M"))
types=c(rep("ESCC_Tumor",21),rep("ESCC_Nonmalignant", 5), rep("EAC/GEJ_Tumor",12), rep("GEJ_Nonmalignantl",7))

fileList=c(ESCC_Tumor_sampleList, ESCC_Nonmalignant_sampleList, EAC_Tumor_sampleList, EAC_Nonmalignant_sampleList)
chrListTarget=paste("chr",1:22,sep="")

annotation_row=data.frame(Type=c(rep("ESCC_Tumor", 21), rep("ESCC_Nonmalignant", 5), rep("EAC/GEJ_Tumor", 12), rep("GEJ_Nonmalignant", 7)), stringsAsFactors = F)
rownames(annotation_row)=fileList
annotation_row$Type=factor(annotation_row$Type, levels=c("ESCC_Nonmalignant", "GEJ_Nonmalignant", "ESCC_Tumor", "EAC/GEJ_Tumor"))
annotation_row=annotation_row[order(annotation_row$Type), ,drop=F]
annotation_col=annotation_row
annotation_col$Sample=rownames(annotation_col)

ann_colors = list(Type = c("#F7CE46", "#75FBFD", "#EA3323", "#0000F5"))
names(ann_colors[[1]])=c("ESCC_Nonmalignant", "GEJ_Nonmalignant", "ESCC_Tumor", "EAC/GEJ_Tumor")

my_theme=theme_classic()+theme(axis.text = element_text(color="black", size=10), 
                               axis.title.y=element_text(color="black", size=12),
                               plot.title=element_text(hjust = 0.5, face="bold", size=14))

################################################PMD calling################################################
###ESCA samples
##Call PMDs by MMSeekR
library(MMSeekR.data)
library(MMSeekR)

changeToTabFile=function(sampleIndex){
  data=read_tsv(paste0("Data/ESCA_rmblackList_cov5/", sampleIndex, ".all.cov5.sorted.rmblackList.bed"), col_names = F)
  data=data[,c(1:3,5,4)]
  colnames(data)=c("Chr","Start","End","T","beta")
  data$M=round(data$T*data$beta,0)
  chrListTarget=paste("chr",1:22,sep="")
  data=data[data$Chr%in%chrListTarget,]
  data$Pos=data$Start+1
  data=data[,c(1,7,4,6)]
  write.table(data, file=paste0("Data/MMSeekR_PMDs/TabFile/", sampleIndex, "methyCov.tab"), sep="\t", row.names = F, col.names = F, quote=F)
}
for(fileIndex in fileList){
  print(fileIndex)
  changeToTabFile(fileIndex)
}
methFileList=dir("Data/MMSeekR_PMDs/TabFile/", pattern=".tab", full.names = T)
for(methFile in methFileList){
  fileIndex=gsub(".tab", "", basename(methFile))
  runPMDs("hg38", methFile, paste0("Data/MMSeekR_PMDs/PMDs/", fileIndex))
}
# results in Data/MMSeekR_PMDs/PMDs/

##Call PMDs by MethylSeekR using default options
# results in Data/MethylSeekR_PMDs/ESCA/

##Call PMDs by MethPipe using default options
# results in Data/MethPipe_PMDs/ESCA/

###BLUEPRINT_Tumor
##Call PMDs by MMSeekR
# results in Data/MMSeekR_PMDs/BLUEPRINT_Tumor/

##Call PMDs by MethylSeekR using default options
# results in Data/MethylSeekR_PMDs/BLUEPRINT_Tumor/

##Call PMDs by MethPipe using default options
# results in Data/MMSeekR_PMDs/BLUEPRINT_Tumor/

###get ESCA_union_PMDs
cancerTypeCommonPMDs=function(PMDfile, cutoff){
  data=read.table(PMDfile,sep="\t", stringsAsFactors = F,header=T)
  data=data[data$num>=cutoff,]
  outputFile=gsub("_PMDs_multiinter.bed", "_commonPMDs.bed",PMDfile)
  write.table(data[,1:3], file=outputFile, sep="\t", row.names = F, col.names=F, quote=F)
  temp=read_bed(outputFile)
  temp=bed_merge(temp)
  temp=as.data.frame(temp)
  ###merge the region and require the region longer than 2kb
  temp=temp[temp$end-temp$start>2000,]
  print(sum(temp$end-temp$start))
  write.table(temp, file=outputFile, sep="\t", row.names = F, col.names=F, quote=F)
}
cancerTypeCommonPMDs("Data/MMSeekR_PMDs/EAC_PMDs_multiinter.bed", 8)
cancerTypeCommonPMDs("Data/MMSeekR_PMDs/ESCC_PMDs_multiinter.bed", 14)

###get EAC and ESCC specific PMDs, shared PMD and HMDs
specificPMDs=function(sharedPMDFile, multiinterFile, cutoff, outputFile){
  sharedPMDs=read_bed(sharedPMDFile)
  sharedPMDs=bed_sort(sharedPMDs)
  sharedPMDs=bed_merge(sharedPMDs)
  sharedPMDs_length=sharedPMDs$end-sharedPMDs$start
  sharedPMDs_length_total=sum(sharedPMDs_length)
  print(sharedPMDs_length_total)
  multiinterPMDsAll=read.table(multiinterFile,sep="\t", stringsAsFactors = F,header=T)
  sampleLength=length(colnames(multiinterPMDsAll)[-1:-5])
  multiinterPMDs=multiinterPMDsAll[,1:4]
  multiinterPMDs=multiinterPMDs[multiinterPMDs$num>=(cutoff*sampleLength),]
  multiinterPMDs=tibble(chrom=multiinterPMDs$chrom, start=multiinterPMDs$start, end=multiinterPMDs$end)
  multiinterPMDs=bed_sort(multiinterPMDs)
  multiinterPMDs=bed_merge(multiinterPMDs)
  specificPMDs=bed_subtract(sharedPMDs, multiinterPMDs)
  specificPMDs=bed_sort(specificPMDs)
  specificPMDs=bed_merge(specificPMDs)
  specificPMDs$length=specificPMDs$end-specificPMDs$start
  specificPMDs=specificPMDs[specificPMDs$length>2000,]
  print(sum(specificPMDs$length))
  write.table(specificPMDs[,1:3], file=outputFile, row.names = F, col.names = F, sep="\t", quote=F)
}
#EAC only PMDs
specificPMDs("Data/MMSeekR_PMDs/EAC_commonPMDs.bed", "Data/MMSeekR_PMDs/ESCC_PMDs_multiinter.bed", (1/3), "Data/MMSeekR_PMDs/EAC_specificPMDs.bed")
#ESCC only PMDs
specificPMDs("Data/MMSeekR_PMDs/ESCC_commonPMDs.bed", "Data/MMSeekR_PMDs/EAC_PMDs_multiinter.bed", (1/3), "Data/MMSeekR_PMDs/ESCC_specificPMDs.bed")

#ESCA sharedHMDs
sharedHMDs=function(multiinterFile1, multiinterFile2, cutoff, outputFile){
  totalGenome=read_bed("meta/hg38_length.bed")
  blackList=read_bed("meta/hg38-blacklist.v2.change.bed",n_fields = 4)
  totalGenome=bed_subtract(totalGenome, blackList)
  
  multiinterPMDsAll1=read.table(multiinterFile1,sep="\t", stringsAsFactors = F,header=T)
  multiinterPMDsAll1_region=tibble(chrom=multiinterPMDsAll1$chrom, start=as.integer(multiinterPMDsAll1$start), 
                                   end=as.integer(multiinterPMDsAll1$end), num=as.integer(multiinterPMDsAll1$num))
  totalGenome=bed_subtract(totalGenome, multiinterPMDsAll1_region)
  sampleLength=length(colnames(multiinterPMDsAll1)[-1:-5])
  multiinterPMDsAll1_region=multiinterPMDsAll1_region[multiinterPMDsAll1_region$num<(cutoff*sampleLength),]
  
  multiinterPMDsAll2=read.table(multiinterFile2,sep="\t", stringsAsFactors = F,header=T)
  multiinterPMDsAll2_region=tibble(chrom=multiinterPMDsAll2$chrom, start=as.integer(multiinterPMDsAll2$start), 
                                   end=as.integer(multiinterPMDsAll2$end), num=as.integer(multiinterPMDsAll2$num))
  totalGenome=bed_subtract(totalGenome, multiinterPMDsAll2_region)
  sampleLength=length(colnames(multiinterPMDsAll2)[-1:-5])
  multiinterPMDsAll2_region=multiinterPMDsAll2_region[multiinterPMDsAll2_region$num<(cutoff*sampleLength),]
  if(nrow(multiinterPMDsAll2_region)>0&nrow(multiinterPMDsAll1_region)>0){
    multiinterPMDsAll12_region=bed_intersect(multiinterPMDsAll1_region, multiinterPMDsAll2_region)
    multiinterPMDsAll12_region$start=rowMaxs(as.matrix(multiinterPMDsAll12_region[,colnames(multiinterPMDsAll12_region)%in%c("start.x", "start.y")]))
    multiinterPMDsAll12_region$end=rowMins(as.matrix(multiinterPMDsAll12_region[,colnames(multiinterPMDsAll12_region)%in%c("end.x", "end.y")]))
    multiinterPMDsAll12_region=multiinterPMDsAll12_region[,colnames(multiinterPMDsAll12_region)%in%c("chrom", "start", "end")]
    totalGenome=rbind(totalGenome, multiinterPMDsAll12_region)
    totalGenome=bed_merge(totalGenome)
  }
  totalGenome$length=totalGenome$end-totalGenome$start
  totalGenome=totalGenome[totalGenome$length>2000,]
  print(sum(totalGenome$length))
  write.table(totalGenome[,1:3], file=outputFile, sep="\t", quote=F, col.names = F, row.names = F)
}
sharedHMDs("Data/MMSeekR_PMDs/ESCC_PMDs_multiinter.bed", "Data/MMSeekR_PMDs/EAC_PMDs_multiinter.bed", 1/3, "Data/MMSeekR_PMDs/ESCA_sharedHMDs.bed")
#ESCA sharedPMDs
sharedESCARegions=function(EACPMDFile, ESCCPMDFile, outputFile, type){
  EAC_PMDs=read_bed(EACPMDFile)
  EAC_PMDs_total=as.numeric(sum(EAC_PMDs$end-EAC_PMDs$start))
  ESCC_PMDs=read_bed(ESCCPMDFile)
  ESCC_PMDs_total=as.numeric(sum(ESCC_PMDs$end-ESCC_PMDs$start))
  ESCA_sharedPMDs=bed_intersect(EAC_PMDs, ESCC_PMDs)
  ESCA_sharedPMDs$start=rowMaxs(as.matrix(ESCA_sharedPMDs[,colnames(ESCA_sharedPMDs)%in%c("start.x","start.y")]))
  ESCA_sharedPMDs$end=rowMins(as.matrix(ESCA_sharedPMDs[,colnames(ESCA_sharedPMDs)%in%c("end.x","end.y")]))
  ESCA_sharedPMDs$length=ESCA_sharedPMDs$end-ESCA_sharedPMDs$start
  ESCA_sharedPMDs=ESCA_sharedPMDs[,colnames(ESCA_sharedPMDs)%in%c("chrom", "start", "end")]
  ESCA_sharedPMDs=bed_merge(bed_sort(ESCA_sharedPMDs))
  ESCA_sharedPMDs$length=ESCA_sharedPMDs$end-ESCA_sharedPMDs$start
  ESCA_sharedPMDs=ESCA_sharedPMDs[ESCA_sharedPMDs$length>2000,]
  ESCA_sharedPMDs_total=sum(ESCA_sharedPMDs$length)
  print(EAC_PMDs_total)
  print(ESCC_PMDs_total)
  print(ESCA_sharedPMDs_total)
  write.table(ESCA_sharedPMDs[,1:3], file=outputFile, sep="\t", row.names = F, col.names=F, quote=F)
}
sharedESCARegions("Data/MMSeekR_PMDs/EAC_commonPMDs.bed", "Data/MMSeekR_PMDs/ESCC_commonPMDs.bed","Data/MMSeekR_PMDs/ESCA_sharedPMDs.bed", "Tumor")


################################################Figure1B################################################
#bash Shell/Figure1B/getGlobalMeanBetaValues.sh global
#bash Shell/Figure1B/getRegionMeanBetaValues.sh meta/commonPMDs_hg38.bed
#bash Shell/Figure1B/getRegionMeanBetaValues.sh meta/gencode.v31.basic.promoter.Takai_Jones.CGI.bed
#bash Shell/Figure1B/getRegionMeanBetaValues.sh meta/hg38_repeat.change.LINE.bed
#bash Shell/Figure1B/getRegionMeanBetaValues.sh meta/hg38_repeat.change.LTR.bed
#bash Shell/Figure1B/getRegionMeanBetaValues.sh meta/hg38_repeat.change.SINE.bed
###combine the result and get plotdata "Figure1B/Methylaton_sixTypes.txt"
plotdata=read.table("Figure1B/Methylaton_sixTypes.txt", sep="\t", stringsAsFactors = F, header=T)
plotdata=plotdata[plotdata$Type%in%c("ESCC_Nonmalignant", "GEJ_Nonmalignant","ESCC_Tumor","EAC/GEJ_Tumor"),]
plotdata$Sample=factor(plotdata$Sample, levels=fileList)
plotdata$Type=factor(plotdata$Type, levels=c("ESCC_Nonmalignant", "GEJ_Nonmalignant","ESCC_Tumor","EAC/GEJ_Tumor"))
plotdata=melt(plotdata, id.vars = c("Sample", "Type"))
colnames(plotdata)=c("Sample", "Type", "DataType", "Methylation")
plotdata$DataType=factor(plotdata$DataType, levels=unique(plotdata$DataType))
plotDotPlot=function(myPlotData, saveFile){
  p<-ggplot(myPlotData, aes(x=Type, y=Methylation, fill=Type, color=Type)) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1.2) + facet_wrap(~ DataType, nrow=1)
  p=p+scale_color_manual(values=c("#EA3323", "#0000F5", "#EA3323", "#0000F5"))
  p=p+scale_fill_manual(values=ann_colors$Type[1:4])
  p=p+ylim(0,1)+xlab("")+ylab("Mean methylation")
  # p=p+ stat_summary(fun=mean, geom="point", shape=18, size=3, color="black")
  p=p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                     geom="pointrange", color="black")
  p=p+theme_bw()+theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(),
                       axis.text = element_text(color="black", size=10), 
                       axis.title.y=element_text(color="black", size=12))
  # print(p)
  pdf(saveFile, width=20, height=3.5)
  print(p)
  dev.off()
}
plotDotPlot(plotdata, "Figure1B/methylation_sixType.pdf")


################################################FigureS2A-D################################################
chr.sel="chr16"
library(MMSeekR.data)
library(MMSeekR)
data("NNscore.hg19")
data("hg19.seqLengths")
data("hg19.blackList")

####FigureS2BC
methFile <-"Data/FigureS2A-D//BP_venous_blood_S01FJZA1_NA_MantleCellLymphoma.methyCov.tab"
meth <- readMethylomeNew(fileName=methFile, NNdat=NNscore.hg19, seqLengths=hg19.seqLengths)
indx <- as.character(seqnames(meth))==chr.sel
meth=meth[indx]
hmm.modelList=trainPMDHMDNew(meth, "chr16", 201, 1, "FigureS2BC/FigureS2BC.pdf")
y.list=PMDviterbiSegNew(meth, hmm.modelList, 201, 1)
seg = createGRangesObjectPMDSegNew(meth, y.list, 1, hg19.seqLengths)
seg = tibble(chrom = as.character(seqnames(seg)), start = as.integer(start(seg)), 
             end = as.integer(end(seg)), type = seg$type, nCG=seg$nCG)
segRmBlackList = bed_subtract(seg, hg19.blackList)
seg=as.data.frame(seg)
targetSeg=seg[c(100,103),]
targetSeg=tibble(chrom=as.character(targetSeg$seqnames), start=targetSeg$start, end=targetSeg$end, type=targetSeg$type)

####FigureS2D
library(zoocat)
chr.sel="chr16"
nCGbin=201
num.cores=1
methTemp=as.data.frame(meth)
T=as.numeric(methTemp$T)
M=as.numeric(methTemp$M)
alphaScore <- MethylSeekR:::calculateAlphaDistr(M, T, nCGbin, num.cores)
methTemp$alphaScore=alphaScore
methTemp2=methTemp[!is.na(methTemp$NNscore),]
NNScore <- as.numeric(methTemp2$NNscore)
methylation=as.numeric(methTemp2$Methylation)
methylationMean=as.vector(runmean(Rle(methylation), k = nCGbin, na.rm = TRUE, endrule = "constant"))
methTemp2$methylationMean=methylationMean
methTemp2$MValue=log2((methTemp2$methylationMean+0.01)/(1-methTemp2$methylationMean+0.01))
corResult=as.vector(rollcor(NNScore, methylation, width = nCGbin,show = F, use="na.or.complete"))
corResult=c(rep(corResult[1],(nCGbin-1)/2), corResult, rep(corResult[length(corResult)],(nCGbin-1)/2))
methTemp2$cor=corResult
save(methTemp2, file="Data/FigureS2A-D/BP_venous_blood_S01FJZA1_NA_MantleCellLymphoma.methySeekR2.RData")

load("Data/FigureS2A-D/BP_venous_blood_S01FJZA1_NA_MantleCellLymphoma.methySeekR2.RData")
write.table(methTemp2, "Data/FigureS2A-D/FigureS2D.txt", row.names = F, col.names = TRUE, sep = "\t", quote=F)
p1=ggplot(methTemp2, aes(x=alphaScore, y=cor)) +theme_classic()+xlab("alpha score")+geom_point(color="darkblue", alpha=0.1, size=0.5)+ylab("PCC")+ggtitle("All CpG sites")
p1=p1 + stat_density_2d(aes(fill = ..level..), geom = "polygon") +gradient_fill("YlOrRd")
p1=p1+xlim(0,1.5)+ylim(-1,0.3)+theme(plot.title = element_text(hjust=0.5))
png("FigureS2D/FigureS2D.png", res=300, width=1300, height = 1000)
print(p1)
dev.off()

####FigureS2A
load("Data/FigureS2A-D/BP_venous_blood_S01FJZA1_NA_MantleCellLymphoma.methySeekR2.RData")
methTemp3=tibble(chrom=as.character(methTemp2$seqnames), start=as.integer(methTemp2$start), end=as.integer(methTemp2$end), 
                 Methylation=methTemp2$Methylation, NNscore=methTemp2$NNscore, cor=methTemp2$cor)
my_theme=theme_classic()+theme(axis.text = element_text(size=12), axis.title = element_text(size=14), plot.title = element_text(size=14, face="bold", hjust = 0.5))

methTemp3PMDs=bed_intersect(methTemp3, targetSeg[targetSeg$type%in%"PMD",])
methTemp3PMDs=as.data.frame(methTemp3PMDs[,1:6])
colnames(methTemp3PMDs)=colnames(methTemp3)
seed=5500
plotdata=methTemp3PMDs[(seed-100):(seed+100),]
region=paste0("PMD region (", plotdata[1,]$chrom, ":",plotdata[1,]$start, "-", plotdata[nrow(plotdata),]$end, ")")
print(paste0(plotdata[1,]$chrom, ":",plotdata[1,]$start, "-", plotdata[nrow(plotdata),]$end))
write.table(plotdata, file="FigureS2A/FigureS2A_PMD.txt", row.names = F, col.names = T, quote=F, sep="\t")
p1=ggplot(plotdata, aes(x=Methylation, y=NNscore)) + geom_point(size=2)+geom_smooth(method=lm)
p1=p1+xlab("CpG methylation")+ylab("CpG NNscore")+ggtitle(region)
p1=p1+annotate(geom="text", x=0.75, y=0.75, label=paste0("PCC = ", round(cor(plotdata$Methylation, plotdata$NNscore), 3)), color="red", size=6)
p1=p1+my_theme

methTemp3HMDs=bed_intersect(methTemp3, targetSeg[targetSeg$type%in%"notPMD",])
methTemp3HMDs=as.data.frame(methTemp3HMDs[,1:6])
colnames(methTemp3HMDs)=colnames(methTemp3)
seed=1500
plotdata=methTemp3HMDs[(seed-100):(seed+100),]
region=paste0("HMD region (", plotdata[1,]$chrom, ":",plotdata[1,]$start, "-", plotdata[nrow(plotdata),]$end, ")")
print(paste0(plotdata[1,]$chrom, ":",plotdata[1,]$start, "-", plotdata[nrow(plotdata),]$end))
write.table(plotdata, file="FigureS2A/FigureS2A_HMD.txt", row.names = F, col.names = T, quote=F, sep="\t")
p2=ggplot(plotdata, aes(x=Methylation, y=NNscore)) + geom_point(size=2)+geom_smooth(method=lm)
p2=p2+xlab("CpG methylation")+ylab("CpG NNscore")+ggtitle(region)
p2=p2+annotate(geom="text", x=0.75, y=0.75, label=paste0("PCC = ", round(cor(plotdata$Methylation, plotdata$NNscore), 3)), color="red", size=6)
p2=p2+my_theme

pdf("FigureS2A/FigureS2A.pdf", width=9, height=4)
print(ggarrange(p1, p2, nrow=1))
dev.off()

################################################Table S2################################################
ESCA_MethPipeList=dir("Data/MethPipe_PMDs/ESCA/", pattern=".bed", full.names = T)
ESCA_MethylSeekRList=dir("Data/MethylSeekR_PMDs/ESCA/", pattern=".bed", full.names = T)
ESCA_MMSeekRList=dir("Data/MMSeekR_PMDs/PMDs/", pattern=".bed", full.names = T)
BPTumor_MethPipeList=dir("Data/MethPipe_PMDs/BLUEPRINT_Tumor/", pattern=".bed", full.names = T)
BPTumor_MethylSeekRList=dir("Data/MethylSeekR_PMDs/BLUEPRINT_Tumor/", pattern=".bed", full.names = T)
BPTumor_MMSeekRPathList=dir("Data/MMSeekR_PMDs/BLUEPRINT_Tumor/", pattern=".bed", full.names = T)

commonPMD=read_bed("meta/commonPMDs_hg38.bed")
commonHMD=read_bed("meta/commonHMDs_hg38.bed")

getRatio.v2=function(fileList, software){
  result=as.data.frame(matrix(numeric(0),ncol=4))
  for(file in fileList){
    comomPMD_size=sum(commonPMD$end-commonPMD$start)
    fileIndex=gsub("\\..*", "", basename(file))
    fileIndex=gsub("_PMDs", "",fileIndex)
    tmpData=read_bed(file)
    tmpData$length=tmpData$end-tmpData$start
    size=sum(tmpData$length)
    tmpPMDData=bed_intersect(tmpData, commonPMD)
    overlapPMDSize=sum(tmpPMDData$.overlap)
    tmpResult=data.frame(sample=fileIndex, precision=overlapPMDSize/size, recall=overlapPMDSize/comomPMD_size)
    tmpResult$F1=2*tmpResult$recall*tmpResult$precision/(tmpResult$recall+tmpResult$precision)
    result=rbind(result, tmpResult)
  }
  colnames(result)=c("sample", paste0(software, "_precision"), paste0(software, "_recall"), paste0(software, "_F1"))
  return(result)
}
ESCA_MethPipeResult=getRatio.v2(ESCA_MethPipeList, "MethPipe")
ESCA_MethylSeekResult=getRatio.v2(ESCA_MethylSeekRList, "MethylSeekR")
ESCA_MMSeekResult=getRatio.v2(ESCA_MMSeekRList, "MMSeekR")
ESCA_result=merge(ESCA_MethPipeResult, ESCA_MethylSeekResult, by.x="sample", by.y="sample")
ESCA_result=merge(ESCA_result, ESCA_MMSeekResult, by.x="sample", by.y="sample")
ESCA_result=merge(annotation_col, ESCA_result, by.x="Sample", by.y="sample")
write.table(ESCA_result, file="TableS2/TableS2_ESCA.precision_recll.txt", row.names = F, col.names = T, sep="\t", quote=F)

BPTumor_MethPipeResult=getRatio.v2(BPTumor_MethPipeList, "MethPipe")
BPTumor_MethylSeekResult=getRatio.v2(BPTumor_MethylSeekRList, "MethylSeekR")
BPTumor_MMSeekResult=getRatio.v2(BPTumor_MMSeekRPathList, "MMSeekR")
BPTumor_result=merge(BPTumor_MethPipeResult, BPTumor_MethylSeekResult, by.x="sample", by.y="sample")
BPTumor_result=merge(BPTumor_result, BPTumor_MMSeekResult, by.x="sample", by.y="sample")
write.table(BPTumor_result, file="TableS2/TableS2_BLUEPRINT_Tumor.precision_recll.txt", row.names = F, col.names = T, sep="\t", quote=F)

####Method comparison for ESCA sample
ESCA_result=read.table("TableS2/TableS2_ESCA.precision_recll.txt", header=T, sep="\t")
ESCA_result$Type=factor(ESCA_result$Type, levels=c("ESCC_Nonmalignant", "GEJ_Nonmalignant", "ESCC_Tumor", "EAC/GEJ_Tumor"))

plotBeeswarm=function(myPlotdata, name, saveFile){
  colnames(myPlotdata)=c("Sample", "Type", "MethPipe", "MethylSeekR", "MMSeekR")
  myPlotdata=melt(myPlotdata, id.vars = c("Sample", "Type"))
  colnames(myPlotdata)[3:4]=c("Method", "Ratio")
  myPlotdata$Method=factor(myPlotdata$Method, levels=c("MethylSeekR", "MethPipe", "MMSeekR"))
  statistic_data=myPlotdata%>%group_by(Type, Method)%>%summarise(sd= sd(Ratio), Ratio= mean(Ratio))
  statistic_data$Method=factor(statistic_data$Method, levels=c("MethylSeekR", "MethPipe", "MMSeekR"))
  
  p=ggplot(myPlotdata, aes(x = Method, y = Ratio, color = Method)) +
    geom_beeswarm(cex = 0.5, alpha=0.4, size=1.5)+facet_grid(. ~ Type, scales = "free")+theme_classic()+xlab("")+
    scale_color_manual(values=c("blue", "#C900B8", "red"))+
    geom_errorbar(aes(ymin = Ratio-sd, ymax = Ratio+sd), width=0.2, data = statistic_data, color="black")
  
  p=p+stat_summary(fun = mean, fun.min = mean, fun.max = mean,
                   geom = "crossbar", width = 0.3, color="black", linewidth=0.3)
  p=p+theme(axis.text.y= element_text(size=8, color="black"), 
            axis.line = element_line(linewidth=0.2),
            axis.ticks = element_line(linewidth = 0.2),
            axis.text.x= element_blank())+xlab("")+ylab(name)
  
  pdf(saveFile, width=4.5, height=1.8)
  print(p)
  dev.off()
  write.table(myPlotdata, gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
}
plotdata1=ESCA_result[,c(1:2,5,8,11)]
plotBeeswarm(plotdata1, "F1","FigureS2FJ/ESCA.F1.beesswarm.pdf")

plotdata1=ESCA_result[,c(1:2,4,7,10)]
plotBeeswarm(plotdata1, "Recall","FigureS2FJ/ESCA.Recall.beesswarm.pdf")

plotdata1=ESCA_result[,c(1:2,3,6,9)]
plotBeeswarm(plotdata1, "Precision","FigureS2FJ/ESCA.Precision.beesswarm.pdf")

#####Method comparison for BP blood tumor samples
bluePrintSamples=read.table("Data/Figure1D_FigureS2EFG/BLUEPRINT_blood_sampleInfo.txt", sep="\t", header=T, stringsAsFactors = F)
bluePrintSamples=bluePrintSamples[order(bluePrintSamples$Type),]
bluePrintSamples_tumor=bluePrintSamples[!bluePrintSamples$Type%in%c("Myeloid","Lymphoid","Others"),]
bluePrintSamples_tumor=bluePrintSamples_tumor[order(bluePrintSamples_tumor$Type),]
bluePrintSamples_tumor$Type=factor(bluePrintSamples_tumor$Type, levels=unique(bluePrintSamples_tumor$Type))
bluePrintSamples_tumor$Subtype=factor(bluePrintSamples_tumor$Subtype, levels=c("AcuteLymphocyticLeukemia", "TcellProlymphocyticLeukemia",
                                                                               "ChronicLymphocyticLeukemia", "MantleCellLymphoma", "MultipleMyeloma", "AcuteMyeloidLeukemia"))
colnames(bluePrintSamples_tumor)[ncol(bluePrintSamples_tumor)]="PMD_methylation"
BLUEPRINT_Tumor_result=read.table("TableS2/TableS2_BLUEPRINT_Tumor.precision_recll.txt", header=T, sep="\t")
BLUEPRINT_Tumor_result=merge(bluePrintSamples_tumor[,c("Sample", "Subtype")], BLUEPRINT_Tumor_result, by.x="Sample", by.y="sample")
colnames(BLUEPRINT_Tumor_result)[2]=c("Type")
BLUEPRINT_Tumor_result$Type=factor(BLUEPRINT_Tumor_result$Type, levels=c("AcuteLymphocyticLeukemia",
                                                                         "ChronicLymphocyticLeukemia", "MultipleMyeloma", "TcellProlymphocyticLeukemia",
                                                                         "MantleCellLymphoma", "AcuteMyeloidLeukemia"))

plotdata1=BLUEPRINT_Tumor_result[,c(1:2,5,8,11)]
plotBeeswarm(plotdata1, "F1","FigureS2FJ/BPTumor.F1.beesswarm.pdf")

plotdata1=BLUEPRINT_Tumor_result[,c(1:2,4,7,10)]
plotBeeswarm(plotdata1, "Recall","FigureS2FJ/BPTumor.Recall.beesswarm.pdf")

plotdata1=BLUEPRINT_Tumor_result[,c(1:2,3,6,9)]
plotBeeswarm(plotdata1, "Precision","FigureS2FJ/BPTumor.Precision.beesswarm.pdf")

################################################Figure1D and FigureS2E-G################################################
#######Figure1D
hg38_commonPMD=read_bed("meta/commonPMDs_hg38.bed")
hg38_commonPMD=bed_merge(bed_sort(hg38_commonPMD))
hg38_commonHMD=read_bed("meta/commonHMDs_hg38.bed")
hg38_commonHMD=bed_merge(bed_sort(hg38_commonHMD))

hg38_regions=read_bed("meta/hg38_length.bed")
blackListResgions=read_bed("meta/hg38-blacklist.v2.change.bed", n_fields = 4)
hg38_regions_rmBlackList=bed_subtract(hg38_regions, blackListResgions)
hg38_regions_rmBlackList=bed_merge(bed_sort(hg38_regions_rmBlackList))
hg38_regions_rmBlackList_win=bed_makewindows(hg38_regions_rmBlackList, win_size=30000)

ESCA_sampleInfo=read.table("meta/ESCA_sampleInfo.txt", header=T, sep="\t", stringsAsFactors = F)
ESCA_PMD_methylation=read.table("Figure1B/Methylaton_sixTypes.txt", header = T, sep="\t")[,c(1,4)]
ESCA_PMD_methylation$Sample=factor(ESCA_PMD_methylation$Sample, levels=ESCA_sampleInfo$Sample)
ESCA_PMD_methylation=ESCA_PMD_methylation[order(ESCA_PMD_methylation$Sample),]
ESCA_sampleInfo$PMD_methylation=ESCA_PMD_methylation$commonPMD

methylSeekRFileList=dir("Data/MethylSeekR_PMDs/ESCA/", pattern = ".MethylSeekR.rmBlackList.bed", full.names = T)
methPipeFileList=dir("Data/MethPipe_PMDs/ESCA/", pattern = ".MethPipe.rmBlackList.bed", full.names = T)
model2D3DFileList=dir("Data/MMSeekR_PMDs/PMDs/", pattern = ".Model2D3D.rmBlackList.bed", full.names = T)

Cluster=function(pmdFileList, num, pattern, fileIndex, outputPath){
  ratioResult=as.data.frame(hg38_regions_rmBlackList_win)
  colnames(ratioResult)=c("chr", "start", "end", "id")
  stateResult=ratioResult
  overlapResult=c()
  for(PMDFile in pmdFileList){
    PMDFileIndex=gsub(pattern, "", basename(PMDFile))
    print(PMDFileIndex)
    PMDRegion=read_bed(PMDFile)
    PMDRegion=bed_merge(bed_sort(PMDRegion))
    temp=bed_intersect(hg38_regions_rmBlackList_win, PMDRegion)
    temp=unique(temp[,c(1:3,7)])
    temp=as.data.frame(temp)
    colnames(temp)=c("chr", "start", "end", "overlap")
    temp$name=paste0(temp$chr, ":", temp$start)
    temp1=temp[!temp$name%in%temp[duplicated(temp$name),]$name,]
    temp2=temp[temp$name%in%temp[duplicated(temp$name),]$name,]
    temp2=do.call(rbind, lapply(unique(temp2$name),function(x){
      tmpData=temp2[temp2$name%in%x,][1,]
      tmpData$overlap=sum(temp2[temp2$name%in%x,]$overlap)
      return(tmpData)
    }))
    temp=rbind(temp1, temp2)
    temp=temp[,1:4]
    overlapResult=c(overlapResult,sum(temp$overlap))
    temp$ratio=temp$overlap/(temp$end-temp$start)
    temp$state=0
    temp[temp$ratio>0.5,]$state=1
    ratioData=temp[,c(1:3,5)]
    colnames(ratioData)[4]=PMDFileIndex
    stateData=temp[,c(1:3,6)]
    colnames(stateData)[4]=PMDFileIndex
    ratioResult=merge(ratioResult, ratioData, by.x=c("chr", "start", "end"), by.y=c("chr", "start", "end"), all.x=T)
    stateResult=merge(stateResult, stateData, by.x=c("chr", "start", "end"), by.y=c("chr", "start", "end"), all.x=T)
  }
  
  overlapResult=data.frame(overlap=overlapResult)
  overlapResult$Sample=colnames(stateResult)[-1:-4]
  overlapResult$overlap=overlapResult$overlap/sum(hg38_regions_rmBlackList$end-hg38_regions_rmBlackList$start)
  
  rownames(ratioResult)=paste0(ratioResult$chr, ":",ratioResult$start, "-",ratioResult$end)
  rownames(stateResult)=paste0(stateResult$chr, ":",stateResult$start, "-",stateResult$end)
  ratioResult=ratioResult[,-1:-4]
  stateResult=stateResult[,-1:-4]
  ratioResult[is.na(ratioResult)]=0
  stateResult[is.na(stateResult)]=0
  
  stateResult2=stateResult[,colnames(stateResult)%in%ESCA_sampleInfo$Sample]
  ratioResult2=ratioResult[,colnames(stateResult)%in%ESCA_sampleInfo$Sample]
  print(ncol(stateResult2))
  overlapResult2=overlapResult[overlapResult$Sample%in%ESCA_sampleInfo$Sample,]
  save(stateResult2, ratioResult2, overlapResult2, file=paste0(outputPath, "/", fileIndex, ".RData"))
}
Cluster(methylSeekRFileList, 5000, ".MethylSeekR.rmBlackList.bed", "ESCA_hg38_MethylSeekR_30kb", "Figure1D_FigureS2EFG/")
Cluster(methPipeFileList, 5000, ".MethPipe.rmBlackList.bed", "ESCA_hg38_MethPipe_30kb", "Figure1D_FigureS2EFG/")
Cluster(model2D3DFileList, 5000, ".Model2D3D.rmBlackList.bed", "ESCA_hg38_MMSeekR_30kb", "Figure1D_FigureS2EFG/")

plotPCA=function(dataSet, num, saveFile){
  varResult=rowVars(as.matrix(dataSet))
  clusterData=dataSet[order(varResult, decreasing = T), ][1:num,]
  
  targetSamples=ESCA_sampleInfo[ESCA_sampleInfo$Sample%in%colnames(dataSet),]
  targetSamples$Type=factor(targetSamples$Type, levels=c("ESCC_Tumor", "ESCC_Nonmalignant", "EAC/GEJ_Tumor", "GEJ_Nonmalignant"))
  targetSamples=targetSamples[order(targetSamples$Type),]
  targetSamples$Color=factor(targetSamples$Color, levels=unique(targetSamples$Color))
  
  targetSamples$Sample=factor(targetSamples$Sample, levels=colnames(dataSet))
  targetSamples=targetSamples[order(targetSamples$Sample),]
  colorsData=targetSamples
  rownames(colorsData)=colorsData$Sample
  colorsData=colorsData[,-1]
  
  plotdata=t(clusterData)
  pca <- prcomp(plotdata,scale = TRUE)
  PCNum=c(1,2,3)
  xlab <- paste("PC", PCNum[1], "(",round((summary(pca))$importance[2,PCNum[1]]*100,1),"%)",sep="")
  ylab1 <- paste("PC", PCNum[2], "(",round((summary(pca))$importance[2,PCNum[2]]*100,1),"%)",sep="")
  ylab2 <- paste("PC", PCNum[3], "(",round((summary(pca))$importance[2,PCNum[3]]*100,1),"%)",sep="")
  x<-paste0("PC", PCNum[1])
  y1<-paste0("PC", PCNum[2])
  y2<-paste0("PC", PCNum[3])
  data_x <- data.frame(varnames=rownames(plotdata), pca$x[,1:5])
  data_x <- cbind(data_x, colorsData)
  
  write.table(data_x, file=gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
  my_theme=theme_bw()+theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                            panel.grid.major=element_blank(),axis.title=element_text(color="black",size=12),
                            axis.text=element_text(size=12), legend.title = element_text(hjust=0.5, face="bold"))
  
  p1 <- ggplot(data_x, aes(PC1,PC2, color=Type, fill=Type))+geom_point(size=4, alpha=0.8, pch=21)+coord_equal(ratio=1)+xlab(xlab)+ylab(ylab1)
  p1 <- p1+scale_color_manual(values =c("#EA3323","#EA3323", "#0000F5", "#0000F5"))+scale_fill_manual(values =levels(targetSamples$Color))+my_theme
  pdf(saveFile,width = 5.5,height =4)
  print(ggarrange(p1))
  dev.off()
}
Cluster2=function(outputPath,num, fileIndex, saveFile){
  load(paste0(outputPath, "/", fileIndex, ".RData"))
  plotPCA(ratioResult2, num, saveFile)
}
Cluster2("Data/Figure1D_FigureS2EFG/", 5000, "ESCA_hg38_MethylSeekR_30kb", "Figure1D/Figure1D_MethylSeekR.pdf")
Cluster2("Data/Figure1D_FigureS2EFG/", 5000, "ESCA_hg38_MethPipe_30kb", "Figure1D/Figure1D_MethPipe.pdf")
Cluster2("Data/Figure1D_FigureS2EFG/", 5000, "ESCA_hg38_MMSeekR_30kb", "Figure1D/Figure1D_MMSeekR.pdf")


###Figure S2HK
CalDistance.v2=function(PCAFile){
  PCAData=read_tsv(PCAFile)
  PCAData=as.data.frame(PCAData)
  PCAData=PCAData[grep("Nonmalignant", PCAData$Type, invert = T),]
  result=as.data.frame(matrix(numeric(0),ncol=3))
  for(type in unique(PCAData$Type)){
    intraSamples=PCAData[PCAData$Type%in%type,]
    av.PC1=mean(intraSamples$PC1)
    av.PC2=mean(intraSamples$PC2)
    intraSamples_distance=mean(do.call(c, lapply(1:nrow(intraSamples), function(x){((intraSamples[x,]$PC1-av.PC1)^2+(intraSamples[x,]$PC2-av.PC2)^2)^0.5})))
    interSamples=PCAData[!PCAData$Type%in%type,]
    interSamples_distance=mean(do.call(c, lapply(1:nrow(interSamples), function(x){((interSamples[x,]$PC1-av.PC1)^2+(interSamples[x,]$PC2-av.PC2)^2)^0.5})))
    result=rbind(result, data.frame(type=type, intraDis=intraSamples_distance, interDis=interSamples_distance))
  }
  result$ratio=result$interDis/result$intraDis
  return(result)
}
MethylSeekR_dis=CalDistance.v2("Figure1D/Figure1D_MethylSeekR.txt")
MethylSeekR_dis$group="MethylSeekR"
MethPipe_dis=CalDistance.v2("Figure1D/Figure1D_MethPipe.txt")
MethPipe_dis$group="MethPipe"
MMSeekR_dis=CalDistance.v2("Figure1D/Figure1D_MMSeekR.txt")
MMSeekR_dis$group="MMSeekR"
plotdata=rbind(MethylSeekR_dis, MethPipe_dis, MMSeekR_dis)
plotdata$type=gsub("_Tumor", "", plotdata$type)
plotdata[plotdata$type%in%"EAC/GEJ",]$type="EAC"
plotdata=as.data.frame(plotdata%>%group_by(group)%>%summarise(ratio=mean(ratio)))
plotdata$group=factor(plotdata$group, levels=c("MethylSeekR", "MethPipe", "MMSeekR"))
p=ggplot(data=plotdata, aes(x=group, y=ratio, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.6)+theme_classic()
p=p+xlab("")+ylab("Inter-tumor/Intra-tumor")+scale_fill_manual(values=c("blue", "#C900B8", "red"))
p=p+ggtitle("PCA distance")
p=p+theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size=10, color = "black"))

pdf("FigureS2HK/ESCA.PCADistance.pdf", width=4.3, height=3.5)
print(p)
dev.off()
write.table(plotdata, "FigureS2HK/ESCA.PCADistance.txt", row.names = F, col.names = T, sep="\t", quote=F)

CalDistance.v3=function(PCAFile){
  PCAData=read_tsv(PCAFile)
  PCAData=as.data.frame(PCAData)
  PCAData=PCAData[grep("Nonmalignant", PCAData$Subtype, invert = T),]
  result=as.data.frame(matrix(numeric(0),ncol=3))
  for(type in unique(PCAData$Subtype)){
    intraSamples=PCAData[PCAData$Subtype%in%type,]
    av.PC1=mean(intraSamples$PC1)
    av.PC2=mean(intraSamples$PC2)
    intraSamples_distance=mean(do.call(c, lapply(1:nrow(intraSamples), function(x){((intraSamples[x,]$PC1-av.PC1)^2+(intraSamples[x,]$PC2-av.PC2)^2)^0.5})))
    interSamples=PCAData[!PCAData$Subtype%in%type,]
    interSamples_distance=mean(do.call(c, lapply(1:nrow(interSamples), function(x){((interSamples[x,]$PC1-av.PC1)^2+(interSamples[x,]$PC2-av.PC2)^2)^0.5})))
    result=rbind(result, data.frame(type=type, intraDis=intraSamples_distance, interDis=interSamples_distance))
  }
  result$ratio=result$interDis/result$intraDis
  return(result)
}
MethylSeekR_dis=CalDistance.v3("FigureS2F/FigureS2F_MethylSeekR.txt")
MethylSeekR_dis$group="MethylSeekR"
MethPipe_dis=CalDistance.v3("FigureS2F/FigureS2F_Methpipe.txt")
MethPipe_dis$group="MethPipe"
MMSeekR_dis=CalDistance.v3("FigureS2F/FigureS2F_MMSeekR.txt")
MMSeekR_dis$group="MMSeekR"
plotdata=rbind(MethylSeekR_dis, MethPipe_dis, MMSeekR_dis)
plotdata=as.data.frame(plotdata%>%group_by(group)%>%summarise(ratio=mean(ratio)))
plotdata$group=factor(plotdata$group, levels=c("MethylSeekR", "MethPipe", "MMSeekR"))

p=ggplot(data=plotdata, aes(x=group, y=ratio, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.6)+theme_classic()
p=p+xlab("")+ylab("Inter-tumor/Intra-tumor")+scale_fill_manual(values=c("blue", "#C900B8", "red"))
p=p+ggtitle("PCA distance")
p=p+theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size=10, color = "black"))
pdf("FigureS2HK/BPTumor.PCADistance.pdf", width=4.3, height=3.5)
print(p)
dev.off()
write.table(plotdata, "FigureS2HK/BPTumor.PCADistance.txt", row.names = F, col.names = T, sep="\t", quote=F)


#######FigureS2I
hg38_regions_win2=bed_makewindows(hg38_regions, win_size=10000)
getPlotDataNew=function(RDataFile, pmdFileList, targetChr, sampleInfo, pattern, saveFile){
  load(RDataFile)
  varResult=rowVars(as.matrix(ratioResult2))
  clusterData=ratioResult2[order(varResult, decreasing = T), ][1:5000,]
  clusterData=t(clusterData)
  
  targetChrRegion=hg38_regions_win2[hg38_regions_win2$chrom%in%targetChr,]
  stateResult=as.data.frame(targetChrRegion)
  colnames(stateResult)=c("chr", "start", "end", "id")
  for(PMDFile in pmdFileList){
    PMDFileIndex=gsub(pattern, "", basename(PMDFile))
    print(PMDFileIndex)
    if(file.size(PMDFile)>0){
      PMDRegion=read_bed(PMDFile)
      PMDRegion=bed_merge(bed_sort(PMDRegion))
      PMDRegion=PMDRegion[PMDRegion$chrom%in%targetChr,]
      temp=bed_intersect(targetChrRegion, PMDRegion)
      temp=unique(temp[,c(1:3,7)])
      temp=as.data.frame(temp)
      colnames(temp)=c("chr", "start", "end", "overlap")
      temp$name=paste0(temp$chr, ":", temp$start)
      temp1=temp[!temp$name%in%temp[duplicated(temp$name),]$name,]
      temp2=temp[temp$name%in%temp[duplicated(temp$name),]$name,]
      temp2=do.call(rbind, lapply(unique(temp2$name),function(x){
        tmpData=temp2[temp2$name%in%x,][1,]
        tmpData$overlap=sum(temp2[temp2$name%in%x,]$overlap)
        return(tmpData)
      }))
      temp=rbind(temp1, temp2)
      temp=temp[,1:4]
      temp$ratio=temp$overlap/(temp$end-temp$start)
      temp$state=0
      temp[temp$ratio>0.5,]$state=1
      stateData=temp[,c(1:3,6)]
      colnames(stateData)[4]=PMDFileIndex
      stateResult=merge(stateResult, stateData, by.x=c("chr", "start", "end"), by.y=c("chr", "start", "end"), all.x=T)
    }else{
      ratioResult[[PMDFileIndex]]=NA
    }
  }
  rownames(stateResult)=paste0(stateResult$chr, ":",stateResult$start, "-",stateResult$end)
  stateResult=stateResult[,-1:-4]
  stateResult[is.na(stateResult)]=0
  stateResult=stateResult[,colnames(stateResult)%in%sampleInfo$Sample]
  stateResult=t(stateResult)
  
  targetSamples=sampleInfo[sampleInfo$Sample%in%rownames(stateResult),]
  annotation_row=targetSamples[,-1]
  rownames(annotation_row)=targetSamples$Sample
  
  annotation_row$commonPMDsMethylation=paste0("[", floor(annotation_row$PMD_methylation/0.05)*0.05, ", " , ceiling(annotation_row$PMD_methylation/0.05)*0.05, ")")
  commonPMDsColor=data.frame(commonPMDsMethylation=paste0("[",seq(0,0.95, by=0.05), ", ", seq(0.05,1, by=0.05), ")"),
                             colors=colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(20), stringsAsFactors = F)
  annotion_color=list(Type=as.character(unique(annotation_row[,c(1,2)])[,2]),
                      commonPMDsMethylation=commonPMDsColor$colors)
  names(annotion_color$Type)=as.character(unique(annotation_row[,c(1,2)])[,1])
  names(annotion_color$commonPMDsMethylation)=commonPMDsColor$commonPMDsMethylation
  annotation_row=annotation_row[,c("commonPMDsMethylation", "Type")]
  
  rowNames=factor(rownames(clusterData), levels=rownames(stateResult))
  clusterData=clusterData[order(rowNames),]
  drows = dist(clusterData)
  p=pheatmap(stateResult, clustering_distance_rows = drows, cluster_cols = F, show_rownames = F, show_colnames = F, legend=F,annotation_names_row = F,
             annotation_row = annotation_row, annotation_colors = annotion_color, color =c("#F7FBFF", "#2171B5"))
  stateResult=data.frame(sample=rownames(stateResult), stateResult)
  write.table(stateResult, file=gsub(".png", ".txt", saveFile), row.names = T, col.names = T, quote=F, sep="\t")
  return(p)
}
p1=getPlotDataNew("Data/Figure1D_FigureS2EFG/ESCA_hg38_MethylSeekR_30kb.RData", methylSeekRFileList, "chr16", ESCA_sampleInfo, ".MethylSeekR.rmBlackList.bed", "FigureS2I/FigureS2I_methylSeekR.chr16.png")
p2=getPlotDataNew("Data/Figure1D_FigureS2EFG/ESCA_hg38_MethPipe_30kb.RData", methPipeFileList, "chr16", ESCA_sampleInfo, ".MethPipe.rmBlackList.bed", "FigureS2I/FigureS2I_methpipe.chr16.png")
p3=getPlotDataNew("Data/Figure1D_FigureS2EFG/ESCA_hg38_MMSeekR_30kb.RData", model2D3DFileList, "chr16", ESCA_sampleInfo, ".Model2D3D.rmBlackList.bed", "FigureS2I/FigureS2I_MMSeekR.chr16.png")
png("FigureS2I/FigureS2I.png", width=3000, height=2100, res=300)
print(ggarrange(p3[[4]], p1[[4]], p2[[4]], ncol=1, align="hv"))
dev.off()


#######FigureS2G
hg19_regions=read_bed("meta/hg19/hg19.chrom.sizes")
blackListResgions=read_bed("meta/hg19/hg19-blacklist.v2.bed", n_fields = 4)
hg19_regions_rmBlackList=bed_subtract(hg19_regions, blackListResgions)
hg19_regions_rmBlackList=bed_merge(bed_sort(hg19_regions_rmBlackList))
hg19_regions_rmBlackList_win=bed_makewindows(hg19_regions_rmBlackList, win_size=30000)

bluePrintSamples=read.table("Data/Figure1D_FigureS2EFG/BLUEPRINT_blood_sampleInfo.txt", sep="\t", header=T, stringsAsFactors = F)
bluePrintSamples=bluePrintSamples[order(bluePrintSamples$Type),]
bluePrintSamples_tumor=bluePrintSamples[!bluePrintSamples$Type%in%c("Myeloid","Lymphoid","Others"),]
bluePrintSamples_tumor=bluePrintSamples_tumor[order(bluePrintSamples_tumor$Type),]
bluePrintSamples_tumor$Type=factor(bluePrintSamples_tumor$Type, levels=unique(bluePrintSamples_tumor$Type))
bluePrintSamples_tumor$Subtype=factor(bluePrintSamples_tumor$Subtype, levels=c("AcuteLymphocyticLeukemia", "TcellProlymphocyticLeukemia",
                                                                               "ChronicLymphocyticLeukemia", "MantleCellLymphoma", "MultipleMyeloma", "AcuteMyeloidLeukemia"))
colnames(bluePrintSamples_tumor)[ncol(bluePrintSamples_tumor)]="PMD_methylation"

methylSeekRFileList=dir("Data/MethylSeekR_PMDs/BLUEPRINT_Tumor/", pattern = ".PMDs.methylSeekR.bed", full.names = T)
methPipeFileList=dir("Data/MethPipe_PMDs/BLUEPRINT_Tumor/", pattern = ".methpipe.rmblaskList.bed", full.names = T)
MMSeekRFileList=dir("Data/MMSeekR_PMDs/BLUEPRINT_Tumor/", pattern = "_PMDs.rmblaskList.bed", full.names = T)

ClusterBPTumor=function(pmdFileList, num, fileIndex, sampleInfo, outputPath, pattern){
  ratioResult=as.data.frame(hg19_regions_rmBlackList_win)
  colnames(ratioResult)=c("chr", "start", "end", "id")
  stateResult=ratioResult
  overlapResult=c()
  for(PMDFile in pmdFileList){
    PMDFileIndex=gsub(pattern, "", basename(PMDFile))
    print(PMDFileIndex)
    if(file.size(PMDFile)>0){
      PMDRegion=read_bed(PMDFile)
      PMDRegion=bed_merge(bed_sort(PMDRegion))
      temp=bed_intersect(hg19_regions_rmBlackList_win, PMDRegion)
      temp=unique(temp[,c(1:3,7)])
      temp=as.data.frame(temp)
      colnames(temp)=c("chr", "start", "end", "overlap")
      overlapResult=c(overlapResult, sum(temp$overlap))
      
      temp$name=paste0(temp$chr, ":", temp$start)
      temp1=temp[!temp$name%in%temp[duplicated(temp$name),]$name,]
      temp2=temp[temp$name%in%temp[duplicated(temp$name),]$name,]
      temp2=do.call(rbind, lapply(unique(temp2$name),function(x){
        tmpData=temp2[temp2$name%in%x,][1,]
        tmpData$overlap=sum(temp2[temp2$name%in%x,]$overlap)
        return(tmpData)
      }))
      temp=rbind(temp1, temp2)
      temp=temp[,1:4]
      temp$ratio=temp$overlap/(temp$end-temp$start)
      temp$state=0
      temp[temp$ratio>0.5,]$state=1
      ratioData=temp[,c(1:3,5)]
      colnames(ratioData)[4]=PMDFileIndex
      stateData=temp[,c(1:3,6)]
      colnames(stateData)[4]=PMDFileIndex
      ratioResult=merge(ratioResult, ratioData, by.x=c("chr", "start", "end"), by.y=c("chr", "start", "end"), all.x=T)
      stateResult=merge(stateResult, stateData, by.x=c("chr", "start", "end"), by.y=c("chr", "start", "end"), all.x=T)
    }else{
      ratioResult[[PMDFileIndex]]=NA
      stateResult[[PMDFileIndex]]=NA
      overlapResult=c(overlapResult, 0)
    }
  }
  overlapResult=data.frame(overlap=overlapResult)
  overlapResult$Sample=colnames(stateResult)[-1:-4]
  overlapResult$overlap=overlapResult$overlap/sum(hg19_regions_rmBlackList$end-hg19_regions_rmBlackList$start)
  
  rownames(ratioResult)=paste0(ratioResult$chr, ":",ratioResult$start, "-",ratioResult$end)
  rownames(stateResult)=paste0(stateResult$chr, ":",stateResult$start, "-",stateResult$end)
  ratioResult=ratioResult[,-1:-4]
  stateResult=stateResult[,-1:-4]
  ratioResult[is.na(ratioResult)]=0
  stateResult[is.na(stateResult)]=0
  save(ratioResult, stateResult, overlapResult, file=paste0(outputPath, "/",fileIndex, ".RData"))
}
ClusterBPTumor(methylSeekRFileList, 5000, "BLUEPRINT_MethylSeekR_30kb", bluePrintSamples_tumor, "Data/Figure1D_FigureS2EFG/",".PMDs.methylSeekR.bed")
ClusterBPTumor(methPipeFileList, 5000, "BLUEPRINT_Methpipe_30kb", bluePrintSamples_tumor, "Data/Figure1D_FigureS2EFG/",".methpipe.rmblaskList.bed")
ClusterBPTumor(MMSeekRFileList, 5000, "BLUEPRINT_MMSeekR_30kb", bluePrintSamples_tumor, "Data/Figure1D_FigureS2EFG/","_PMDs.rmblaskList.bed")

plotBPTumorPCA=function(dataSet, sampleInfo, saveFile){
  varResult=rowVars(as.matrix(dataSet))
  clusterData=dataSet[order(varResult, decreasing = T), ][1:5000,]
  targetSamples=sampleInfo[sampleInfo$Sample%in%colnames(dataSet),]
  targetSamples$Type=factor(targetSamples$Type, levels=unique(targetSamples$Type))
  targetSamples$Subtype=factor(targetSamples$Subtype, levels=unique(targetSamples$Subtype))
  targetSamples$Color=factor(targetSamples$Color2, levels=unique(targetSamples$Color2))
  targetSamples$Sample=factor(targetSamples$Sample, levels=colnames(dataSet))
  targetSamples=targetSamples[order(targetSamples$Sample),]
  colorsData=targetSamples
  rownames(colorsData)=colorsData$Sample
  colorsData=colorsData[,-1]
  
  plotdata=t(clusterData)
  pca <- prcomp(plotdata,scale = TRUE)
  xlab <- paste("PC1","(",round((summary(pca))$importance[2,1]*100,1),"%)",sep="")
  ylab <- paste("PC2","(",round((summary(pca))$importance[2,2]*100,1),"%)",sep="")
  x<-"PC1"
  y<-"PC2"
  data_x <- data.frame(varnames=rownames(plotdata), pca$x[,1:2])
  data_x <- cbind(data_x, colorsData)
  my_theme=theme_bw()+theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                            panel.grid.major=element_blank(),axis.title=element_text(color="black",size=12),
                            axis.text=element_text(size=12), legend.title = element_text(hjust=0.5, face="bold"))
  write.table(data_x, gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote = F)
  p1 <- ggplot(data_x, aes(PC1,PC2))+geom_point(aes(color=Subtype, shape=Type),size=3, alpha=1)+coord_equal(ratio=1)+xlab(xlab)+ylab(ylab)
  p1 <- p1+scale_color_manual(values = levels(targetSamples$Color))+scale_shape_manual(values =c(15,17, 18))+my_theme
  pdf(saveFile, width = 5.5,height =4)
  print(p1)
  dev.off()
}
ClusterBPTumor2=function(fileIndex, sampleInfo, RDataPath, saveFile){
  load(paste0(RDataPath, "/",fileIndex, ".RData"))
  plotBPTumorPCA(ratioResult, sampleInfo, saveFile)
}
ClusterBPTumor2("BLUEPRINT_MMSeekR_30kb", bluePrintSamples_tumor, "Data/Figure1D_FigureS2EFG/", "FigureS2G/FigureS2G_MMSeekR.pdf")
ClusterBPTumor2("BLUEPRINT_Methpipe_30kb", bluePrintSamples_tumor, "Data/Figure1D_FigureS2EFG/", "FigureS2G/FigureS2G_Methpipe.pdf")
ClusterBPTumor2("BLUEPRINT_MethylSeekR_30kb", bluePrintSamples_tumor, "Data/Figure1D_FigureS2EFG/", "FigureS2G/FigureS2G_MethylSeekR.pdf")


#######FigureS2E
hg19_regions_win2=bed_makewindows(hg19_regions, win_size=10000)
getPlotDataNewBPTumor=function(RDataFile, pmdFileList, targetChr, sampleInfo, pattern, saveFile){
  load(RDataFile)
  varResult=rowVars(as.matrix(ratioResult))
  clusterData=ratioResult[order(varResult, decreasing = T), ][1:5000,]
  clusterData=t(clusterData)
  
  targetChrRegion=hg19_regions_win2[hg19_regions_win2$chrom%in%targetChr,]
  stateResult=as.data.frame(targetChrRegion)
  colnames(stateResult)=c("chr", "start", "end", "id")
  for(PMDFile in pmdFileList){
    PMDFileIndex=gsub(pattern, "", basename(PMDFile))
    print(PMDFileIndex)
    if(file.size(PMDFile)>0){
      PMDRegion=read_bed(PMDFile)
      PMDRegion=bed_merge(bed_sort(PMDRegion))
      PMDRegion=PMDRegion[PMDRegion$chrom%in%targetChr,]
      temp=bed_intersect(targetChrRegion, PMDRegion)
      temp=unique(temp[,c(1:3,7)])
      temp=as.data.frame(temp)
      colnames(temp)=c("chr", "start", "end", "overlap")
      temp$name=paste0(temp$chr, ":", temp$start)
      temp1=temp[!temp$name%in%temp[duplicated(temp$name),]$name,]
      temp2=temp[temp$name%in%temp[duplicated(temp$name),]$name,]
      temp2=do.call(rbind, lapply(unique(temp2$name),function(x){
        tmpData=temp2[temp2$name%in%x,][1,]
        tmpData$overlap=sum(temp2[temp2$name%in%x,]$overlap)
        return(tmpData)
      }))
      temp=rbind(temp1, temp2)
      temp=temp[,1:4]
      temp$ratio=temp$overlap/(temp$end-temp$start)
      temp$state=0
      temp[temp$ratio>0.5,]$state=1
      stateData=temp[,c(1:3,6)]
      colnames(stateData)[4]=PMDFileIndex
      stateResult=merge(stateResult, stateData, by.x=c("chr", "start", "end"), by.y=c("chr", "start", "end"), all.x=T)
    }else{
      ratioResult[[PMDFileIndex]]=NA
    }
  }
  rownames(stateResult)=paste0(stateResult$chr, ":",stateResult$start, "-",stateResult$end)
  stateResult=stateResult[,-1:-4]
  stateResult[is.na(stateResult)]=0
  stateResult=stateResult[,colnames(stateResult)%in%sampleInfo$Sample]
  stateResult=t(stateResult)
  
  targetSamples=sampleInfo[sampleInfo$Sample%in%rownames(stateResult),]
  annotation_row=targetSamples[,-1]
  annotation_row$commonPMDsMethylation=paste0("[", floor(annotation_row$PMD_methylation/0.05)*0.05, ", " , ceiling(annotation_row$PMD_methylation/0.05)*0.05, ")")
  commonPMDsColor=data.frame(commonPMDsMethylation=paste0("[",seq(0,0.95, by=0.05), ", ", seq(0.05,1, by=0.05), ")"),
                             colors=colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(20), stringsAsFactors = F)
  
  rownames(annotation_row)=targetSamples$Sample
  annotion_color=list(Type=as.character(unique(annotation_row[,c(1,3)])[,2]), 
                      Subtype=as.character(unique(annotation_row[,c(2,4)])[,2]),
                      commonPMDsMethylation=commonPMDsColor$colors)
  names(annotion_color$Type)=as.character(unique(annotation_row[,c(1,3)])[,1])
  names(annotion_color$Subtype)=as.character(unique(annotation_row[,c(2,4)])[,1])
  names(annotion_color$commonPMDsMethylation)=commonPMDsColor$commonPMDsMethylation
  annotation_row=annotation_row[,c(6,2,1)]
  
  rowNames=factor(rownames(clusterData), levels=rownames(stateResult))
  clusterData=clusterData[order(rowNames),]
  drows = dist(clusterData)
  
  callback = function(hc, mat){
    sv = svd(t(mat))$v[,2]
    dend = reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
  }
  
  p=pheatmap(stateResult, clustering_distance_rows = drows, cluster_cols = F, clustering_callback = callback, show_rownames = F, show_colnames = F, legend=F,annotation_names_row = F,
             annotation_row = annotation_row, annotation_colors = annotion_color, color =c("#F7FBFF", "#2171B5"))
  write.table(stateResult, file=gsub(".png", ".txt", saveFile), row.names = T, col.names = T, quote=F, sep="\t")
  return(p)
}
p1=getPlotDataNewBPTumor("Data/Figure1D_FigureS2EFG/BLUEPRINT_MethylSeekR_30kb.RData", methylSeekRFileList, "chr16", bluePrintSamples_tumor, ".PMDs.methylSeekR.bed", "FigureS2E/FigureS2E_methylSeekR.chr16.png")
p2=getPlotDataNewBPTumor("Data/Figure1D_FigureS2EFG/BLUEPRINT_Methpipe_30kb.RData", methPipeFileList, "chr16", bluePrintSamples_tumor, ".methpipe.rmblaskList.bed", "FigureS2E/FigureS2E_Methpipe.chr16.png")
p3=getPlotDataNewBPTumor("Data/Figure1D_FigureS2EFG/BLUEPRINT_MMSeekR_30kb.RData", MMSeekRFileList, "chr16", bluePrintSamples_tumor, "_PMDs.rmblaskList.bed", "FigureS2E/FigureS2E_MMSeekR.chr16.png")

png("FigureS2E/FigureS2E.png", width=3000, height=2100, res=300)
print(ggarrange(p3[[4]], p1[[4]], p2[[4]], ncol=1, align="hv"))
dev.off()

################################################FigureS1B################################################F
load("Data/FigureS1B/mergeESCCNonmalignantCov.RData")
load("Data/FigureS1B/mergeESCCNonmalignantBeta.RData")

##calculate the correlation
cor_result=as.data.frame(matrix(numeric(0), ncol=3))
for(i in 4:(8-1)){
  for(j in (i+1):8){
    group1=colnames(ESCCNonmalignantCovMatrix)[i]
    group2=colnames(ESCCNonmalignantCovMatrix)[j]
    if(group1!=group2){
      targetCovMatrix=ESCCNonmalignantCovMatrix[,colnames(ESCCNonmalignantCovMatrix)%in%c(group1, group2, "region")]
      targetRegion=targetCovMatrix[rowSums(targetCovMatrix[,-1]>=cut1&targetCovMatrix[,-1]<cut2)==2,]$region
      plotdata=ESCCNonmalignantBetaMatrix[,colnames(ESCCNonmalignantBetaMatrix)%in%c(group1, group2, "region")]
      plotdata=plotdata[plotdata$region%in%targetRegion,]
      cor=cor(plotdata[[group1]],plotdata[[group2]])
      cor_result=rbind(cor_result, data.frame(group1=group1, group2=group2, cor= cor))
    }
  }
}
write.table(cor_result, file="FigureS1B/cor_result.txt", row.names = F, col.names = T, sep="\t", quote=F)

##change FigureS1B/cor_result.txt to FigureS1B/FigureS1B.txt
data=read.table("FigureS1B/FigureS1B.txt", sep="\t", row.names = 1, header = T, stringsAsFactors = F)
breaks=seq(0,1,by=0.01)
colorList=colorRampPalette(brewer.pal(9, "PuRd"))(length(breaks))
pdf("FigureS1B/FigureS1B.pdf")
print(pheatmap(data, breaks= breaks, color=colorList, cluster_rows = F, cluster_cols = F,display_numbers=T,fontsize = 14))
dev.off()

################################################FigureS3A################################################
totalGenomeLength=read_bed("meta/hg38_length.bed")
blackList=read_bed("meta/hg38-blacklist.v2.change.bed",n_fields = 4)
totalGenomeLength=as.data.frame(bed_subtract(totalGenomeLength, blackList))
totalGenomeLength$length=totalGenomeLength$end-totalGenomeLength$start
totalGenomeLength=sum(totalGenomeLength$length)

plotGenomeCov=function(PMDfilePath, pattern, saveFile){
  fileList=dir(PMDfilePath, pattern=pattern, full.names = T)
  fileList=fileList[grep("Nonmalignant", fileList, invert = T)]
  pmd_result=as.data.frame(matrix(numeric(0),ncol=2))
  for(file in fileList){
    fileIndex=gsub(pattern, "", basename(file))
    print(fileIndex)
    data=read.table(file, sep="\t", stringsAsFactors = F)
    data$length=data$V3-data$V2
    totalLength=sum(data$length)
    rate=round(totalLength/totalGenomeLength,4)
    pmd_result=rbind(pmd_result, data.frame(sample=fileIndex, PMD_fraction=rate, stringsAsFactors = F))
  }
  pmd_result$Type=""
  pmd_result[pmd_result$sample%in%EAC_Tumor_sampleList,]$Type="EAC/GEJ tumor"
  pmd_result[pmd_result$sample%in%ESCC_Tumor_sampleList,]$Type="ESCC tumor"
  pmd_result$Type=factor(pmd_result$Type, levels=c("ESCC tumor", "EAC/GEJ tumor"))
  pmd_result=pmd_result[order(pmd_result$Type),]
  
  pmd_result$PMD_fraction=pmd_result$PMD_fraction*100
  write.table(pmd_result, file=gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
  p=ggplot2::qplot(Type, PMD_fraction, data = pmd_result) + geom_beeswarm(priority='density',size=2)+theme_classic()
  p=p+xlab("")+ylab("Fraction of genome \ncovered by PMDs (%)")
  p=p+theme(plot.title = element_text(hjust = 0.5, size=14, color="black", face="bold"),
            axis.text = element_text(size=14, color="black"),
            axis.title = element_text(size=16, color="black", face="bold"),
            legend.title =element_text(size=13, color="black", face="bold"),
            legend.text =element_text(size=13, color="black"),
            legend.position = "none")
  p=p+scale_y_continuous(limits = c(0, 75), labels = seq(0, 75, 25), breaks = seq(0, 75, 25))
  tTest=t.test(pmd_result[pmd_result$Type%in%"ESCC tumor",]$PMD_fraction, pmd_result[pmd_result$Type%in%"EAC/GEJ tumor",]$PMD_fraction)
  print(tTest$p.value)
  pdf(saveFile, width=4, height=5)
  print(p)
  dev.off()
}
plotGenomeCov("Data/MMSeekR_PMDs/PMDs/", ".Model2D3D.rmBlackList.bed", "FigureS3A/FigureS3A.pdf")


################################################FigureS3BC################################################
##bash Shell/FigureS3BC/runMerge.sh
plotCumSumPlot=function(PMDfile,title, width, saveFile){
  data=read.table(PMDfile,sep="\t", stringsAsFactors = F,header=T)
  data=data[1:4]
  data$length=data$end-data$start
  write.table(data, file=gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
  plotdata=data.frame(num= unique(data$num), length=do.call(c, lapply(unique(data$num), function(x){sum(data[data$num==x,]$length)})))
  plotdata=plotdata[order(plotdata$num),]
  plotdata$temp=c(0,cumsum(as.numeric(plotdata$length)))[1:nrow(plotdata)]
  plotdata$cumsum=sum(as.numeric(plotdata$length))-plotdata$temp
  p<-ggplot(data=plotdata, aes(x=num, y=cumsum)) + geom_bar(stat="identity", fill="black")+ggtitle(title)
  p=p+scale_x_continuous(name = ">= # samples", breaks = unique(plotdata$num), labels = unique(plotdata$num))
  p=p+scale_y_continuous(name="# bps in PMDs (MB)", limits=c(0, 2000000000), breaks = seq(0,2000000000,400000000), labels = paste0(format(seq(0,2000,400),big.mark=",",scientific=FALSE)))+theme_classic()
  p = p+theme(axis.title = element_text(size=14, face="bold", color = "black"), axis.text.y = element_text(size=12, color = "black"), axis.text.x = element_text(size=9, color = "black"), plot.title = element_text(size=16, hjust = 0.5, face="bold", color = "black"))
  pdf(saveFile, width=width, height=5)
  print(p)
  dev.off()
}
plotCumSumPlot("Data/MMSeekR_PMDs/ESCC_PMDs_multiinter.bed", "ESCC Tumor", 8, "FigureS3BC/FigureS3B.pdf")
plotCumSumPlot("Data/MMSeekR_PMDs/EAC_PMDs_multiinter.bed", "EAC/GEJ Tumor", 8, "FigureS3BC/FigureS3C.pdf")


################################################FigureS4A################################################
load("Data/FigureS4A/AllSample_beta_cov7.RData")
plotPCA=function(plotdata,annotation_col, colorInfo1, colorInfo2, saveFile){
  plotdata=t(plotdata)
  pca <- prcomp(plotdata,scale = TRUE)
  xlab <- paste("PC1","(",round((summary(pca))$importance[2,1]*100,1),"%)",sep="")
  ylab <- paste("PC2","(",round((summary(pca))$importance[2,2]*100,1),"%)",sep="")
  x<-"PC1"
  y<-"PC2"
  data_x <- data.frame(varnames=rownames(plotdata), pca$x)
  data_x=data_x[,colnames(data_x)%in%c("varnames", "PC1", "PC2")]
  data_x <- merge(data_x, annotation_col, by.x="varnames", by.y="Sample")
  data_x=merge(data_x, colorInfo1[,c(1,3)], by.x="varnames", by.y="Sample")
  data_x=merge(data_x, colorInfo2[,c(1,3)], by.x="varnames", by.y="Sample")
  colnames(data_x)[5:6]=c("global", "commonPMDs")
  my_theme=theme_bw()+ theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                             panel.grid.major=element_blank(),axis.title=element_text(color="black",size=12),
                             axis.text=element_text(size=12), legend.title = element_text(face="bold"))

  write.table(data_x, file=gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
  p1 <- ggplot(data_x, aes(PC1,PC2, label = varnames,color=Type, fill=Type))+geom_point(size=4, alpha=1, pch=21)+coord_equal(ratio=1)+xlab(xlab)+ylab(ylab)
  p1=p1+scale_color_manual(values=c("#EA3323", "#0000F5", "#EA3323", "#0000F5"))+scale_fill_manual(values=ann_colors$Type[1:4])+my_theme
  
  p3 <- ggplot(data_x, aes(PC1,PC2, label = varnames,color=commonPMDs))+geom_point(size=4, alpha=1)+coord_equal(ratio=1)+xlab(xlab)+ylab(ylab)
  p3= p3+scale_color_gradientn(name="commonPMDs methylation", colours = c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"), limits=c(0,1))+my_theme
  

  pdf(saveFile,width = 12,height =5)
  print(ggarrange(p1, p3, align = "hv", nrow=1))
  dev.off()
}
getPCAPlot=function(plotdata, saveFile, num=8000){
  rowVar=rowVars(as.matrix(plotdata[,-1:-3]))
  plotdata=plotdata[order(rowVar,decreasing = T),]
  if(nrow(plotdata)>8000){
    plotdata=plotdata[1:8000,-1:-3]
  }else{
    plotdata=plotdata[1:num,-1:-3]
  }
  plotPCA(plotdata,annotation_col, colorInfo1, colorInfo2, saveFile)
}

allMethylationInfo=read.table("Figure1B/Methylaton_sixTypes.txt", sep="\t", stringsAsFactors = F, header=T)
colorInfo1=allMethylationInfo[,colnames(allMethylationInfo)%in%c("Sample", "Type", "Gobal")]
colorInfo2=allMethylationInfo[,colnames(allMethylationInfo)%in%c("Sample", "Type", "commonPMD")]
getPCAPlot(combine_betaMatrix, "FigureS4A/FigureS4A.pdf")

################################################FigureS4B################################################
totalGenomeLength=read_bed("meta/hg38_length.bed")
blackList=read_bed("meta/hg38-blacklist.v2.change.bed",n_fields = 4)
totalGenomeLength=as.data.frame(bed_subtract(totalGenomeLength, blackList))
totalGenomeLength$length=totalGenomeLength$end-totalGenomeLength$start
totalGenomeLength=sum(totalGenomeLength$length)

vennTwogroup=function(Agroup_Total, Bgroup_Total, AB_intersect, Agroup_name, Bgroup_name, Agroup_color, Bgroup_color, saveFile){
  Agroup=seq(1,Agroup_Total, by=1)
  Bgroup=c(seq(1,AB_intersect, by=1), seq(Agroup_Total+1, Bgroup_Total+Agroup_Total-AB_intersect))
  print(length(unique(c(Agroup, Bgroup))))
  print(length(intersect(Agroup, Bgroup)))
  venn.plot <- venn.diagram(
    x = list(Agroup, Bgroup),category.names=c(Agroup_name, Bgroup_name),
    filename = NULL,
    col = "transparent",
    fill = c(Agroup_color, Bgroup_color),
    label.col = c("black", "black", "black"),
    cex = 1.5,
    alpha=0.8,
    fontface = "bold",
    height = 1500,
    width = 1500,
    resolution = 1500,
    cat.col = c(Agroup_color, Bgroup_color),
    cat.cex = 1.2,
    cat.pos = 0,
    cat.dist = 0.07,
    cat.fontfamily = "serif",
    margin = 0.02
  )
  pdf(saveFile)
  grid.draw(venn.plot)
  dev.off()
}
maskedPMDRegions=function(EACPMDFile, ESCCPMDFile, outputFile, saveVennPlot){
  EAC_PMDs=read_bed(EACPMDFile)
  EAC_PMDs_length=sum(EAC_PMDs$end-EAC_PMDs$start)
  EAC_PMDs_length_ratio=EAC_PMDs_length/totalGenomeLength
  print(paste0(EAC_PMDs_length, "    ", EAC_PMDs_length_ratio))
  
  ESCC_PMDs=read_bed(ESCCPMDFile)
  ESCC_PMDs_length=sum(ESCC_PMDs$end-ESCC_PMDs$start)
  ESCC_PMDs_length_ratio=ESCC_PMDs_length/totalGenomeLength
  print(paste0(ESCC_PMDs_length, "    ",ESCC_PMDs_length_ratio))
  
  union_PMDs=rbind(EAC_PMDs, ESCC_PMDs)
  union_PMDs=bed_sort(union_PMDs)
  union_PMDs=bed_merge(union_PMDs)
  union_PMDs=as.data.frame(union_PMDs)
  write.table(union_PMDs, file=outputFile, row.names = F, col.names = F, sep="\t", quote=F)
  length=sum(union_PMDs$end-union_PMDs$start)
  ratio=length/totalGenomeLength
  print(paste0(length, "    ", ratio))
  
  intersectLength=ESCC_PMDs_length-(length-EAC_PMDs_length)
  vennTwogroup(round(EAC_PMDs_length/1000000), round(ESCC_PMDs_length/1000000), round(intersectLength/1000000),
               "EAC_commonPMDs", "ESCC_commonPMDs", "#0000F5", "#EA3323", saveVennPlot)
}
maskedPMDRegions("Data/MMSeekR_PMDs/EAC_commonPMDs.bed", "Data/MMSeekR_PMDs/ESCC_commonPMDs.bed", "Data/MMSeekR_PMDs/ESCA_commonPMDs_union.bed",
                 "FigureS4B/FigureS4B.pdf")

################################################Figure2AB, FigureS1A and Figure4B################################################
##Figure2A 
##Prepare the data 
#bash Shell/Figure2A/get_Hg38_10kb_meanMethylation.sh
get_hg38_10k=function(inputPath, suffix){
  hg38_10k=read.table("meta/hg38_10k.bed", sep="\t", stringsAsFactors = F)
  colnames(hg38_10k)=c("chr", "start", "end")
  hg38_10k$region=paste0(hg38_10k$chr, "_", hg38_10k$start, "_", hg38_10k$end)
  for(file in fileList){
    print(file)
    data=read.table(paste0(inputPath, "/","hg38_10k.", file, suffix,".rmCGI.txt"), sep="\t", stringsAsFactors = F)
    colnames(data)=c("region", paste0(file, "_Mean"))
    hg38_10k=merge(hg38_10k, data, by.x="region", by.y ="region", all.x=T)
  }
  save(hg38_10k, file=paste0(inputPath, "/", "hg38_10k", suffix, ".rmCGI.RData"))
}
get_hg38_10k("Data/Figure2A/", "")

load("Data/Figure2A/hg38_10k.rmblackListCGI.RData")
hg38_10k$chr=factor(hg38_10k$chr, levels=paste0("chr", 1:22))
hg38_10k=hg38_10k[order(hg38_10k$chr, hg38_10k$start),]
hg38_10k_mean=hg38_10k[,-1:-4]
rownames(hg38_10k_mean)=hg38_10k$region
hg38_10k_mean_temp=hg38_10k_mean[rowSums(is.na(hg38_10k_mean))==0,]
chr1_methy=hg38_10k_mean_temp[grep("chr1_",rownames(hg38_10k_mean_temp)),]
ESCC_methy=chr1_methy[,colnames(chr1_methy)%in%paste0("ESCC_", c(1:17,19:22))]
ESCC_methy_sample=colnames(ESCC_methy)[hclust(dist(t(ESCC_methy)))$order]
ESCC_nonmalignant_methy=chr1_methy[,colnames(chr1_methy)%in%paste0("ESCC_Nonmalignant_", c(1:3, "53F", "54M"))]
ESCC_nonmalignant_methy_sample=colnames(ESCC_nonmalignant_methy)[hclust(dist(t(ESCC_nonmalignant_methy)))$order]
EAC_methy=chr1_methy[,colnames(chr1_methy)%in%c(paste0("EAC_", c(1:4,6)),paste0("GEJ_", 1:7))]
EAC_methy_sample=colnames(EAC_methy)[hclust(dist(t(EAC_methy)))$order]
EAC_nonmalignant_methy=chr1_methy[,colnames(chr1_methy)%in%paste0("GEJ_Nonmalignant_", 1:7)]
EAC_nonmalignant_methy_sample=colnames(EAC_nonmalignant_methy)[hclust(dist(t(EAC_nonmalignant_methy)))$order]
colnamesOrder=factor(colnames(hg38_10k_mean), levels=c(ESCC_nonmalignant_methy_sample, EAC_nonmalignant_methy_sample, ESCC_methy_sample, EAC_methy_sample))
hg38_10k_mean=hg38_10k_mean[,order(colnamesOrder)]
hg38_10k_mean=as.data.frame(t(hg38_10k_mean))

for(chrInfo in paste0("chr", 1:22, "_")){
  methyBreaksList=seq(0, 0.8, by = 0.01)
  print(chrInfo)
  chr_methy=hg38_10k_mean[,grep(chrInfo,colnames(hg38_10k_mean))]
  annotation_row2=annotation_row[rownames(annotation_row)%in%rownames(chr_methy),,drop=F]
  if(chrInfo=="chr1_"){
    pdf(paste0("FigureS1A/", gsub("_", "", chrInfo),"_test.pdf"), width=7, height=6)
    print(pheatmap(chr_methy[,1:500], breaks = methyBreaksList, show_colnames = F,cluster_rows = F, cluster_cols = F, annotation_row=annotation_row2, annotation_colors = ann_colors,
                   color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(n = length(methyBreaksList))))
    dev.off()
  }
  png(paste0("FigureS1A/", gsub("_", "", chrInfo), ".png"), res=300, width=1500, height=1500)
  print(pheatmap(chr_methy, breaks = methyBreaksList, show_colnames = F,cluster_rows = F, cluster_cols = F, annotation_row=annotation_row2, annotation_colors = ann_colors,fontsize_row=8,
                 color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(n = length(methyBreaksList))))
  dev.off()
  if(chrInfo=="chr3_"){
    chr_methy1=chr_methy[,c(grep("chr3_102000000_102010000",colnames(chr_methy)):grep("chr3_119990000_120000000",colnames(chr_methy)))]
    pdf("Figure2A/Figure2A.pdf",width=8, height=4)
    print(pheatmap(chr_methy1, breaks = methyBreaksList, show_colnames = F,cluster_rows = F, cluster_cols = F, 
                   annotation_row=annotation_row, annotation_colors = ann_colors,fontsize_row=8, border_color = NA,
                   color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(n = length(methyBreaksList))))
    dev.off()
    chr_methy1=t(chr_methy1)
    chr_methy1=data.frame(region=rownames(chr_methy1), chr_methy1)
    write.table(chr_methy1, "Figure2A/Figure2A.txt", row.names = F, col.names = T, sep="\t", quote=F)
    
    chr_methy1=chr_methy[,c(grep("chr3_185040000_185050000",colnames(chr_methy)):grep("chr3_185140000_185150000",colnames(chr_methy)))]
    pdf("Figure4B/Figure4B_chr3.pdf",width=8, height=4)
    print(pheatmap(chr_methy1, breaks = methyBreaksList, show_colnames = F,cluster_rows = F, cluster_cols = F, 
                   annotation_row=annotation_row, annotation_colors = ann_colors,fontsize_row=8, border_color = NA,
                   color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(n = length(methyBreaksList))))
    dev.off()
    chr_methy1=t(chr_methy1)
    chr_methy1=data.frame(region=rownames(chr_methy1), chr_methy1)
    write.table(chr_methy1, "Figure4B/Figure4B_chr3.txt", row.names = F, col.names = T, sep="\t", quote=F)
  }
  
  if(chrInfo=="chr12_"){
    chr_methy1=chr_methy[,c(grep("chr12_57880000_57890000",colnames(chr_methy)):grep("chr12_58240000_58250000",colnames(chr_methy)))]
    pdf("Figure4B/Figure4B_chr12.pdf",width=8, height=4)
    print(pheatmap(chr_methy1, breaks = methyBreaksList, show_colnames = F,cluster_rows = F, cluster_cols = F, 
                   annotation_row=annotation_row, annotation_colors = ann_colors,fontsize_row=8, border_color = NA,
                   color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(n = length(methyBreaksList))))
    dev.off()
    chr_methy1=t(chr_methy1)
    chr_methy1=data.frame(region=rownames(chr_methy1), chr_methy1)
    write.table(chr_methy1, "Figure4B/Figure4B_chr12.txt", row.names = F, col.names = T, sep="\t", quote=F)
  }
}

##Figure2B
getPMDForIGV=function(multiinterFile, saveIndex){
  totalGenomeLength=read_bed("meta/hg38_length.bed")
  blackList=read_bed("meta/hg38-blacklist.v2.change.bed",n_fields = 4)
  otherGenome=bed_subtract(totalGenomeLength, blackList)
  multiinterPMDsAll=read.table(multiinterFile,sep="\t", stringsAsFactors = F,header=T)
  sampleLength=length(colnames(multiinterPMDsAll)[-1:-5])
  multiinterPMDsAll_region=tibble(chrom=multiinterPMDsAll$chrom, start=as.integer(multiinterPMDsAll$start), 
                                  end=as.integer(multiinterPMDsAll$end))
  otherGenome=bed_subtract(otherGenome, multiinterPMDsAll_region)
  
  multiinterPMDs=multiinterPMDsAll[,1:4]
  multiinterPMDs1=multiinterPMDs[multiinterPMDs$num>=(2/3*sampleLength),]
  multiinterPMDs1=tibble(chrom=multiinterPMDs1$chrom, start=multiinterPMDs1$start, end=multiinterPMDs1$end)
  multiinterPMDs1=bed_sort(multiinterPMDs1)
  multiinterPMDs1=bed_merge(multiinterPMDs1)
  write.table(multiinterPMDs1, paste0(saveIndex, "_2_3.bed"), row.names = F, col.names = F, sep="\t", quote = F)
  
  multiinterPMDs2=multiinterPMDs[(multiinterPMDs$num>=(1/3*sampleLength))&(multiinterPMDs$num<(2/3*sampleLength)),]
  multiinterPMDs2=tibble(chrom=multiinterPMDs2$chrom, start=multiinterPMDs2$start, end=multiinterPMDs2$end)
  multiinterPMDs2=bed_sort(multiinterPMDs2)
  multiinterPMDs2=bed_merge(multiinterPMDs2)
  write.table(multiinterPMDs2, paste0(saveIndex, "_1_3.bed"), row.names = F, col.names = F, sep="\t", quote = F)
  
  multiinterPMDs3=multiinterPMDs[multiinterPMDs$num<(1/3*sampleLength),]
  multiinterPMDs3=tibble(chrom=multiinterPMDs3$chrom, start=multiinterPMDs3$start, end=multiinterPMDs3$end)
  multiinterPMDs3=rbind(multiinterPMDs3, otherGenome)
  multiinterPMDs3=bed_sort(multiinterPMDs3)
  multiinterPMDs3=bed_merge(multiinterPMDs3)
  write.table(multiinterPMDs3, paste0(saveIndex, "_others.bed"), row.names = F, col.names = F, sep="\t", quote = F)
  
}
getPMDForIGV("Data/MMSeekR_PMDs/ESCC_PMDs_multiinter.bed", "Figure2B/ESCC_PMDs")
getPMDForIGV("Data/MMSeekR_PMDs/EAC_PMDs_multiinter.bed", "Figure2B/EAC_PMDs")


################################################Figure2C_7A################################################
##get region level PMD methylation
#./Shell/Figure2C/getRegionMeanBetaValues.solo.sh Data/MMSeekR_PMDs/EAC_specificPMDs.bed
#./Shell/Figure2C/getRegionMeanBetaValues.solo.sh Data/MMSeekR_PMDs/ESCC_specificPMDs.bed
#./Shell/Figure2C/getRegionMeanBetaValues.solo.sh Data/MMSeekR_PMDs/ESCA_sharedPMDs.bed
#./Shell/Figure2C/getRegionMeanBetaValues.solo.sh Data/MMSeekR_PMDs/ESCA_sharedHMDs.bed
##merge data into regionLevels_methylation.soloCpGs.txt

plotRegionLineplot=function(plotdata, type, ylab){
  plotdata$RegionType=factor(plotdata$RegionType, levels=c("SharedPMDs", "EAC_specificPMDs", "ESCC_specificPMDs", "SharedHMDs"))
  p=ggplot(data=plotdata, aes(x=RegionType, y=Methylation, group=Sample)) + geom_line(alpha=0.2, color="darkblue")+ geom_point(alpha=0.5, color="darkblue", size=0.8)+theme_classic()
  p=p+xlab("")+ylab(ylab)+ggtitle(type)+ylim(0, 1)
  p=p+theme(plot.title = element_text(hjust = 0.5, size=14, color="black", face="bold"),
            axis.text.x = element_text(size=11, color="black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size=11, color="black", angle = 0),
            axis.title = element_text(size=13, color="black", face="bold"),
            legend.title =element_text(size=12, color="black", face="bold"),
            legend.text =element_text(size=12, color="black"),legend.position = "none")
  return(p)
}
data=read.table("Figure2C_7A/regionLevels_methylation.soloCpGs.txt", sep="\t", stringsAsFactors = F, header=T, check.names = F)
data=melt(data, id.vars = c("Sample", "Type"))
colnames(data)=c("Sample", "SampleType", "RegionType", "Methylation")
p1=plotRegionLineplot(data[data$SampleType%in%"EAC/GEJ_Tumor",], "EAC/GEJ_Tumor", "Domain mean methylation")
p2=plotRegionLineplot(data[data$SampleType%in%"ESCC_Tumor",], "ESCC_Tumor", "Domain mean methylation")
p3=plotRegionLineplot(data[data$SampleType%in%"GEJ_Nonmalignant",], "GEJ_Nonmalignant", "Domain mean methylation")
p4=plotRegionLineplot(data[data$SampleType%in%"ESCC_Nonmalignant",], "ESCC_Nonmalignant", "Domain mean methylation")
pdf("Figure2C_7A/Figure2C.pdf", width=7, height=4)
print(ggarrange(p1,p2, nrow=1))
dev.off()

pdf("Figure2C_7A/Figure7A.pdf", width=7, height=4)
print(ggarrange(p3,p4, nrow=1))
dev.off()

################################################Figure2D################################################
##TCGA data
load("meta/TCGA-ESCA.clinc.RData")
load("meta/TCGA-ESCA_methylation_hg38.v2.RData")
TCGAprobeMethylationMatrix=result
TCGAprobeSampleInfo=match.file.cases
TCGAprobeSampleInfo$Subtype=""
TCGAprobeSampleInfo[TCGAprobeSampleInfo$Sample%in%sampleInformation[sampleInformation$Label%in%"EAC_Tumor",]$Sample,]$Subtype="EAC_Tumor"
TCGAprobeSampleInfo[TCGAprobeSampleInfo$Sample%in%sampleInformation[sampleInformation$Label%in%"EAC_Normal",]$Sample,]$Subtype="EAC_Normal"
TCGAprobeSampleInfo[TCGAprobeSampleInfo$Sample%in%sampleInformation[sampleInformation$Label%in%"ESCC_Tumor",]$Sample,]$Subtype="ESCC_Tumor"
TCGAprobeSampleInfo=TCGAprobeSampleInfo[,c(4,5)]
colnames(TCGAprobeSampleInfo)=c("Sample", "Type")
TCGAprobeSampleInfo=TCGAprobeSampleInfo[TCGAprobeSampleInfo$Type!="",]
getMethylationHT450k=function(regionFile, probesFile, type, methylationMatrix){
  tempData=read_bed(regionFile,n_fields = 5)
  probesMatrix=read_bed(probesFile, n_fields = 4)
  tempData_probes=bed_intersect(probesMatrix,tempData)
  tempData_probes=as.data.frame(tempData_probes)
  tempData_probes=unique(tempData_probes[tempData_probes$.overlap==2,]$name.x)
  print(length(tempData_probes))
  tempData=methylationMatrix[rownames(methylationMatrix)%in%tempData_probes,, drop=F]
  tempData=data.frame(colMeans(tempData, na.rm = T))
  colnames(tempData)=type
  return(tempData)
}

plotHM450kMethylationLinePlot=function(probesFile, methylationMatrix, sampleInfo, saveFile, height, width){
  temp1=getMethylationHT450k("Data/MMSeekR_PMDs/ESCC_specificPMDs.bed", probesFile, "ESCC_specificPMDs",  methylationMatrix)
  temp2=getMethylationHT450k("Data/MMSeekR_PMDs/EAC_specificPMDs.bed", probesFile, "EAC_specificPMDs", methylationMatrix)
  temp3=getMethylationHT450k("Data/MMSeekR_PMDs/ESCA_sharedPMDs.bed", probesFile, "SharedPMDs", methylationMatrix)
  temp4=getMethylationHT450k("Data/MMSeekR_PMDs/ESCA_sharedHMDs.bed", probesFile, "SharedHMDs", methylationMatrix)
  plotdata=cbind(temp1, temp2, temp3, temp4)
  plotdata1=plotdata
  plotdata$SampleName=rownames(plotdata)
  plotdata=merge(plotdata, sampleInfo, by.x="SampleName", by.y="Sample")
  write.table(plotdata, gsub("pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
  my_theme=theme_classic()+ theme(plot.title = element_text(hjust = 0.5, size=14, color="black", face="bold"),
                                  axis.text.x = element_text(size=11, color="black", angle = 90, hjust = 1, vjust = 0.5),
                                  axis.text.y = element_text(size=11, color="black", angle = 0),
                                  axis.title = element_text(size=13, color="black", face="bold"),
                                  legend.title =element_text(size=12, color="black", face="bold"),
                                  legend.text =element_text(size=12, color="black"),legend.position = "none")
  plotList=list()
  for(type in unique(plotdata$Type)){
    data=plotdata[plotdata$Type%in%type,]
    data=melt(data, id.vars = c("SampleName", "Type"))
    colnames(data)=c("Sample", "Type", "GeneType", "Methylation")
    data$GeneType=factor(data$GeneType, levels=c("SharedPMDs", "EAC_specificPMDs", "ESCC_specificPMDs", "SharedHMDs"))
    p=ggplot(data, aes(x=GeneType, y=Methylation, group=Sample)) + geom_line(alpha=0.2, color="darkblue")+ geom_point(alpha=0.5, size=0.8, color="darkblue")
    p=p+xlab("")+ylab("methylation of each sample")+ggtitle(paste(type, " samples"))+ylim(0, 1)+my_theme
    plotList[[type]]=p
  }
  for(type in unique(plotdata$Type)){
    data=plotdata[plotdata$Type%in%type,]
    data2=data.frame(mean=colMeans(data[,2:5]), sd=colSds(as.matrix(data[,2:5])))
    data2$Type=rownames(data2)
    data2$Type=factor(data2$Type, levels=c("SharedPMDs", "EAC_specificPMDs", "ESCC_specificPMDs", "SharedHMDs"))
    if(length(grep("Normal", type))>0){
      ymin=0.3
      ymax=0.8
    }else{
      ymin=0
      ymax=1
    }
    p1<- ggplot(data2, aes(x=Type, y=mean, group=1)) + geom_point(shape=15) + geom_line() +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+ ylim(ymin, ymax)+my_theme
    p1=p1+theme(axis.ticks.x = element_blank())+xlab("")
    plotList[[paste0(type, "_trend")]]=p1
  }
  for(type in unique(plotdata$Type)){
    targetPlotdata=plotdata1[rownames(plotdata1)%in%plotdata[plotdata$Type%in%type,]$Sample,]
    targetPlotdata=targetPlotdata[order(rowSums(targetPlotdata)),]
    targetPlotdata=targetPlotdata[,c(3,2,1,4)]
    if(length(grep("Normal", type))>0){
      breaks=seq(0.3,0.8,by=0.001)
    }else{
      breaks=seq(0,1,by=0.001)
    }
    print(paste0(type, ":", nrow(targetPlotdata)))
    p1=pheatmap(targetPlotdata, show_rownames = F, cluster_rows = F, cluster_cols = F, breaks=breaks, border_color = NA, main=type,
                color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(length(breaks)))
    plotList[[paste0(type,"_heatmap")]]=p1[[4]]
  }
  pdf(saveFile, height = height, width=width)
  print(ggarrange(plotlist = plotList))
  dev.off()
}
plotHM450kMethylationLinePlot("meta/HT450k.probe.rmCGI.rmblackList.solo.bed", TCGAprobeMethylationMatrix, 
                              TCGAprobeSampleInfo, "Figure2D/Figure2D.TCGA.pdf", 13, 11)

################################################FigureS3D###############################################
load("Data/FigureS3D/EPIC.RData")
load("Data/FigureS3D/EPIC_probes.RData")
EPICprobeMethylationMatrix=result
EPICprobeSampleInfo=sampleInfo
EPICprobeSampleInfo[EPICprobeSampleInfo$Type%in%"EAC_Normal",]$Type="ESCC_Normal"
getMethylationEPIC=function(regionFile, probesInfo, type, methylationMatrix){
  tempData=read_bed(regionFile, n_fields = 5)
  tempData=bed_intersect(tempData, probesInfo)
  tempData=unique(tempData[tempData$.overlap==2,]$name.y)
  tempData=result[rownames(result)%in%tempData,]
  tempData=data.frame(colMeans(tempData, na.rm = T))
  colnames(tempData)=type
  return(tempData)
}
plotEPICMethylationLinePlot=function(probesInfo, methylationMatrix, sampleInfo, saveFile, height, width){
  temp1=getMethylationEPIC("Data/MMSeekR_PMDs/ESCC_specificPMDs.bed", probesInfo, "ESCC_specificPMDs", methylationMatrix)
  temp2=getMethylationEPIC("Data/MMSeekR_PMDs/EAC_specificPMDs.bed", probesInfo, "EAC_specificPMDs", methylationMatrix)
  temp3=getMethylationEPIC("Data/MMSeekR_PMDs/ESCA_sharedPMDs.bed", probesInfo, "SharedPMDs", methylationMatrix)
  temp4=getMethylationEPIC("Data/MMSeekR_PMDs/ESCA_sharedHMDs.bed", probesInfo, "SharedHMDs", methylationMatrix)
  plotdata=cbind(temp1, temp2, temp3, temp4)
  plotdata1=plotdata
  plotdata$SampleName=rownames(plotdata)
  plotdata=merge(plotdata, sampleInfo, by.x="SampleName", by.y="SampleName")
  write.table(plotdata, gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
  
  my_theme=theme_classic()+ theme(plot.title = element_text(hjust = 0.5, size=14, color="black", face="bold"),
                                  axis.text.x = element_text(size=11, color="black", angle = 90, hjust = 1, vjust = 0.5),
                                  axis.text.y = element_text(size=11, color="black", angle = 0),
                                  axis.title = element_text(size=13, color="black", face="bold"),
                                  legend.title =element_text(size=12, color="black", face="bold"),
                                  legend.text =element_text(size=12, color="black"),legend.position = "none")
  plotList=list()
  for(type in unique(plotdata$Type)){
    data=plotdata[plotdata$Type%in%type,]
    data=melt(data, id.vars = c("SampleName", "Type"))
    colnames(data)=c("Sample", "Type", "GeneType", "Methylation")
    data$GeneType=factor(data$GeneType, levels=c("SharedPMDs", "EAC_specificPMDs", "ESCC_specificPMDs", "SharedHMDs"))
    p=ggplot(data, aes(x=GeneType, y=Methylation, group=Sample)) + geom_line(alpha=0.2, color="darkblue")+ geom_point(alpha=0.5, size=0.8, color="darkblue")
    p=p+xlab("")+ylab("methylation of each sample")+ggtitle(type)+ylim(0, 1)+my_theme
    plotList[[type]]=p
  }
  
  for(type in unique(plotdata$Type)){
    data=plotdata[plotdata$Type%in%type,]
    data2=data.frame(mean=colMeans(data[,2:5]), sd=colSds(as.matrix(data[,2:5])))
    data2$Type=rownames(data2)
    data2$Type=factor(data2$Type, levels=c("SharedPMDs", "EAC_specificPMDs", "ESCC_specificPMDs", "SharedHMDs"))
    if(length(grep("Normal", type))>0){
      ymin=0.5
      ymax=0.8
    }else{
      ymin=0.2
      ymax=0.8
    }
    p1<- ggplot(data2, aes(x=Type, y=mean, group=1)) + geom_point(shape=15) + geom_line() +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2, position=position_dodge(.9))+ ylim(ymin, ymax)+my_theme+ggtitle(type)
    p1=p1+theme(axis.ticks.x = element_blank())+xlab("")
    plotList[[paste0(type, "_trend")]]=p1
  }
  
  for(type in unique(plotdata$Type)){
    targetPlotdata=plotdata1[rownames(plotdata1)%in%plotdata[plotdata$Type%in%type,]$Sample,]
    targetPlotdata=targetPlotdata[order(rowSums(targetPlotdata)),]
    targetPlotdata=targetPlotdata[,c(3,2,1,4)]
    if(length(grep("Normal", type))>0){
      breaks=seq(0.5,0.8,by=0.001)
    }else{
      breaks=seq(0.2,0.8,by=0.001)
    }
    print(paste0(type, ":", nrow(targetPlotdata)))
    p1=pheatmap(targetPlotdata, show_rownames = F, cluster_rows = F, cluster_cols = F, breaks=breaks, border_color = NA, main=type,
                color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(length(breaks)))
    plotList[[paste0(type,"_heatmap")]]=p1[[4]]
  }
  pdf(saveFile, height = height, width=width)
  print(ggarrange(plotlist = plotList))
  dev.off()
}
plotEPICMethylationLinePlot(EPIC_probes_rmblackList_rmCGI_solo, EPICprobeMethylationMatrix, 
                            EPICprobeSampleInfo, "FigureS3D/FigureS3D.EPIC.pdf", 9, 11)

################################################FigureS3E################################################
###Validate by WGBS
## data in Data/OtherESCC
#Shell FigureS3E/getRegionMeanBetaValuesForGSE149608.sh ESCA_sharedHMDs.bed
#Shell FigureS3E/getRegionMeanBetaValuesForGSE149608.sh ESCA_sharedPMDs.bed
#Shell FigureS3E/getRegionMeanBetaValuesForGSE149608.sh EAC_sharedPMDs_specificRegion.bed
#Shell FigureS3E/getRegionMeanBetaValuesForGSE149608.sh ESCC_sharedPMDs_specificRegion.bed
#merge the result "FigureS3E/GSE149608_regionLevels_methylation.soloCpGs.txt"
plotGeneLineplot=function(plotdata, type){
  plotdata$GeneType=factor(plotdata$GeneType, levels=c("SharedPMDs", "EAC_specificPMDs", "ESCC_specificPMDs", "SharedHMDs"))
  p=ggplot(data=plotdata, aes(x=GeneType, y=Methylation, group=Sample)) + geom_line(alpha=0.2, color="darkblue")+ geom_point(alpha=0.5, color="darkblue", size=0.8)+theme_classic()
  p=p+xlab("")+ylab("methylation of each sample")+ggtitle(type)+ylim(0, 1)
  p=p+theme(plot.title = element_text(hjust = 0.5, size=14, color="black", face="bold"),
            axis.text.x = element_text(size=11, color="black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size=11, color="black", angle = 0),
            axis.title = element_text(size=13, color="black", face="bold"),
            legend.title =element_text(size=12, color="black", face="bold"),
            legend.text =element_text(size=12, color="black"),legend.position = "none")
  return(p)
}
GSEWGBS=read.table("FigureS3E/GSE149608_regionLevels_methylation.soloCpGs.txt",header=T, stringsAsFactors =F, sep="\t")
GSEWGBS=melt(GSEWGBS, id.vars = c("Sample", "Type"))
colnames(GSEWGBS)=c("Sample", "SampleType", "GeneType", "Methylation")
p1=plotGeneLineplot(GSEWGBS[GSEWGBS$SampleType%in%"ESCC_Tumor",], "ESCC_Tumor samples")
p2=plotGeneLineplot(GSEWGBS[GSEWGBS$SampleType%in%"ESCC_Normal",], "ESCC_Normal samples")
pdf("FigureS3E/FigureS3E.pdf", width=7, height=4.5)
print(ggarrange(p1,p2, nrow=1))
dev.off()

################################################Figure2E and FigureS3F################################################
#####region compartmentB (predicted using TCGA 450k array)
domainFile1="Data/MMSeekR_PMDs/EAC_specificPMDs.bed"
domainFile2="Data/MMSeekR_PMDs/ESCC_specificPMDs.bed"
domainFile3="Data/MMSeekR_PMDs/ESCA_sharedPMDs.bed"
domainData1=read_bed(domainFile1)
domainData1=bed_merge(bed_sort(domainData1))
domainData2=read_bed(domainFile2)
domainData2=bed_merge(bed_sort(domainData2))
domainData3=read_bed(domainFile3)
domainData3=bed_merge(bed_sort(domainData3))

##origin data from /Volumes/YUAN_2T/WGBS_project/20200312_Unmask/PMD/CompartmentsA_B
comBFile1="Data/Figure2E_S3F/EACTumor_result.compartmentB.bed"
comBFile2="Data/Figure2E_S3F/ESCCTumor_result.compartmentB.bed"
comBData1=read_bed(comBFile1, n_fields = 4)
comBData1=comBData1[,1:3]
comBData1=bed_merge(bed_sort(comBData1))
comBData2=read_bed(comBFile2, n_fields = 4)
comBData2=comBData2[,1:3]
comBData2=bed_merge(bed_sort(comBData2))

domainData1_comB1=bed_intersect(domainData1, comBData1)
rate1=sum(domainData1_comB1$.overlap)/sum(domainData1$end-domainData1$start)
domainData1_comB2=bed_intersect(domainData1, comBData2)
rate2=sum(domainData1_comB2$.overlap)/sum(domainData1$end-domainData1$start)
domainData2_comB1=bed_intersect(domainData2, comBData1)
rate3=sum(domainData2_comB1$.overlap)/sum(domainData2$end-domainData2$start)
domainData2_comB2=bed_intersect(domainData2, comBData2)
rate4=sum(domainData2_comB2$.overlap)/sum(domainData2$end-domainData2$start)

domainData3_comB1=bed_intersect(domainData3, comBData1)
rate5=sum(domainData3_comB1$.overlap)/sum(domainData3$end-domainData3$start)
domainData3_comB2=bed_intersect(domainData3, comBData2)
rate6=sum(domainData3_comB2$.overlap)/sum(domainData3$end-domainData3$start)

plotdata=data.frame(type=c("EAC_tumor_PMDs", "EAC_tumor_PMDs", "ESCC_tumor_PMDs", "ESCC_tumor_PMDs"),
                    compartment=c("EAC_tumor_compB", "ESCC_tumor_compB", "EAC_tumor_compB", "ESCC_tumor_compB"),
                    rate=c(rate1, rate2, rate3, rate4))
plotdata$type=factor(plotdata$type, levels=c("EAC_tumor_PMDs", "ESCC_tumor_PMDs"))
plotdata$compartment=factor(plotdata$compartment, levels=c("EAC_tumor_compB", "ESCC_tumor_compB"))
write.table(plotdata, file="Figure2E/Figure2E_comparementB.txt")
p<-ggplot(data=plotdata, aes(x=type, y=rate, fill=compartment)) + geom_bar(stat="identity", position=position_dodge(), color="white")+theme_classic()
p=p+xlab("")+ylab("Fraction of PMDs overlapping with compartentB")
p=p+theme(axis.text =element_text(size=12), axis.title=element_text(size=14))+scale_fill_manual(values=c("#0000F5","#EA3323"))
pdf("Figure2E/Figure2E_comparementB.pdf")
print(p)
dev.off()

plotdata=data.frame(type=c("SharedPMDs", "SharedPMDs"), compartment=c("EAC_tumor_compB", "ESCC_tumor_compB"), rate=c(rate5, rate6))
plotdata$compartment=factor(plotdata$compartment, levels=c("EAC_tumor_compB", "ESCC_tumor_compB"))
write.table(plotdata, file="FigureS3F/FigureS3F_comparementB.txt")
p<-ggplot(data=plotdata, aes(x=type, y=rate, fill=compartment)) + geom_bar(stat="identity", position=position_dodge(), color="white")+theme_classic()
p=p+xlab("")+ylab("Fraction of PMDs overlapping with compartentB")
p=p+theme(axis.text =element_text(size=12), axis.title=element_text(size=14))+scale_fill_manual(values=c("#0000F5","#EA3323"))
pdf("FigureS3F/FigureS3F_comparementB.pdf")
print(p)
dev.off()


###Figure S3H and Figure S4I
##DNA methylation data from TCGA
load("meta/TCGA-ESCA_methylation_hg38.v2.RData")
library(TCGAbiolinks)
clinical <- GDCquery_clinic(project = "TCGA-ESCA", type = "clinical")
clinical_change=clinical[,c("submitter_id", "ajcc_pathologic_stage", "primary_diagnosis", "ajcc_pathologic_n",
                            "alcohol_history", "vital_status", "pack_years_smoked", "age_at_index")]
clinical_change$type="other"
clinical_change[clinical_change$primary_diagnosis%in%c("Adenocarcinoma, NOS"),]$type="EAC"
clinical_change[clinical_change$primary_diagnosis%in%c("Squamous cell carcinoma, NOS"),]$type="ESCC"
clinical_change=merge(clinical_change, match.file.cases[match.file.cases$Type%in%"Tumor",c(1,4)], by.x="submitter_id", by.y="Sample")
clinical_change=rbind(clinical_change[clinical_change$type%in%"EAC", ], clinical_change[clinical_change$type%in%"ESCC", ])

getMethylationHT450k=function(regionFile, probesFile, type, methylationMatrix){
  tempData=read_bed(regionFile,n_fields = 5)
  probesMatrix=read_bed(probesFile, n_fields = 4)
  tempData_probes=bed_intersect(probesMatrix,tempData)
  tempData_probes=as.data.frame(tempData_probes)
  tempData_probes=unique(tempData_probes[tempData_probes$.overlap==2,]$name.x)
  print(length(tempData_probes))
  tempData=methylationMatrix[rownames(methylationMatrix)%in%tempData_probes,, drop=F]
  tempData=data.frame(colMeans(tempData, na.rm = T))
  colnames(tempData)=type
  return(tempData)
}

ESCA_methy=result[ ,clinical_change$SampleName]
ESCA_methy_ESCCPMD=getMethylationHT450k("Data/MMSeekR_PMDs/ESCC_specificPMDs.bed", "meta/HT450k.probe.rmCGI.rmblackList.solo.bed", "ESCC",  ESCA_methy)
ESCA_methy_EACPMD=getMethylationHT450k("Data/MMSeekR_PMDs/EAC_specificPMDs.bed", "meta/HT450k.probe.rmCGI.rmblackList.solo.bed", "EAC",  ESCA_methy)
ESCA_methy_PMD=cbind(ESCA_methy_ESCCPMD, ESCA_methy_EACPMD)
ESCA_methy_PMD$delta=ESCA_methy_PMD$ESCC-ESCA_methy_PMD$EAC
ESCA_methy_PMD=data.frame(sample=rownames(ESCA_methy_PMD), ESCA_methy_PMD)
ESCA_methy_PMD=merge(ESCA_methy_PMD, clinical_change[,c(2,4,5,6,7,8,9,10)], by.x="sample", by.y="SampleName")
write.table(ESCA_methy_PMD, file="FigureS3H_S4I/FigureS3H.PMD.clinical.txt", row.names = F, col.names = T, sep="\t", quote=F)

ESCA_methy_ESCCDMR=getMethylationHT450k("Data/MaskUnionPMDs_DMRs/hypoESCC_Tumor.bed",
                                        "meta/HT450k.probe.rmblackList.bed", "ESCC",  ESCA_methy)
ESCA_methy_EACDMR=getMethylationHT450k("Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed",
                                       "meta/HT450k.probe.rmblackList.bed", "EAC",  ESCA_methy)
ESCA_methy_DMR=cbind(ESCA_methy_ESCCDMR, ESCA_methy_EACDMR)
ESCA_methy_DMR$delta=ESCA_methy_DMR$ESCC-ESCA_methy_DMR$EAC
ESCA_methy_DMR=data.frame(sample=rownames(ESCA_methy_DMR), ESCA_methy_DMR)
ESCA_methy_DMR=merge(ESCA_methy_DMR, clinical_change[,c(2,4,5,6,7,8,9,10)], by.x="sample", by.y="SampleName")
write.table(ESCA_methy_DMR, file="FigureS3H_S4I/FigureS4I.DMR.clinical.txt", row.names = F, col.names = T, sep="\t", quote=F)

plotClinicalFigures=function(type, ylab, ESCA_methy_data){
  ESCA_methy_data=ESCA_methy_data[ESCA_methy_data$type%in%type,]
  
  plotdata1=ESCA_methy_data[,c("sample", type, "age_at_index")]
  colnames(plotdata1)=c("sample", "methylation", "Age")
  plotdata1$group=">=60"
  plotdata1[plotdata1$Age<60,]$group="<60"
  plotdata1$group=factor(plotdata1$group, levels=c("<60", ">=60"))
  print(table(plotdata1$group))
  p1=ggplot(plotdata1, aes(x = group, y = methylation)) +geom_boxplot(width=0.5)+ geom_beeswarm(cex = 3)+theme_classic()
  p1=p1+xlab("")+ylab(ylab)+ylim(0,0.8)
  p1=p1+theme(axis.text = element_text(size=12, color="black"), axis.title.y = element_text(size=14))
  test=t.test(plotdata1[plotdata1$group%in%"<60",]$methylation, plotdata1[plotdata1$group%in%">=60",]$methylation)
  print(test$p.value)

  plotdata3=ESCA_methy_data[,c("sample", type, "ajcc_pathologic_n")]
  colnames(plotdata3)=c("sample", "methylation", "Lymph_node_metastasis")
  plotdata3$group=NA
  plotdata3[plotdata3$Lymph_node_metastasis%in%c("N0"),]$group="No"
  plotdata3[plotdata3$Lymph_node_metastasis%in%c("N1","N2","N3"),]$group="Yes"
  plotdata3=plotdata3[!is.na(plotdata3$group),]
  print(table(plotdata3$group))
  plotdata3$group=factor(plotdata3$group, levels=c("No", "Yes"))
  p3=ggplot(plotdata3, aes(x = group, y = methylation)) +geom_boxplot(width=0.5)+ geom_beeswarm(cex = 3)+theme_classic()
  p3=p3+xlab("Lymph node metastasis")+ylab(ylab)+ylim(0,0.8)
  p3=p3+theme(axis.text = element_text(size=12, color="black"), axis.title = element_text(size=14))
  test=t.test(plotdata3[plotdata3$group%in%"No",]$methylation, plotdata3[plotdata3$group%in%"Yes",]$methylation)
  print(test$p.value)
  
  plotdata4=ESCA_methy_data[,c("sample", type, "ajcc_pathologic_stage")]
  colnames(plotdata4)=c("sample", "methylation", "stage")
  plotdata4$group=NA
  plotdata4[plotdata4$stage%in%c("Stage IA", "Stage IB", "Stage IC", "Stage I",
                                 "Stage IIA", "Stage IIB", "Stage IIC", "Stage II"),]$group="Low stage (I-II)"
  plotdata4[plotdata4$stage%in%c("Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage III",
                                 "Stage IVA", "Stage IVB", "Stage IVC", "Stage IV"),]$group="High stage (III-IV)"
  plotdata4=plotdata4[!is.na(plotdata4$group),]
  plotdata4$group=factor(plotdata4$group, levels=c("Low stage (I-II)", "High stage (III-IV)"))
  print(table(plotdata4$group))
  p4=ggplot(plotdata4, aes(x = group, y = methylation)) +geom_boxplot(width=0.5)+ geom_beeswarm(cex = 3)+theme_classic()
  p4=p4+xlab("")+ylab(ylab)+ylim(0,0.8)
  p4=p4+theme(axis.text = element_text(size=12, color="black"), axis.title = element_text(size=14))
  test=t.test(plotdata4[plotdata4$group%in%"Low stage (I-II)",]$methylation, plotdata4[plotdata4$group%in%"High stage (III-IV)",]$methylation)
  print(test$p.value)
  
  plotdata5=ESCA_methy_data[,c("sample", type, "alcohol_history")]
  colnames(plotdata5)=c("sample", "methylation", "alcohol_history")
  plotdata5=plotdata5[!is.na(plotdata5$alcohol_history)&!plotdata5$alcohol_history%in%"Not Reported",]
  plotdata5$alcohol_history=factor(plotdata5$alcohol_history, levels=c("No", "Yes"))
  print(table(plotdata5$alcohol_history))
  p5=ggplot(plotdata5, aes(x = alcohol_history, y = methylation)) +geom_boxplot(width=0.5)+ geom_beeswarm(cex = 3)+theme_classic()
  p5=p5+xlab("alcohol history")+ylab(ylab)+ylim(0,0.8)
  p5=p5+theme(axis.text = element_text(size=12, color="black"), axis.title = element_text(size=14))
  test=t.test(plotdata5[plotdata5$alcohol_history%in%"No",]$methylation, plotdata5[plotdata5$alcohol_history%in%"Yes",]$methylation)
  print(test$p.value)
  
  plotdata6=ESCA_methy_data[,c("sample", type, "pack_years_smoked")]
  colnames(plotdata6)=c("sample", "methylation", "smoking_years")
  plotdata6$group=NA
  plotdata6[is.na(plotdata6$smoking_years),]$group="No"
  plotdata6[!is.na(plotdata6$smoking_years),]$group="Yes"
  plotdata6$group=factor(plotdata6$group, levels=c("No", "Yes"))
  print(table(plotdata6$group))
  p6=ggplot(plotdata6, aes(x = group, y = methylation)) +geom_boxplot(width=0.5)+ geom_beeswarm(cex = 3)+theme_classic()
  p6=p6+xlab("smoking status")+ylab(ylab)+ylim(0,0.8)
  p6=p6+theme(axis.text = element_text(size=12, color="black"), axis.title = element_text(size=14))
  test=t.test(plotdata6[plotdata6$group%in%"No",]$methylation, plotdata6[plotdata6$group%in%"Yes",]$methylation)
  print(test$p.value)
  
  plitList=list()
  plitList[["age"]]=p1
  # plitList[["vital"]]=p2
  plitList[["metastasis"]]=p3
  plitList[["stage"]]=p4
  plitList[["alcohol"]]=p5
  plitList[["smoking"]]=p6
  return(plitList)
}
ESCCDMR=plotClinicalFigures("ESCC", "ESCC specific DMRs methylation", ESCA_methy_DMR)
pdf("FigureS3H_S4I/FigureS4I.ESCC.DMR.clinical.pdf", width=15, height=3)
ggarrange(ESCCDMR$age, ESCCDMR$metastasis, ESCCDMR$stage, ESCCDMR$alcohol, ESCCDMR$smoking, nrow=1)
dev.off()
ESCCPMD=plotClinicalFigures("ESCC", "ESCC specific PMDs methylation", ESCA_methy_PMD)
pdf("FigureS3H_S4I/FigureS3H.ESCC.PMD.clinical.pdf", width=15, height=3)
ggarrange(ESCCPMD$age, ESCCPMD$metastasis, ESCCPMD$stage, ESCCPMD$alcohol, ESCCPMD$smoking, nrow=1)
dev.off()

EACDMR=plotClinicalFigures("EAC", "EAC specific DMRs methylation", ESCA_methy_DMR)
pdf("FigureS3H_S4I/FigureS4I.EAC.DMR.clinical.pdf", width=15, height=3)
ggarrange(EACDMR$age, EACDMR$metastasis, EACDMR$stage, EACDMR$alcohol, EACDMR$smoking, nrow=1)
dev.off()
EACPMD=plotClinicalFigures("EAC", "EAC specific PMDs methylation", ESCA_methy_PMD)
pdf("FigureS3H_S4I/FigureS3H.EAC.PMD.clinical.pdf", width=15, height=3)
ggarrange(EACPMD$age, EACPMD$metastasis, EACPMD$stage, EACPMD$alcohol, EACPMD$smoking, nrow=1)
dev.off()

###FigS3I and S4J###
##WGS data from 32194853
plotClinicalFigures=function(ylab, sampleFile, methylationFile, saveFile){
  ESCC_Pair42_methy=read_tsv(methylationFile)
  colnames(ESCC_Pair42_methy)=c("sample", "methylation")
  ESCC_Pair42_methy=as.data.frame(ESCC_Pair42_methy)
  data_info=read_tsv(sampleFile)
  data_info=as.data.frame(data_info)
  ESCC_Pair42_methy=merge(ESCC_Pair42_methy, data_info[,c(2:4,8:10)], by.x="sample", by.y="Name")
  
  
  plotdata1=ESCC_Pair42_methy[,c("sample", "methylation", "Age")]
  plotdata1$group=">=60"
  plotdata1[plotdata1$Age<60,]$group="<60"
  plotdata1$group=factor(plotdata1$group, levels=c("<60", ">=60"))
  p1=ggplot(plotdata1, aes(x = group, y = methylation)) +geom_boxplot(width=0.5)+ geom_beeswarm(cex = 4)+theme_classic()
  p1=p1+xlab("Age")+ylab(ylab)+ylim(0,0.85)
  p1=p1+theme(axis.text = element_text(size=12, color="black"), axis.title.y = element_text(size=14))
  test=t.test(plotdata1[plotdata1$group%in%"<60",]$methylation, plotdata1[plotdata1$group%in%">=60",]$methylation)
  print(test$p.value)
  write.table(plotdata1, gsub(".clinical.pdf",".age.txt", saveFile), row.names = F, col.names = T,sep="\t", quote=F)
  
  plotdata3=ESCC_Pair42_methy[,c("sample", "methylation", "Lymph_node_metastasis")]
  plotdata3$Lymph_node_metastasis=factor(plotdata3$Lymph_node_metastasis, levels=c("No", "Yes"))
  p3=ggplot(plotdata3, aes(x = Lymph_node_metastasis, y = methylation)) +geom_boxplot(width=0.5)+ geom_beeswarm(cex = 4)+theme_classic()
  p3=p3+xlab("Lymph node metastasis")+ylab(ylab)+ylim(0,0.85)
  p3=p3+theme(axis.text = element_text(size=12, color="black"), axis.title = element_text(size=14))
  test=t.test(plotdata3[plotdata3$Lymph_node_metastasis%in%"No",]$methylation, plotdata3[plotdata3$Lymph_node_metastasis%in%"Yes",]$methylation)
  print(test$p.value)
  write.table(plotdata3, gsub(".clinical.pdf",".metastasis.txt", saveFile), row.names = F, col.names = T,sep="\t", quote=F)
  
  plotdata4=ESCC_Pair42_methy[,c("sample", "methylation", "TNM_stage")]
  rstatix::anova_test(data=plotdata4, methylation ~ TNM_stage)
  print(table(ESCC_Pair42_methy$TNM_stage))
  p4=ggplot(plotdata4, aes(x = TNM_stage, y = methylation)) +geom_boxplot(width=0.5)+ geom_beeswarm(cex = 2)+theme_classic()
  p4=p4+ylab(ylab)+ylim(0,0.85)
  p4=p4+theme(axis.text = element_text(size=12, color="black"), axis.title = element_text(size=14))
  stageType=unique(plotdata4$TNM_stage)
  for(i in 1:(length(stageType)-1)){
    for(j in (i+1):length(stageType)){
      test=t.test(plotdata4[plotdata4$TNM_stage%in%stageType[i],]$methylation, plotdata4[plotdata4$TNM_stage%in%stageType[j],]$methylation)
      print(paste0(stageType[i], " ", stageType[j], " ", test$p.value))
    }
  }
  write.table(plotdata4, gsub(".clinical.pdf",".TNM_stage.txt", saveFile), row.names = F, col.names = T,sep="\t", quote=F)
  
  plitList=list()
  plitList[["age"]]=p1
  plitList[["metastasis"]]=p3
  plitList[["TNM_stage"]]=p4
  pdf(saveFile, width=10, height=3)
  print(ggarrange(ESCCPMD$age, ESCCPMD$metastasis, ESCCPMD$TNM_stage, nrow=1))
  dev.off()
}
plotClinicalFigures("ESCC-specific\nPMDs methylation", "Data/FigureS3I_S4J/Pair42_sample_info.txt",
                    "Data/FigureS3I_S4J/ESCC_PMD.result.txt", "FigureS3I_S4J/FigureS3I.PMD.clinical.pdf")

plotClinicalFigures("ESCC-specific\nDMRs methylation", "Data/FigureS3I_S4J/Pair42_sample_info.txt",
                    "Data/FigureS3I_S4J/ESCC_DMR.result.txt", "FigureS3I_S4J/FigureS4J.DMR.clinical.pdf")

################################################Figure2F################################################
########Figure2F EBI EAC mutation
plotMutationBoxplot=function(mutationData, saveFile){
  p <- ggplot(mutationData, aes(x=Type, y=Mutation)) +  geom_boxplot(width=0.8)+ylim(0,20)+theme_classic()+xlab("")
  p=p+theme(axis.title=element_text(size=12, color="black"), axis.text = element_text(size=10, color="black"), axis.line = element_line(size=0.5))
  pdf(saveFile,height=6)
  print(p)
  dev.off()
  
  mutation_pvalue=as.data.frame(matrix(numeric(0), ncol=3))
  sharedPMDs_mutant=mutationData[mutationData$Type%in%"sharedPMDs",]$Mutation
  EAC_specificPMDs_mutant=mutationData[mutationData$Type%in%"EAC_specificPMDs",]$Mutation
  ESCC_specificPMDs_mutant=mutationData[mutationData$Type%in%"ESCC_specificPMDs",]$Mutation
  sharedHMDs_mutant=mutationData[mutationData$Type%in%"sharedHMDs",]$Mutation
  mutation_pvalue=rbind(mutation_pvalue, data.frame(Type1="sharedPMDs", Type2="EAC_specificPMDs", pvalue=t.test(sharedPMDs_mutant, EAC_specificPMDs_mutant)$p.value))
  mutation_pvalue=rbind(mutation_pvalue, data.frame(Type1="sharedPMDs", Type2="ESCC_specificPMDs", pvalue=t.test(sharedPMDs_mutant, ESCC_specificPMDs_mutant)$p.value))
  mutation_pvalue=rbind(mutation_pvalue, data.frame(Type1="sharedPMDs", Type2="sharedHMDs", pvalue=t.test(sharedPMDs_mutant, sharedHMDs_mutant)$p.value))
  mutation_pvalue=rbind(mutation_pvalue, data.frame(Type1="EAC_specificPMDs", Type2="ESCC_specificPMDs", pvalue=t.test(EAC_specificPMDs_mutant, ESCC_specificPMDs_mutant)$p.value))
  mutation_pvalue=rbind(mutation_pvalue, data.frame(Type1="EAC_specificPMDs", Type2="sharedHMDs", pvalue=t.test(EAC_specificPMDs_mutant, sharedHMDs_mutant)$p.value))
  mutation_pvalue=rbind(mutation_pvalue, data.frame(Type1="ESCC_specificPMDs", Type2="sharedHMDs", pvalue=t.test(ESCC_specificPMDs_mutant, sharedHMDs_mutant)$p.value))
  write.table(mutation_pvalue, gsub(".pdf", "_pvalue.txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
}
EAC_specificPMDs=read_bed("Data/MMSeekR_PMDs/EAC_specificPMDs.hg19.bed")
ESCC_specificPMDs=read_bed("Data/MMSeekR_PMDs/ESCC_specificPMDs.hg19.bed")
ESCA_sharePMDs=read_bed("Data/MMSeekR_PMDs/ESCA_sharedPMDs.hg19.bed")
ESCA_shareHMDs=read_bed("Data/MMSeekR_PMDs/ESCA_sharedHMDs.hg19.bed")
vcfFileList=dir("../EAC_EBI_vcf/", pattern = "vcf", full=T)
mutationResult=as.data.frame(matrix(numeric(0), ncol=5))
getDomainMutantRatio=function(mutantData, domainData){
  mutant_domainData=bed_intersect(mutantData, domainData)
  mutant_domainData=mutant_domainData[,1:3]
  ratio=nrow(mutant_domainData)/sum(domainData$end-domainData$start)*(1e6)
}
for(vcfFile in vcfFileList){
  fileIndex=strsplit(basename(vcfFile), "__")[[1]][1]
  print(fileIndex)
  data=read_tsv(vcfFile, comment = "##",col_types=cols(
    `#CHROM` = col_character(),
    POS = col_integer(),
    ID = col_character(),
    REF = col_character(),
    ALT = col_character(),
    QUAL = col_character(),
    FILTER = col_character(),
    INFO = col_character(),
    FORMAT = col_character(),
    NORMAL = col_character(),
    TUMOR = col_character()
  ))
  data=unique(tibble(chrom=data$`#CHROM`, start=data$POS-1, end=data$POS))
  data_EACspecificPMDs=getDomainMutantRatio(data, EAC_specificPMDs)
  data_ESCCspecificPMDs=getDomainMutantRatio(data, ESCC_specificPMDs)
  data_sharedPMDs=getDomainMutantRatio(data, ESCA_sharePMDs)
  data_sharedHMDs=getDomainMutantRatio(data, ESCA_shareHMDs)
  
  temp=data.frame(sample=fileIndex, EAC_specificPMDs=data_EACspecificPMDs, ESCC_specificPMDs=data_ESCCspecificPMDs,
                  sharedPMDs=data_sharedPMDs, sharedHMDs= data_sharedHMDs)
  mutationResult=rbind(mutationResult, temp)
}
mutationResult=melt(mutationResult, id.vars = "sample")
colnames(mutationResult)=c("Sample", "Type", "Mutation")
mutationResult$Type=factor(as.character(mutationResult$Type), levels=c("sharedPMDs", "EAC_specificPMDs", "ESCC_specificPMDs", "sharedHMDs"))
write.table(mutationResult, "Figure2F/Figure2F_EAC.txt", row.names = F, col.names = T, sep="\t", quote=F)
plotMutationBoxplot(mutationResult, "Figure2F/Figure2F_EAC.pdf")

####Cell research: PMID32398863 ESCC mutation
mutationData=read.csv("Data/Figure2F/PMID32398863_41422_2020_333_MOESM23_ESM.csv")
mutationResult=as.data.frame(matrix(numeric(0), ncol=5))
for(tumorSample in unique(as.character(mutationData$Tumor_Sample_Barcode))){
  tumorMutant=mutationData[mutationData$Tumor_Sample_Barcode%in%tumorSample,]
  tumorMutant=tumorMutant[tumorMutant$Variant_Type%in%"SNP",]
  tumorMutant=tibble(chrom=as.character(tumorMutant$Chromosome), start=tumorMutant$Start_Position-1, end=tumorMutant$End_Position)
  getDomainMutantRatio=function(mutantData, domainData){
    mutant_domainData=bed_intersect(mutantData, domainData)
    mutant_domainData=mutant_domainData[,1:3]
    ratio=nrow(mutant_domainData)/sum(domainData$end-domainData$start)*(1e6)
  }
  data_EACspecificPMDs=getDomainMutantRatio(tumorMutant, EAC_specificPMDs)
  data_ESCCspecificPMDs=getDomainMutantRatio(tumorMutant, ESCC_specificPMDs)
  data_sharedPMDs=getDomainMutantRatio(tumorMutant, ESCA_sharePMDs)
  data_sharedHMDs=getDomainMutantRatio(tumorMutant, ESCA_shareHMDs)
  temp=data.frame(sample=tumorSample, EAC_specificPMDs=data_EACspecificPMDs, ESCC_specificPMDs=data_ESCCspecificPMDs,
                  sharedPMDs=data_sharedPMDs, sharedHMDs= data_sharedHMDs)
  mutationResult=rbind(mutationResult, temp)
}
mutationResult=melt(mutationResult, id.vars = "sample")
colnames(mutationResult)=c("Sample", "Type", "Mutation")
mutationResult$Type=factor(as.character(mutationResult$Type), levels=c("sharedPMDs", "EAC_specificPMDs", "ESCC_specificPMDs", "sharedHMDs"))
write.table(mutationResult, "Figure2F/Figure2F_ESCC.txt", row.names = F, col.names = T, sep="\t", quote=F)
plotMutationBoxplot(mutationResult, "Figure2F/Figure2F_ESCC.pdf")

##########Figure S3G signitures
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
genome <- BSgenome::getBSgenome(ref_genome)
EAC_specificPMDs=read_bed("Data/MMSeekR_PMDs/EAC_specificPMDs.hg19.bed")
ESCC_specificPMDs=read_bed("Data/MMSeekR_PMDs/ESCC_specificPMDs.hg19.bed")
ESCA_sharePMDs=read_bed("Data/MMSeekR_PMDs/ESCA_sharedPMDs.hg19.bed")
ESCA_shareHMDs=read_bed("Data/MMSeekR_PMDs/ESCA_sharedHMDs.hg19.bed")

vcfFileList=dir("../EAC_EBI_vcf/", pattern = "vcf", full=T)
vcfFileSample=gsub("__pv.1.[1-9]__rg.grch37_g1k__al.bwa_mem__.snp.pass.vcf", "",basename(vcfFileList))
vcfFileData=data.frame(sample=vcfFileSample, vcf=vcfFileList)

grl_EACspecificPMDList=list()
gral_allList=list()
for(i in 1:nrow(vcfFileData)){
  print(i)
  grl <- purrr::map(vcfFileData[i,2], MutationalPatterns:::.read_single_vcf_as_grange, 
                    genome, group="auto", change_seqnames=T) %>% GenomicRanges::GRangesList()
  grl=get_mut_type(grl)
  ###all SNV
  gral_allList[[vcfFileData[i,1]]]=grl[[1]]
  
  grl_data=gr_to_bed(grl[[1]])
  grl_data$name=names(grl[[1]])
  
  ###EAC_specificPMD_SNVs
  grl_data2=bed_intersect(grl_data, EAC_specificPMDs)
  grl_EACspecificPMD=grl[[1]][names(grl[[1]])%in%grl_data2$name.x,]
  grl_EACspecificPMDList[[vcfFileData[i,1]]]=grl_EACspecificPMD
}
grl_EACspecificPMD_mat <- mut_matrix(vcf_list = grl_EACspecificPMDList, ref_genome = ref_genome)
grl_all_mat <- mut_matrix(vcf_list = gral_allList, ref_genome = ref_genome)
save(grl_EACspecificPMD_mat, grl_all_mat, file="FigureS3G/EAC_mutation96.RData")
grl_all_mat2=data.frame("Mutation Types"=rownames(grl_all_mat), grl_all_mat, check.names = F)
write.table(grl_all_mat2, file="FigureS3G/EAC_all_mutation.txt", row.names = F, col.names = T, sep="\t", quote=F)
grl_EACspecificPMD_mat2=data.frame("Mutation Types"=rownames(grl_EACspecificPMD_mat), grl_EACspecificPMD_mat, check.names = F)
write.table(grl_EACspecificPMD_mat2, file="FigureS3G/EAC_specificPMD_mutation.txt", row.names = F, col.names = T, sep="\t", quote=F)

####ESCC mutation
mutationData=read.csv("Data/Figure2F/PMID32398863_41422_2020_333_MOESM23_ESM.csv")
mutationData=mutationData[mutationData$Variant_Type%in%"SNP",]
mutationData$lable=paste0(mutationData$Tumor_Sample_Barcode,":", mutationData$Matched_Norm_Sample_Barcode)
mutationData2=split(mutationData, mutationData$lable)
grl_ESCCspecificPMDList=list()
gral_allList=list()
for(i in 1:length(mutationData2)){
  print(i)
  grlData=mutationData2[[i]]
  grlData=grlData[,c(2:4, 6, 7, 9)]
  colnames(grlData)[4:6]=c("TYPE", "REF", "ALT")
  grl= makeGRangesFromDataFrame(grlData,
                                keep.extra.columns=T,
                                ignore.strand=T,
                                seqinfo=Seqinfo(genome="hg19"),
                                seqnames.field="Chromosome",
                                start.field="Start_Position",
                                end.field="End_Position",
                                strand.field="strand",
                                starts.in.df.are.0based=FALSE)
  names(grl)=paste0(grlData$Chromosome, ":", grlData$Start_Position, "_", grlData$REF, "/", grlData$ALT)
  grl$QUAL=NA
  grl$FILTER="PASS"
  grl=get_mut_type(grl)
  ###all SNV
  gral_allList[[names(mutationData2)[i]]]=grl
  
  grlData=gr_to_bed(grl)
  grlData$name=names(grl)
  
  ###ESCC_specificPMD_SNVs
  grl_data2=bed_intersect(grlData, ESCC_specificPMDs)
  grl_ESCCspecificPMD=grl[names(grl)%in%grl_data2$name.x,]
  grl_ESCCspecificPMDList[[names(mutationData2)[i]]]=grl_ESCCspecificPMD
}
grl_ESCCspecificPMD_mat <- mut_matrix(vcf_list = grl_ESCCspecificPMDList, ref_genome = ref_genome)
grl_ESCAsharePMDs_mat <- mut_matrix(vcf_list = grl_ESCAsharePMDsList, ref_genome = ref_genome)
grl_all_mat <- mut_matrix(vcf_list = gral_allList, ref_genome = ref_genome)
save(grl_ESCCspecificPMD_mat, grl_all_mat, file="FigureS3G/ESCC_mutation96.RData")
grl_all_mat2=data.frame("Mutation Types"=rownames(grl_all_mat), grl_all_mat, check.names = F)
write.table(grl_all_mat2, file="FigureS3G/ESCC_all_mutation.txt", row.names = F, col.names = T, sep="\t", quote=F)
grl_ESCCspecificPMD_mat2=data.frame("Mutation Types"=rownames(grl_ESCCspecificPMD_mat), grl_ESCCspecificPMD_mat, check.names = F)
write.table(grl_ESCCspecificPMD_mat2, file="FigureS3G/ESCC_specificPMD_mutation.txt", row.names = F, col.names = T, sep="\t", quote=F)


getSigitureSNV=function(mut_mat, signatures, max_delta){
  mut_mat=mut_mat+0.001
  fit_res <- fit_to_signatures_strict(mut_mat, signatures, max_delta = max_delta)
  myPlotdata=fit_res$fit_res$contribution
  myPlotdata=myPlotdata[rowSums(myPlotdata>0)>=4,]
  result=as.data.frame(fit_res$fit_res$contribution)
  return(result)
}
signatures = get_known_signatures()
signatures = signatures[,!colnames(signatures)%in%
                          c("SBS85", "SBS84","SBS85","SBS86","SBS87","SBS88","SBS89","SBS90")]
load("FigureS3G/EAC_mutation96.RData")
EAC_all_mutation_0.05=getSigitureSNV(grl_all_mat, signatures, 0.05)
EAC_specific_mutation_0.05=getSigitureSNV(grl_EACspecificPMD_mat, signatures, 0.05)
load("FigureS3G/ESCC_mutation96.RData")
ESCC_all_mutation_0.05=getSigitureSNV(grl_all_mat, signatures, 0.05)
ESCC_specific_mutation_0.05=getSigitureSNV(grl_ESCCspecificPMD_mat, signatures, 0.05)

combindResult=function(dataset1, dataset2, dataset3, dataset4, 
                       dataset1.length, dataset2.length, dataset3.length, dataset4.length,
                       dataset1.name, dataset2.name, dataset3.name, dataset4.name, saveFile){
  dataset1_sig=data.frame(dataset1=rowSums(dataset1>0)/ncol(dataset1))
  dataset2_sig=data.frame(dataset2=rowSums(dataset2>0)/ncol(dataset2))
  dataset3_sig=data.frame(dataset3=rowSums(dataset3>0)/ncol(dataset3))
  dataset4_sig=data.frame(dataset4=rowSums(dataset4>0)/ncol(dataset4))
  mergeSig=cbind(dataset1_sig, dataset2_sig, dataset3_sig, dataset4_sig)
  colnames(mergeSig)=c(dataset1.name, dataset2.name, dataset3.name, dataset4.name)
  mergeSig=mergeSig[rowSums(mergeSig<0.1)!=4,]
  
  dataset1_TMB= melt(t(dataset1))
  dataset1_TMB=dataset1_TMB[dataset1_TMB$value>0,]
  dataset1_TMB=data.frame(dataset1_TMB%>%group_by(Var2)%>%summarise(median(value)))
  colnames(dataset1_TMB)=c("sig", "dataset1")
  dataset1_TMB$dataset1=dataset1_TMB$dataset1/dataset1.length*1e6
  
  dataset2_TMB= melt(t(dataset2))
  dataset2_TMB=dataset2_TMB[dataset2_TMB$value>0,]
  dataset2_TMB=data.frame(dataset2_TMB%>%group_by(Var2)%>%summarise(median(value)))
  colnames(dataset2_TMB)=c("sig", "dataset2")
  dataset2_TMB$dataset2=dataset2_TMB$dataset2/dataset2.length*1e6
  
  dataset3_TMB= melt(t(dataset3))
  dataset3_TMB=dataset3_TMB[dataset3_TMB$value>0,]
  dataset3_TMB=data.frame(dataset3_TMB%>%group_by(Var2)%>%summarise(median(value)))
  colnames(dataset3_TMB)=c("sig", "dataset3")
  dataset3_TMB$dataset3=dataset3_TMB$dataset3/dataset3.length*1e6
  
  dataset4_TMB= melt(t(dataset4))
  dataset4_TMB=dataset4_TMB[dataset4_TMB$value>0,]
  dataset4_TMB=data.frame(dataset4_TMB%>%group_by(Var2)%>%summarise(median(value)))
  colnames(dataset4_TMB)=c("sig", "dataset4")
  dataset4_TMB$dataset4=dataset4_TMB$dataset4/dataset4.length*1e6
  
  mergeTMB=merge(dataset1_TMB, dataset2_TMB, all=T)
  mergeTMB=merge(mergeTMB, dataset3_TMB, all=T)
  mergeTMB=merge(mergeTMB, dataset4_TMB, all=T)
  rownames(mergeTMB)=mergeTMB$sig
  mergeTMB=mergeTMB[rownames(mergeSig),]
  colnames(mergeTMB)[-1]=c(dataset1.name, dataset2.name, dataset3.name, dataset4.name)
  mergeTMB=melt(mergeTMB, id.vars = "sig")
  colnames(mergeTMB)=c("sig", "type", "TMB")
  
  mergeSig=data.frame(sig=rownames(mergeSig), mergeSig)
  mergeSig=melt(mergeSig, id.vars = "sig")
  colnames(mergeSig)=c("sig", "type", "proportion")
  mergeSig[mergeSig$proportion==0,]$proportion=NA
  
  plotdata=merge(mergeTMB, mergeSig)
  plotdata$TMB_type=NA
  tmpData=quantile(plotdata$TMB, probs = seq(0, 1, 0.1), na.rm = T)
  tmpData=round(tmpData, 2)
  for(i in 2:length(tmpData)){
    print(i)
    plotdata[!is.na(plotdata$TMB)&plotdata$TMB>tmpData[[i-1]]&plotdata$TMB<=tmpData[[i]],]$TMB_type=as.character(tmpData[[i]])
  }
  plotdata$TMB_type=factor(plotdata$TMB_type, levels=as.character(tmpData))
  plotdata$sig=factor(plotdata$sig, levels=rev(rownames(dataset1)))
  plotdata$type=factor(plotdata$type, levels=c(dataset1.name, dataset2.name, dataset3.name, dataset4.name))
  
  p=ggplot(plotdata, aes(x=type, y=sig, size = proportion, color=TMB_type)) +
    geom_point()+scale_size(breaks = seq(0,1,0.1), range = c(1, 12)) +
    scale_color_manual(name="TMB", values=c("#FFC74E","#FFAA5B","#FF8B68",
                                            "#FF6A75","#DC6FA5","#AA74D1","#7478FD","#4C68D8","#1A59B1", "#004889"))
  #    scale_color_manual(name="TMB", values=colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(10))
  #    scale_color_manual(name="TMB", values=colorRampPalette(c("green","blue","red"))(10))
  
  
  p=p+theme_classic()+xlab("")+ylab("")
  p=p+theme(axis.text.x = element_text(size=12, angle = 30, hjust = 1, color="black"),
            axis.text.y = element_text(size=12, color="black"))
  pdf(saveFile, width=5, height=6.5)
  print(p)
  dev.off()
  write.table(plotdata, file=gsub(".pdf", "_data.txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
}
combindResult(EAC_all_mutation_0.05, EAC_specific_mutation_0.05, ESCC_all_mutation_0.05, ESCC_specific_mutation_0.05,
              3*1e9, sum(EAC_specificPMDs$end-EAC_specificPMDs$start), 3*1e9, sum(ESCC_specificPMDs$end-ESCC_specificPMDs$start),
              "EAC_all", "EAC_specificPMD", "ESCC_all", "ESCC_specificPMD", 
              "FigureS3G/combine.sig.0.05.pdf")

################################################Figure3AB################################################
load("meta/TCGA-ESCA_Expression_hg38.v2.RData")
load("meta/TCGA-ESCA.clinc.RData")
EACTumor_expression=expression[,colnames(expression)%in%sampleInformation[sampleInformation$Label%in%"EAC_Tumor",]$SampleName]
ESCCTumor_expression=expression[,colnames(expression)%in%sampleInformation[sampleInformation$Label%in%"ESCC_Tumor",]$SampleName]
gene_region=read_bed("meta/gencode.v31.basic.gene.bed", n_fields = 4)
gene_region$length=gene_region$end-gene_region$start
getRegionLevelExp=function(domainFile, name){
  domainData=read_bed(domainFile)
  domainData$length=domainData$end-domainData$start
  domain_relatedGenes=bed_intersect(gene_region, domainData)
  domain_relatedGenes=domain_relatedGenes[domain_relatedGenes$.overlap>0.5*domain_relatedGenes$length.x,]
  domain_relatedGenes=unique(domain_relatedGenes$name.x)
  
  relatedGenes_EACTumor=EACTumor_expression[rownames(EACTumor_expression)%in%domain_relatedGenes,]
  relatedGenes_EACTumor=log2(relatedGenes_EACTumor+0.1)
  relatedGenes_EACTumor=data.frame(Gene=rowMeans(relatedGenes_EACTumor))
  relatedGenes_EACTumor$Type="EAC_Tumor"
  
  relatedGenes_ESCCTumor=ESCCTumor_expression[rownames(ESCCTumor_expression)%in%domain_relatedGenes,]
  relatedGenes_ESCCTumor=log2(relatedGenes_ESCCTumor+0.1)
  relatedGenes_ESCCTumor=data.frame(Gene=rowMeans(relatedGenes_ESCCTumor))
  relatedGenes_ESCCTumor$Type="ESCC_Tumor"
  
  relatedGenes_exp=rbind(relatedGenes_EACTumor, relatedGenes_ESCCTumor)
  colnames(relatedGenes_exp)=c("Expression", "Type")
  relatedGenes_exp$Domain=name
  return(relatedGenes_exp)
}

EAC_PMDs_exp=getRegionLevelExp("Data/MMSeekR_PMDs/EAC_specificPMDs.bed", "EAC_specificPMDs")
ESCC_PMDs_exp=getRegionLevelExp("Data/MMSeekR_PMDs/ESCC_specificPMDs.bed", "ESCC_specificPMDs")
ESCA_sharePMDs_exp=getRegionLevelExp("Data/MMSeekR_PMDs/ESCA_sharedPMDs.bed", "sharedPMDs")
ESCA_shareHMDs_exp=getRegionLevelExp("Data/MMSeekR_PMDs/ESCA_sharedHMDs.bed", "sharedHMDs")
domainExpResult=rbind(EAC_PMDs_exp, ESCC_PMDs_exp)
domainExpResult=rbind(domainExpResult, ESCA_sharePMDs_exp)
domainExpResult=rbind(domainExpResult, ESCA_shareHMDs_exp)
domainExpResult$Domain=factor(domainExpResult$Domain, levels=c("sharedPMDs", "EAC_specificPMDs", "ESCC_specificPMDs", "sharedHMDs"))
EACTumor_domainExpResult=domainExpResult[domainExpResult$Type%in%"EAC_Tumor",]
ESCCTumor_domainExpResult=domainExpResult[domainExpResult$Type%in%"ESCC_Tumor",]

plotBoxplot=function(expressionData, saveFile){
  p <- ggplot(expressionData, aes(x=Domain, y=Expression)) +  geom_boxplot(width=0.8)+ylim(-4,4)+theme_classic()+xlab("")
  p=p+theme(axis.title=element_text(size=12, color="black"), axis.text = element_text(size=10, color="black"), axis.line = element_line(size=0.5))
  pdf(saveFile,height=6)
  print(p)
  dev.off()
  expressionData=data.frame(gene=rownames(expressionData), expressionData)
  write.table(expressionData, gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
  expression_pvalue=as.data.frame(matrix(numeric(0), ncol=3))
  A=expressionData[expressionData$Domain%in%"sharedPMDs",]$Expression
  B=expressionData[expressionData$Domain%in%"EAC_specificPMDs",]$Expression
  C=expressionData[expressionData$Domain%in%"ESCC_specificPMDs",]$Expression
  D=expressionData[expressionData$Domain%in%"sharedHMDs",]$Expression
  expression_pvalue=rbind(expression_pvalue, data.frame(Type1="sharedPMDs", Type2="EAC_specificPMDs", pvalue=t.test(A, B)$p.value))
  expression_pvalue=rbind(expression_pvalue, data.frame(Type1="sharedPMDs", Type2="ESCC_specificPMDs", pvalue=t.test(A, C)$p.value))
  expression_pvalue=rbind(expression_pvalue, data.frame(Type1="sharedPMDs", Type2="sharedHMDs", pvalue=t.test(A, D)$p.value))
  expression_pvalue=rbind(expression_pvalue, data.frame(Type1="EAC_specificPMDs", Type2="ESCC_specificPMDs", pvalue=t.test(B, C)$p.value))
  expression_pvalue=rbind(expression_pvalue, data.frame(Type1="EAC_specificPMDs", Type2="sharedHMDs", pvalue=t.test(B, D)$p.value))
  expression_pvalue=rbind(expression_pvalue, data.frame(Type1="ESCC_specificPMDs", Type2="sharedHMDs", pvalue=t.test(C, D)$p.value))
  write.table(expression_pvalue, gsub(".pdf", "_pvalue.txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
}
plotBoxplot(EACTumor_domainExpResult, "Figure3AB/Figure3A_EAC.pdf")
plotBoxplot(ESCCTumor_domainExpResult, "Figure3AB/Figure3B_ESCC.pdf")


################################################Figure3CD################################################
load("meta/TCGA_ESCCvsEAC.DESeq2.RData")
ESCC_EAC_result=ESCC_EAC_result[rowMeans(ESCC_EAC_result[,-1:-8])>0.1,]
ESCC_EAC_result=ESCC_EAC_result[!is.na(ESCC_EAC_result$padj),]
ESCC_EAC_result=ESCC_EAC_result[,1:8]
ESCC_EAC_result=data.frame(GeneID=rownames(ESCC_EAC_result), ESCC_EAC_result)
write.table(ESCC_EAC_result, file="TCGA_ESCCvsEAC.DESeq2.txt", row.names = F, col.names = T, sep="\t", quote=F)

load("meta/TCGA_EACvsNormal.DESeq2.RData")
EAC_TN_result=EAC_TN_result[rowMeans(EAC_TN_result[,-1:-8])>0.1,]
EAC_TN_result=EAC_TN_result[!is.na(EAC_TN_result$padj),]
EAC_TN_result=EAC_TN_result[,1:8]
EAC_TN_result=data.frame(GeneID=rownames(EAC_TN_result), EAC_TN_result)
write.table(EAC_TN_result, file="TCGA_EACvsNormal.DESeq2.txt", row.names = F, col.names = T, sep="\t", quote=F)

load("meta/ESCC_Paired_result.RData")
ESCC_Paired_resuts=ESCC_Paired_resuts[rowMeans(ESCC_Paired_resuts[,3:22])>0.1,]
ESCC_Paired_resuts=ESCC_Paired_resuts[!is.na(ESCC_Paired_resuts$padj),]
ESCC_Paired_resuts=ESCC_Paired_resuts[,c(1,25:30,24,23)]
write.table(ESCC_Paired_resuts, file="ESCCTumorvsNorml.DESeq2.txt", row.names = F, col.names = T, sep="\t", quote=F)

##using cistrom-GO website
##EAC specific PMDs ("Data/MMSeekR_PMDs/EAC_specificPMDs.bed" and EAC down-regulated genes in "meta/TCGA_ESCCvsEAC.DESeq2.txt")
##ESCC specific PMDs ("Data/MMSeekR_PMDs/ESCC_specificPMDs.bed" and ESCC down-regulated genes in "meta/TCGA_ESCCvsEAC.DESeq2.txt")

plotCistomePathway=function(cistromePathwayFile){
  cistromePathwayResult=read_tsv(cistromePathwayFile, col_names=T)
  colnames(cistromePathwayResult)=c("GO_term", "detail", "enrichment", "pValue", "FDR", "Genes")
  cistromePathwayResult=cistromePathwayResult[order(cistromePathwayResult$FDR),]
  cistromePathwayResult=cistromePathwayResult[cistromePathwayResult$FDR<0.05,]
  cistromePathwayResult$FDR=-log10(cistromePathwayResult$FDR)
  if(nrow(cistromePathwayResult)>15){
    cistromePathwayResult=cistromePathwayResult[1:15,]
  }
  cistromePathwayResult$GO_term=gsub("\\(.*", "", cistromePathwayResult$GO_term)
  cistromePathwayResult$GO_term=factor(cistromePathwayResult$GO_term, levels=rev(unique(cistromePathwayResult$GO_term)))
  
  p<-ggplot(data=cistromePathwayResult, aes(x=GO_term, y=FDR)) + geom_bar(stat="identity", fill="black")+ coord_flip()+theme_classic()
  p=p+ylab("-log10(FDR)")+xlab("")
  p=p+theme(axis.text = element_text(colour = "black"))
  pdf(gsub(".txt", ".pdf", cistromePathwayFile),width=10, height=2+nrow(cistromePathwayResult)/5)
  print(p)
  dev.off()
}
plotCistomePathway("Figure3CD/Figure3C_EAC_specificPMDs_GO_BP.txt")
plotCistomePathway("Figure3CD/Figure3D_ESCC_specificPMDs_GO_BP.txt")

################################################Figure3E################################################
load("Data/Figure2A/hg38_10k.rmblackListCGI.RData")
hg38_10k_data=hg38_10k[,c(2:4,1)]
hg38_10k_data=tibble(chrom=as.character(hg38_10k_data$chr), start=hg38_10k_data$start, end=hg38_10k_data$end, region=hg38_10k_data$region)
targetRegion=read_bed("Data/Figure2A/KRT74_SLC2A2.bed", n_fields = 4)
targetRegion=bed_intersect(targetRegion, hg38_10k_data)
targetRegion=targetRegion[targetRegion$.overlap==10000,]

KRT74_regions=hg38_10k[hg38_10k$region%in%targetRegion[targetRegion$name.x%in%"KRT74",]$region.y,]
write.table(KRT74_regions, file="Figure3E/Figure3E_KRT74_methylation.txt", row.names = F, col.names = T, sep="\t", quote=F)
SLC2A2_regions=hg38_10k[hg38_10k$region%in%targetRegion[targetRegion$name.x%in%"SLC2A2",]$region.y,]
write.table(SLC2A2_regions, file="Figure3E/Figure3E_SLC2A2_methylation.txt", row.names = F, col.names = T, sep="\t", quote=F)
plotHeaytmap=function(plotdata, saveFile){
  rownames(plotdata)=plotdata$region
  plotdata=plotdata[,-1:-4]
  plotdata=plotdata[,c(39:45,27:38, 22:26, 1:21)]
  plotdata=data.frame(t(plotdata))
  methyBreaksList=seq(0, 0.8, by = 0.01)
  pdf(saveFile, width=6, height=5)
  print(pheatmap(plotdata, breaks = methyBreaksList, show_colnames = F,cluster_rows = F, cluster_cols = F, fontsize_row=8, border_color = NA,
                 color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(n = length(methyBreaksList))))
  dev.off()
}
plotHeaytmap(KRT74_regions, "Figure3E/Figure3E_KRT74_methylation.pdf")
plotHeaytmap(SLC2A2_regions, "Figure3E/Figure3E_SLC2A2_methylation.pdf")

################################################Figure3F################################################
load("meta/TCGA_ESCCvsEAC.DESeq2.RData")
ESCC_EAC_result=ESCC_EAC_result[rowMeans(ESCC_EAC_result[,-1:-8])>0.1,]
ESCC_EAC_result=ESCC_EAC_result[!is.na(ESCC_EAC_result$padj),]
ESCC_EAC_targetResults=ESCC_EAC_result[,c(2,6,7,8)]
ESCC_EAC_targetResults$GeneID=rownames(ESCC_EAC_targetResults)

geneInformation=read.table("meta/gencode.v22.all.20190716.txt", sep="\t", header=T, stringsAsFactors = F)
geneInformation=unique(geneInformation[,colnames(geneInformation)%in%c("GeneName","GeneId","Genebiotype")])
geneInformation$GeneId=gsub("\\.[0-9]*", "", geneInformation$GeneId)

ESCC_EAC_targetResults=merge(ESCC_EAC_targetResults, geneInformation[,c(2,1)], by.x="GeneID", by.y="GeneId")
write.table(ESCC_EAC_targetResults, "Figure3F/Figure3F.txt", row.names = F, col.names = T, sep="\t", quote=F)
plotVolcanoPlot=function(plotdata, targets, saveFile){
  plotdata$padj=-log10(plotdata$padj)
  plotdata=as.data.frame(plotdata)
  plotdata1=plotdata[!plotdata$GeneName%in%targets,]
  plotdata2=plotdata[plotdata$GeneName%in%targets,]
  
  p <- ggplot(data=plotdata1,aes(x=log2FoldChange,y=padj))+geom_point(size=1)
  p=p+geom_point(data=plotdata2, aes(x=log2FoldChange,y=padj), size=1.5, color="red")
  p=p+geom_text(data=plotdata2, aes(x=log2FoldChange,y=padj, label=GeneName), size=3, color="red", hjust=0.5)
  p <- p+labs(x = "log2(ESCC tumor v.s. EAC tumor)", y = "-log10 FDR")
  p <- p +theme_classic() + theme(axis.text = element_text(colour = "black",size=9),axis.title=element_text(colour="black",size=12))
  p= p + geom_vline(xintercept = c(-1, 1), linetype="dotted", color = "grey", size=0.8)
  p= p + geom_hline(yintercept = -log10(0.05), linetype="dotted", color = "grey", size=0.8)
  pdf(saveFile, width=5, height=4)
  print(p)
  dev.off()
  
  p <- ggplot(data=plotdata1,aes(x=log2FoldChange,y=padj))+geom_point(size=1)
  p=p+geom_point(data=plotdata2, aes(x=log2FoldChange,y=padj), size=1.5, color="red")
  p <- p+labs(x = "log2(ESCC tumor v.s. EAC tumor)", y = "-log10 FDR")
  p <- p +theme_classic() + theme(axis.text = element_text(colour = "black",size=9),axis.title=element_text(colour="black",size=12))
  p= p + geom_vline(xintercept = c(-1, 1), linetype="dotted", color = "grey", size=0.8)
  p= p + geom_hline(yintercept = -log10(0.05), linetype="dotted", color = "grey", size=0.8)
  png(gsub(".pdf", ".png", saveFile), res=300, width=1500, height=1200)
  print(p)
  dev.off()
}
plotVolcanoPlot(ESCC_EAC_targetResults, c("SLC2A2", "KRT71", "KRT74", "KRT72","KRT73", "KRT2", "KRT1", "KRT77", "KRT79", "KRT78"), "Figure3F/Figure3F.pdf")

################################################Figure4A and Figure 4C################################################
########Figure4A
###analyze ChIP-seq data
##Mapping: Shell/Figure4AC/run.01.bowtie2.Paired.sh (Pair-end sample);  Shell/Figure4AC/run.01.bowtie2.Single.sh (Singel-end sample)
##Call peaks and generate the bw file: Shell/Figure4AC/run.02macs2H3K36me2.sh
##Get H3K36me2 signal for PMD regions: Shell/Figure4AC/getPMDH3K36me2SignalExtend.sh
###result were in "Data/Figure4A/"
plotHistoneProfileExtend=function(inputFileIndex, ymin, ymax, ylab, main, saveFile){
  sortRegions=read.table(paste0(inputFileIndex, ".extend100k.sort.bed"), sep="\t", stringsAsFactors = F)
  sortRegions=sortRegions[,c(1:4,13)]
  colnames(sortRegions)=c("chrom", "start", "end", "region", "type")
  deeptoolsHeatmapHeatmap=read.table(paste(inputFileIndex, ".extend100k.tab",sep=""), nrows=1, comment.char ="",stringsAsFactors = F, sep="\t")
  num1=as.numeric(strsplit(deeptoolsHeatmapHeatmap[1,1],":")[[1]][2])
  num2=as.numeric(strsplit(deeptoolsHeatmapHeatmap[1,2],":")[[1]][2])
  num3=as.numeric(strsplit(deeptoolsHeatmapHeatmap[1,3],":")[[1]][2])
  num4=as.numeric(strsplit(deeptoolsHeatmapHeatmap[1,4],":")[[1]][2])
  start=c(1, (num1+1), (num1+num2+1), (num1+num2+num3+1))
  end=c(num1, (num1+num2), (num1+num2+num3), (num1+num2+num3+num4))
  sortRegions[start[1]:end[1],]$type="EAC_specificPMDs"
  sortRegions[start[2]:end[2],]$type="ESCC_specificPMDs"
  sortRegions[start[3]:end[3],]$type="ESCA_sharedPMDs"
  sortRegions[start[4]:end[4],]$type="ESCA_sharedHMDs"
  
  deeptoolsHeatmapHeatmapHeader=read.table(paste(inputFileIndex,".extend100k.tab",sep=""), skip =2, nrows=1, comment.char ="",stringsAsFactors = F, sep="\t")
  deeptoolsHeatmapHeatmapHeader=deeptoolsHeatmapHeatmapHeader[,-1:-4]
  deeptoolsHeatmapHeatmapHeader=data.frame(type=t(deeptoolsHeatmapHeatmapHeader),stringsAsFactors = F)
  sampleNames=unique(deeptoolsHeatmapHeatmapHeader$type)
  print(sampleNames)
  deeptoolsHeatmapHeatmap=read.table(paste(inputFileIndex,".extend100k.tab",sep=""), skip =3, comment.char ="",stringsAsFactors = F, sep="\t")
  
  pdf(saveFile, width=7, height=7)
  par(mfrow=c(1,length(sampleNames)))
  plotdata=as.data.frame(matrix(numeric(0),ncol=ncol(deeptoolsHeatmapHeatmap)))
  for(i in 1:4){
    temp=data.frame(deeptoolsHeatmapHeatmap[start[i]:end[i], ])
    temp=colMeans(temp, na.rm = T)
    plotdata=rbind(plotdata, temp)
  }
  colnames(plotdata)=paste0("V", 1:ncol(deeptoolsHeatmapHeatmap))
  plotdata=data.frame(t(plotdata))
  par(mar = c(3, 5, 2, 2))
  # plot(plotdata[,1], type="l", col = "#0000F5", lwd = 3,axes=F,  xlab="", ylab=ylab, xaxt = "n", 
  #      ylim=c(ymin, ymax),cex.lab=1.2, cex.axis=1.2)
  # lines(plotdata[,2], type="l", col = "#EA3323", lwd = 3)
  # lines(plotdata[,3], type="l", col = "#FF00FD", lwd = 3)
  # lines(plotdata[,4], type="l", col = "#000000", lwd = 3)
  plot(plotdata[,1], type="l", col = "#ED1D24", lwd = 3,axes=F,  xlab="", ylab=ylab, xaxt = "n", 
       ylim=c(ymin, ymax),cex.lab=1.2, cex.axis=1.2)
  lines(plotdata[,2], type="l", col = "#0827F5", lwd = 3)
  lines(plotdata[,3], type="l", col = "#F6BE00", lwd = 3)
  lines(plotdata[,4], type="l", col = "#00873E", lwd = 3)
  
  axis(1, c(0,100, 200,300), c("+100kb","start","End","-100kb"), cex.axis=1.3)
  axis(2, cex.axis=1.5)
  title(main=main)
  
  # legend(10,ymax,legend=c("EAC_specificPMDs","ESCC_specificPMDs","ESCA_sharedPMDs", "ESCA_sharedHMDs"),pch=16,
  #          col=c("#0000F5","#EA3323","#FF00FD", "#000000"),bty="n", lwd = 4,cex =1.5)
  legend(10,ymax,legend=c("EAC_specificPMDs","ESCC_specificPMDs","ESCA_sharedPMDs", "ESCA_sharedHMDs"),pch=16,
         col=c("#ED1D24","#0827F5","#F6BE00", "#00873E"),bty="n", lwd = 6,cex =1.5)
  dev.off()
  write.table(plotdata, gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
}
plotHistoneProfileExtend("Data/Figure4A/OE19_H3K36me2vsinput.ratio_5k", 0.5, 1.5, "H3K36me2 CPM ratio", "OE19", "Figure4A/Figure4A_OE19.pdf")
plotHistoneProfileExtend("Data/Figure4A/TE5_H3K36me2vsInput.ratio_5k", 0.5, 1.5, "H3K36me2 CPM ratio", "TE5", "Figure4A/Figure4A_TE5.pdf")
plotHistoneProfileExtend("Data/Figure4A/KYSE70_H3K36me2vsInput.ratio_5k", 0.5, 1.5, "H3K36me2 CPM ratio", "KYSE70", "Figure4A/Figure4A_KYSE70.pdf")

########Figure4C
#SRR11657281	 Cal27_H3K36me2 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11657281/SRR11657281
#SRR11657285  Cal27_Input    https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11657285/SRR11657285
#SRR11657318	 Det562_Input   https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11657318/SRR11657318
#SRR11657314  Det562_H3K36me2 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11657314/SRR11657314
#SRR11657346  FaDu_H3K36me2   https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11657346/SRR11657346
#SRR11657348  FaDu_Input     https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR11657348/SRR11657348
###analyze ChIP-seq data
##Mapping: Shell/Figure4AC/run.01.bowtie2.Single.sh (Singel-end sample)
##Call peaks and generate the bw file: Shell/Figure4AC/run.02macs2H3K36me2.sh
##Get H3K36me2 signal for PMD regions: Shell/Figure4AC/getPMDH3K36me2SignalExtend.sh
###result were in "Data/Figure4C/
plotHistoneProfileExtend("Data/Figure4C/Cal27_H3K36me2vsinput.ratio_5k", 0.5, 1.5, "H3K36me2 CPM ratio", "Cal27", "Figure4C/Figure4C_Cal27.pdf")
plotHistoneProfileExtend("Data/Figure4C/Det562_H3K36me2vsinput.ratio_5k", 0.5, 1.5, "H3K36me2 CPM ratio", "Det562", "Figure4C/Figure4C_Det562.pdf")
plotHistoneProfileExtend("Data/Figure4C/FaDu_H3K36me2vsinput.ratio_5k", 0.5, 1.5, "H3K36me2 CPM ratio", "FaDu", "Figure4C/Figure4C_FaDu.pdf")

#########################################################################################################################
#########################################################################################################################
########################################################DMR analysis#####################################################
##Prepare for DMR analysis##
###files in Data/ESCA_cov3_maskUnionPMDs are CpGs removing the union PMDs (Data/MMSeekR_PMDs/ESCA_commonPMDs_union.bed)
step1_merge=function(sampleStart, sampleEnd){
  options(scipen = 20)
  samples=c(paste0("ESCC_",c(1:17,19:22)),paste0("ESCC_Nonmalignant_", 1:3), "ESCC_Nonmalignant_54M", "ESCC_Nonmalignant_53F", 
            paste0("EAC_",c(1:4,6)), paste0("GEJ_",1:7), paste0("GEJ_Nonmalignant_",1:7))
  resultList1=list()
  resultList2=list()
  resultList3=list()
  resultList4=list()
  targetSamples=samples[sampleStart:sampleEnd]
  for(fileName in targetSamples){
    methData=read.table(paste("Data/ESCA_cov3_maskUnionPMDs/", fileName, ".all.sorted.rmblackList.rmPMDs.bed", sep=""),sep="\t",stringsAsFactors=F)
    print(paste("Data/ESCA_cov3_maskUnionPMDs/", fileName, ".all.sorted.rmblackList.rmPMDs.bed", sep=""))
    methData$V2=methData$V2+1
    methData$region=paste(methData$V1,":",methData$V2,sep="")
    methData$methy=round(methData$V5*methData$V4,0)
    methData=methData[,c(6,1:2,7,5,4)]
    colnames(methData)=c("region","chr","pos",paste("methy.", fileName, sep=""), paste("coverage.", fileName, sep=""), paste("beta.", fileName,sep=""))
    methData1=methData[methData$chr%in%paste0("chr",1:5),]
    resultList1[[fileName]]=methData1
    methData2=methData[methData$chr%in%paste0("chr",6:10),]
    resultList2[[fileName]]=methData2
    methData3=methData[methData$chr%in%paste0("chr",11:15),]
    resultList3[[fileName]]=methData3
    methData4=methData[methData$chr%in%paste0("chr",16:22),]
    resultList4[[fileName]]=methData4
  }
  save(resultList1, file=paste("Data/MaskUnionPMDs_DMRs/Sample_", sampleStart, "_",sampleEnd, ".chr1_5.RData",sep=""))
  save(resultList2, file=paste("Data/MaskUnionPMDs_DMRs/Sample_", sampleStart, "_",sampleEnd, ".chr6_10.RData",sep=""))
  save(resultList3, file=paste("Data/MaskUnionPMDs_DMRs/Sample_", sampleStart, "_",sampleEnd, ".chr11_15.RData",sep=""))
  save(resultList4, file=paste("Data/MaskUnionPMDs_DMRs/Sample_", sampleStart, "_",sampleEnd, ".chr16_22.RData",sep=""))
  
  mergeMatrix = resultList1[[targetSamples[1]]]
  for(i in 2:length(targetSamples)){
    mergeMatrix=merge(mergeMatrix, resultList1[[targetSamples[i]]], by.x=c("region","chr","pos"), by.y=c("region","chr","pos"), all=TRUE)
  }
  save(mergeMatrix, file=paste("Data/MaskUnionPMDs_DMRs/Sample_", sampleStart,"_",sampleEnd, "_merge.chr1_5.RData",sep=""))
  
  mergeMatrix = resultList2[[targetSamples[1]]]
  for(i in 2:length(targetSamples)){
    mergeMatrix=merge(mergeMatrix, resultList2[[targetSamples[i]]], by.x=c("region","chr","pos"), by.y=c("region","chr","pos"), all=TRUE)
  }
  save(mergeMatrix, file=paste("Data/MaskUnionPMDs_DMRs/Sample_", sampleStart,"_",sampleEnd, "_merge.chr6_10.RData",sep=""))
  
  mergeMatrix = resultList3[[targetSamples[1]]]
  for(i in 2:length(targetSamples)){
    mergeMatrix=merge(mergeMatrix, resultList3[[targetSamples[i]]], by.x=c("region","chr","pos"), by.y=c("region","chr","pos"), all=TRUE)
  }
  save(mergeMatrix, file=paste("Data/MaskUnionPMDs_DMRs/Sample_", sampleStart,"_",sampleEnd, "_merge.chr11_15.RData",sep=""))
  
  mergeMatrix = resultList4[[targetSamples[1]]]
  for(i in 2:length(targetSamples)){
    mergeMatrix=merge(mergeMatrix, resultList4[[targetSamples[i]]], by.x=c("region","chr","pos"), by.y=c("region","chr","pos"), all=TRUE)
  }
  save(mergeMatrix, file=paste("Data/MaskUnionPMDs_DMRs/Sample_", sampleStart,"_",sampleEnd, "_merge.chr16_22.RData",sep=""))
}
step1_merge(1, 10)
step1_merge(11, 21)
step1_merge(22, 26)
step1_merge(39, 45)
step1_merge(27, 38)

step2_merge=function(chrStart, chrEnd){
  options(scipen = 20)
  load(paste0("Data/MaskUnionPMDs_DMRs/Sample_1_10_merge.chr", chrStart, "_", chrEnd, ".RData"))
  print(paste0("Data/MaskUnionPMDs_DMRs/Sample_1_10_merge.chr", chrStart, "_", chrEnd, ".RData"))
  mergeMatrixAll=mergeMatrix
  load(paste0("Data/MaskUnionPMDs_DMRs/Sample_11_21_merge.chr", chrStart, "_", chrEnd, ".RData"))
  print(paste0("Data/MaskUnionPMDs_DMRs/Sample_11_21_merge.chr", chrStart, "_", chrEnd, ".RData"))
  mergeMatrixAll=merge(mergeMatrixAll, mergeMatrix, by.x=c("region","chr","pos"), by.y=c("region","chr","pos"), all=TRUE)
  load(paste0("Data/MaskUnionPMDs_DMRs/Sample_22_26_merge.chr",chrStart, "_", chrEnd, ".RData"))
  print(paste0("Data/Sample_22_26_merge.chr",chrStart, "_", chrEnd, ".RData"))
  mergeMatrixAll=merge(mergeMatrixAll, mergeMatrix, by.x=c("region","chr","pos"), by.y=c("region","chr","pos"), all=TRUE)
  load(paste0("Data/MaskUnionPMDs_DMRs/Sample_27_38_merge.chr", chrStart, "_", chrEnd, ".RData"))
  print(paste0("Data/MaskUnionPMDs_DMRs/Sample_27_38_merge.chr", chrStart, "_", chrEnd, ".RData"))
  mergeMatrixAll=merge(mergeMatrixAll, mergeMatrix, by.x=c("region","chr","pos"), by.y=c("region","chr","pos"), all=TRUE)
  load(paste0("Data/MaskUnionPMDs_DMRs/Sample_39_45_merge.chr", chrStart, "_", chrEnd, ".RData"))
  print(paste0("Data/MaskUnionPMDs_DMRs/Sample_39_45_merge.chr", chrStart, "_", chrEnd, ".RData"))
  mergeMatrixAll=merge(mergeMatrixAll, mergeMatrix, by.x=c("region","chr","pos"), by.y=c("region","chr","pos"), all=TRUE)
  save(mergeMatrixAll, file=paste0("Data/MaskUnionPMDs_DMRs/AllSample_merge.chr", chrStart, "_", chrEnd, ".RData"))
}
step2_merge(1, 5)
step2_merge(6, 10)
step2_merge(11, 15)
step2_merge(16, 22)

step3_splitForDMR=function(chrStart, chrEnd){
  options(scipen = 20)
  targetChromosome=paste("chr",chrStart:chrEnd, sep="")
  load(paste0("Data/MaskUnionPMDs_DMRs/AllSample_merge.chr", chrStart, "_", chrEnd, ".RData"))
  for(chrInfo in targetChromosome){
    print(chrInfo)
    targetMatrix=mergeMatrixAll[mergeMatrixAll$chr%in%chrInfo,]
    targetMatrix[is.na(targetMatrix)]=0
    chr=targetMatrix$chr
    pos=targetMatrix$pos
    allCov=targetMatrix[,grep("coverage.",colnames(targetMatrix))]
    colnames(allCov)=gsub("coverage.","",colnames(allCov))
    allBeta=targetMatrix[,grep("beta.",colnames(targetMatrix))]
    colnames(allBeta)=gsub("beta.","",colnames(allBeta))
    allMethy=targetMatrix[,grep("methy.",colnames(targetMatrix))]
    colnames(allMethy)=gsub("methy.","",colnames(allMethy))
    allSamples=colnames(allBeta)
    allType=c(rep("ESCC_Tumor",21),rep("ESCC_Nonmalignant",5),rep("EAC_Tumor",5), rep("GEJ_Tumor",7), rep("GEJ_Nonmalignant",7))
    
    ####ESCC Tumor vs EAC Tumor######
    cov=allCov[,which(allType%in%c("ESCC_Tumor","EAC_Tumor","GEJ_Tumor"))]
    methy=allMethy[,which(allType%in%c("ESCC_Tumor","EAC_Tumor","GEJ_Tumor"))]
    replicate=c(seq(1:21),seq(1:12))
    type=c(rep("ESCC_Tumor",21),rep("EAC_Tumor",12))
    samples=colnames(methy)
    print("Saving File....")
    save(chr, pos, cov, methy, samples, type, replicate, file=paste("MaskUnionPMDs_DMRs/ESCCTumorvsEACTumor/", chrInfo,"_ESCCTumorvsEACTumor.Mask.RData",sep=""))
  }
}
step3_splitForDMR(1, 5)
step3_splitForDMR(6, 10)
step3_splitForDMR(11, 15)
step3_splitForDMR(16, 22)

step4_splitForDMR=function(){
  options(scipen = 20)
  load("MaskUnionPMDs_DMRs/AllSample_merge.chr1_5.RData")
  allCov=mergeMatrixAll[,grep("coverage.",colnames(mergeMatrixAll))]
  colnames(allCov)=gsub("coverage.","",colnames(allCov))
  allBeta=mergeMatrixAll[,grep("beta.",colnames(mergeMatrixAll))]
  colnames(allBeta)=gsub("beta.","",colnames(allBeta))
  betaMatrixTemp=data.frame(mergeMatrixAll[,1:3], allBeta,stringsAsFactors = F)
  betaMatrix=betaMatrixTemp
  keep7=mergeMatrixAll[rowSums(allCov>=7, na.rm = T)==ncol(allBeta),]$region
  
  load("MaskUnionPMDs_DMRs/AllSample_merge.chr6_10.RData")
  allCov=mergeMatrixAll[,grep("coverage.",colnames(mergeMatrixAll))]
  colnames(allCov)=gsub("coverage.","",colnames(allCov))
  allBeta=mergeMatrixAll[,grep("beta.",colnames(mergeMatrixAll))]
  colnames(allBeta)=gsub("beta.","",colnames(allBeta))
  betaMatrixTemp=data.frame(mergeMatrixAll[,1:3], allBeta,stringsAsFactors = F)
  betaMatrix=rbind(betaMatrix, betaMatrixTemp)
  keep7=c(keep7, mergeMatrixAll[rowSums(allCov>=7, na.rm = T)==ncol(allBeta),]$region)
  
  load("MaskUnionPMDs_DMRs/AllSample_merge.chr11_15.RData")
  allCov=mergeMatrixAll[,grep("coverage.",colnames(mergeMatrixAll))]
  colnames(allCov)=gsub("coverage.","",colnames(allCov))
  allBeta=mergeMatrixAll[,grep("beta.",colnames(mergeMatrixAll))]
  colnames(allBeta)=gsub("beta.","",colnames(allBeta))
  betaMatrixTemp=data.frame(mergeMatrixAll[,1:3], allBeta,stringsAsFactors = F)
  betaMatrix=rbind(betaMatrix, betaMatrixTemp)
  keep7=c(keep7, mergeMatrixAll[rowSums(allCov>=7, na.rm = T)==ncol(allBeta),]$region)
  
  load("MaskUnionPMDs_DMRs/AllSample_merge.chr16_22.RData")
  allCov=mergeMatrixAll[,grep("coverage.",colnames(mergeMatrixAll))]
  colnames(allCov)=gsub("coverage.","",colnames(allCov))
  allBeta=mergeMatrixAll[,grep("beta.",colnames(mergeMatrixAll))]
  colnames(allBeta)=gsub("beta.","",colnames(allBeta))
  betaMatrixTemp=data.frame(mergeMatrixAll[,1:3], allBeta,stringsAsFactors = F)
  betaMatrix=rbind(betaMatrix, betaMatrixTemp)
  keep7=c(keep7, mergeMatrixAll[rowSums(allCov>=7, na.rm = T)==ncol(allBeta),]$region)
  
  save(betaMatrix, file="MaskUnionPMDs_DMRs/AllSample_beta.maskPMD.RData")
  betaMatrix_cov7=betaMatrix[betaMatrix$region%in%keep7,]
  save(betaMatrix_cov7, file="MaskUnionPMDs_DMRs/AllSample_beta_cov7.maskPMD.RData")
}
step4_splitForDMR()

###runDMR for ESCCTumor vs EACTumor
runDMR=function(seed, maxPerms, chrInfo, cutoff){
  options(scipen = 20)
  options(MulticoreParam=quote(MulticoreParam(workers=1)))
  library(dmrseq)
  library(bsseq)
  load(paste("Data/MaskUnionPMDs_DMRs/ESCCTumorvsEACTumor/",chrInfo, "_ESCCTumorvsEACTumor.Mask.RData", sep=""))
  bs = BSseq(chr = chr, pos = pos, M = as.matrix(methy), Cov = as.matrix(cov), sampleNames = samples)
  type=c(rep("ESCC_Tumor",21),rep("EAC_Tumor",12))
  pData(bs)$Type <- type
  pData(bs)$Replicate <- replicate
  loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs, type="Cov")>=3) >= (length(samples)*0.95))
  print(length(loci.idx))
  sample.idx <- which(pData(bs)$Type %in% c("EAC_Tumor", "ESCC_Tumor"))
  bs.filtered <- bs[loci.idx, sample.idx]
  bs.filtered=sort(bs.filtered)
  set.seed(seed)
  regions <- dmrseq(bs=bs.filtered, testCovariate="Type", cutoff =cutoff, bpSpan=1000, minInSpan=30, maxPerms=maxPerms)
  save(regions, file=paste("Data/MaskUnionPMDs_DMRs/ESCCTumorvsEACTumor/Results_20201121_500_0.1/",chrInfo, "_ESCCTumorvsEACTumor_DMR_Mask.RData", sep=""))
}
for(chr in paste0("chr", 1:22)){
  runDMR(20201121, 500, chr, 0.1)
}

###Analysis DMR results
dmrType="Data/MaskUnionPMDs_DMRs/ESCCTumorvsEACTumor/"
name="ESCCTumorvsEACTumor_DMR_Mask"
name2="ESCCTumorvsEACTumor.Mask"
group1Name="EAC_Tumor"
group2Name="ESCC_Tumor"
resultPath="Results_20201121_500_0.1"
DMRPlotPath="DMRPlot_20201121_500_0.1"

####
DMRList=as.data.frame(matrix(numeric(0),ncol=16))
DMRProbeList=as.data.frame(matrix(numeric(0),ncol=35))
DMRAllProbeList=as.data.frame(matrix(numeric(0),ncol=35))
DMRAllList=as.data.frame(matrix(numeric(0),ncol=16))
for(chrInfo in chrListTarget){
  filename=paste(dmrType, "/", resultPath ,"/",chrInfo,"_",name,".RData",sep="")
  if(file.exists(filename)){
    print(chrInfo)
    load(filename)
    regionsMatrix=as.data.frame(regions)
    
    load(paste(dmrType, "/",chrInfo, "_", name2, ".RData", sep=""))
    bs = BSseq(chr = chr, pos = pos, M = as.matrix(methy), Cov = as.matrix(cov), sampleNames = samples)
    pData(bs)$Type <- type
    pData(bs)$Replicate <- replicate
    loci.idx <- which(rowSums2(getCoverage(bs, type="Cov")>=3) >= (length(samples)*0.95))
    print(length(loci.idx))
    sample.idx <- which(pData(bs)$Type %in% c(group1Name, group2Name))
    bs.filtered <- bs[loci.idx, ]
    bs.filtered=sort(bs.filtered)
    
    rawDiff <- meanDiff(bs.filtered, dmrs=regions, testCovariate="Type")
    regionsMatrix$deltaMean=rawDiff
    number = sum(regionsMatrix$qval < 0.05)
    print(paste("Total regions:", nrow(regionsMatrix), "   Significant regions:", number, sep=""))
    
    result=as.data.frame(regionsMatrix)
    result$DMR=paste(result$seqnames,":",result$start,"-",result$end,sep="")
    result_methylation = lapply(1:nrow(result), function(x){
      term=result[x,]
      target = as.data.frame(bs.filtered@rowRanges[(term$index.start):(term$index.end),])
      target$DMR=term$DMR
      return(target)})
    result_methylation=do.call(rbind, result_methylation)
    DMRAllProbeList=rbind(DMRAllProbeList, result_methylation)
    if(number>0){
      sigRegions <- regionsMatrix[regionsMatrix$qval < 0.05,]
      ###the proportion of regions with hyper-methylation
      result=as.data.frame(sigRegions)
      result$DMR=paste(result$seqnames,":",result$start,"-",result$end,sep="")
      DMRList=rbind(DMRList,result)
      result_methylation = lapply(1:nrow(result), function(x){
        term=result[x,]
        target = as.data.frame(bs.filtered@rowRanges[(term$index.start):(term$index.end),])
        target$DMR=term$DMR
        return(target)})
      result_methylation=do.call(rbind, result_methylation)
      DMRProbeList=rbind(DMRProbeList, result_methylation)
    }
    DMRAllList=rbind(DMRAllList, regionsMatrix)
  }
}
DMRProbeList$region=paste(DMRProbeList$seqnames, ":", DMRProbeList$start, sep="")
DMRProbeList=DMRProbeList[,colnames(DMRProbeList)%in%c("DMR","region")]
DMRAllProbeList$region=paste(DMRAllProbeList$seqnames, ":", DMRAllProbeList$start, sep="")
DMRAllProbeList=DMRAllProbeList[,colnames(DMRAllProbeList)%in%c("DMR","region")]
save(DMRAllList, DMRAllProbeList, DMRList, DMRProbeList, file=paste(dmrType, "/", DMRPlotPath,"/", basename(dmrType), ".DMR.RData",sep=""))

load("Data/MaskUnionPMDs_DMRs/AllSample_beta.maskPMD.RData")
DMRProbeMatrix=betaMatrix[betaMatrix$region%in%DMRProbeList$region,]
DMRProbeMatrix$region=factor(DMRProbeMatrix$region, levels=DMRProbeList$region)
DMRProbeMatrix=DMRProbeMatrix[order(DMRProbeMatrix$region),]
DMRProbeMatrix$DMR=DMRProbeList$DMR
DMRProbeMatrix=DMRProbeMatrix[,-1:-3]
DMRMethyMatrix=lapply(unique(DMRProbeMatrix$DMR),function(x){
  colMeans(DMRProbeMatrix[DMRProbeMatrix$DMR%in%x,-ncol(DMRProbeMatrix)], na.rm = T)
})
DMRMethyMatrix=do.call(rbind,DMRMethyMatrix)
rownames(DMRMethyMatrix)=unique(DMRProbeMatrix$DMR)
DMRMethyMatrix[is.na(DMRMethyMatrix)]=0
save(DMRMethyMatrix,file=paste(dmrType, "/", DMRPlotPath, "/", basename(dmrType), ".DMRBetaMatrix.RData",sep=""))

DMRAllProbeMatrix=betaMatrix[betaMatrix$region%in%DMRAllProbeList$region,]
rm(betaMatrix)
DMRAllProbeMatrix$region=factor(DMRAllProbeMatrix$region, levels=DMRAllProbeList$region)
DMRAllProbeMatrix=DMRAllProbeMatrix[order(DMRAllProbeMatrix$region),]
DMRAllProbeMatrix$DMR=DMRAllProbeList$DMR
DMRAllProbeMatrix=DMRAllProbeMatrix[,-1:-3]
AllRegionMethyMatrix=lapply(unique(DMRAllProbeMatrix$DMR),function(x){
  colMeans(DMRAllProbeMatrix[DMRAllProbeMatrix$DMR%in%x,-ncol(DMRAllProbeMatrix)], na.rm=T)
})
AllRegionMethyMatrix=do.call(rbind,AllRegionMethyMatrix)
rownames(AllRegionMethyMatrix)=unique(DMRAllProbeMatrix$DMR)
colnames(AllRegionMethyMatrix)[39:45]=paste("GEJ_Nonmalignant_",1:7, sep="")
AllRegionMethyMatrix[is.na(AllRegionMethyMatrix)]=0
save(AllRegionMethyMatrix,file=paste(dmrType, "/", DMRPlotPath, "/", basename(dmrType), ".AllRegionBetaMatrix.RData",sep=""))


########################################################Figure5A, 6A, S4A#####################################################
load(paste(dmrType, "/", DMRPlotPath,"/", basename(dmrType), ".DMR.RData",sep=""))
load(paste(dmrType, "/", DMRPlotPath, "/", basename(dmrType), ".DMRBetaMatrix.RData",sep=""))
load(paste(dmrType, "/", DMRPlotPath, "/", basename(dmrType), ".AllRegionBetaMatrix.RData",sep=""))

if(group1Name%in%"EAC_Tumor"){
  sampleType1="EAC/GEJ_Tumor"
}else{
  sampleType1=group1Name
}
if(group2Name%in%"EAC_Tumor"){
  sampleType2="EAC/GEJ_Tumor"
}else{
  sampleType2=group2Name
}

AllRegionMethyMatrix_delta=rowMeans(AllRegionMethyMatrix[,colnames(AllRegionMethyMatrix)%in%rownames(annotation_col[annotation_col$Type%in%sampleType1,])])-
  rowMeans(AllRegionMethyMatrix[,colnames(AllRegionMethyMatrix)%in%rownames(annotation_col[annotation_col$Type%in%sampleType2,])])
DMRAllList$DMR=paste0(DMRAllList$seqnames, ":", DMRAllList$start, "-",DMRAllList$end)
DMRAllList$deltaMean2=AllRegionMethyMatrix_delta
DMRPlot=function(dmrResultList, path, group1Name, group2Name, saveFileIndex){
  print("plot methylation volcano plot")
  plotdata = data.frame(chr=dmrResultList$seqnames, stat=dmrResultList$stat, beta=dmrResultList$deltaMean2, FDR=(-log10(dmrResultList$qval)), region=DMRAllList$DMR)
  plotdata$Type="No significant change"
  group1=paste("Hypomethylation in ", group1Name, sep="")
  group2=paste("Hypomethylation in ", group2Name, sep="")
  plotdata[plotdata$stat>0&plotdata$FDR>(-log10(0.05))&abs(plotdata$beta)>0.2,]$Type=group1
  plotdata[plotdata$stat<0&plotdata$FDR>(-log10(0.05))&abs(plotdata$beta)>0.2,]$Type=group2
  plotdata$Type=factor(plotdata$Type,levels=c("No significant change",group1,group2))
  print(table(plotdata[,c(1,6)]))
  write.table(plotdata, paste0(saveFileIndex, ".txt"), row.names = F, col.names = T, sep="\t", quote=F)
  p <- ggplot(data=plotdata,aes(x=beta,y=FDR,colour=Type))+geom_point(size=0.5)
  p <- p+labs(x = "Methylation difference", y = "-log10 FDR")+
    scale_color_manual(values=c('black','#0000F5',"#EA3323"),
                       labels=c(paste("No-change (",format(table(plotdata$Type)[[1]], big.mark=",",scientific=FALSE),")",sep=""),
                                paste(group1, " (",format(table(plotdata$Type)[names(table(plotdata$Type))%in%group1][[1]],big.mark=",",scientific=FALSE),")",sep=""),
                                paste(group2, " (",format(table(plotdata$Type)[names(table(plotdata$Type))%in%group2][[1]],big.mark=",",scientific=FALSE),")",sep="")))
  p <- p +theme_classic() + theme(axis.text = element_text(colour = "black",size=9),axis.title=element_text(colour="black",size=12))
  p = p + geom_hline(yintercept=(-log10(0.05)), linetype="dashed", color = "red", size=0.5)
  p = p + geom_vline(xintercept = c(-0.2, 0.2), linetype="dashed", color = "red", size=0.5)
  p=p+scale_x_continuous(limits = c(-0.7,0.7), label=round(seq(-0.7,0.7,by=0.1),1), breaks=seq(-0.7,0.7,by=0.1))
  pdf(paste(saveFileIndex, "_VolcanoPlot.pdf",sep=""),height=4, width=8)
  print(p)
  dev.off()
  png(paste(saveFileIndex, "_VolcanoPlot.png",sep=""),height=1500, width=2600, res=300)
  print(p)
  dev.off()
}
DMRPlot(DMRAllList, paste(dmrType,"/", DMRPlotPath, "/",sep=""), group1Name, group2Name, "Figure5A/Figure5A_ESCCvsEAC")

####get EAC and ESCC tumor DMRs
DMRList=DMRAllList[DMRAllList$qval<0.05&abs(DMRAllList$deltaMean2)>0.2,]
hypoGroup1=DMRList[DMRList$stat>0&abs(DMRList$deltaMean2)>0.2,]
write.table(hypoGroup1[,1:3], file=paste("Data/MaskUnionPMDs_DMRs/hypo", group1Name, ".bed", sep=""),sep="\t",quote=F, row.names = F, col.names = F)
hypoGroup2=DMRList[DMRList$stat<0&abs(DMRList$deltaMean2)>0.2,]
write.table(hypoGroup2[,1:3], file=paste("Data/MaskUnionPMDs_DMRs/hypo", group2Name, ".bed", sep=""),sep="\t",quote=F, row.names = F, col.names = F)

###get subtype-specific DMRs
hypogroup1BetaMatrix=DMRMethyMatrix[rownames(DMRMethyMatrix)%in%hypoGroup1$DMR,]
hypogroup2BetaMatrix=DMRMethyMatrix[rownames(DMRMethyMatrix)%in%hypoGroup2$DMR,]
calculatePvalue=function(dataMatrix, type1, type2){
  TtestResult1= lapply(1:nrow(dataMatrix), function(x){
    t.test(dataMatrix[x,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%type1,])],
           dataMatrix[x,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%type2,])], alternative ="less")$p.value
  })
  TtestResult1=do.call(c, TtestResult1)
  
  TtestResult2=lapply(1:nrow(dataMatrix), function(x){
    t.test(dataMatrix[x,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%type1,])],
           dataMatrix[x,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%type2,])], alternative ="great")$p.value
  })
  TtestResult2=do.call(c, TtestResult2)
  
  TtestResult=data.frame(lessPvalue=TtestResult1, greatPvale=TtestResult2)
  TtestResult$lessFDR=p.adjust(TtestResult$lessPvalue, method="BH")
  TtestResult$greatFDR=p.adjust(TtestResult$greatPvale, method="BH")
  rownames(TtestResult)=rownames(dataMatrix)
  return(TtestResult)
}
hypogroup1BetaMatrix_Ttest=calculatePvalue(hypogroup1BetaMatrix, sampleType1, "GEJ_Nonmalignant")
hypogroup2BetaMatrix_Ttest=calculatePvalue(hypogroup2BetaMatrix, sampleType2, "ESCC_Nonmalignant")
specific_hypoEAC=hypoGroup1[hypoGroup1$DMR%in%rownames(hypogroup1BetaMatrix_Ttest[hypogroup1BetaMatrix_Ttest$lessFDR<0.05,]),1:3]
specific_hypoESCC=hypoGroup2[hypoGroup2$DMR%in%rownames(hypogroup2BetaMatrix_Ttest[hypogroup2BetaMatrix_Ttest$lessFDR<0.05,]),1:3]
write.table(specific_hypoEAC, file=paste("Data/MaskUnionPMDs_DMRs/specific_hypo_",group1Name, ".bed", sep=""), quote=F, sep="\t", row.names = F, col.names = F)
write.table(specific_hypoESCC, file=paste("Data/MaskUnionPMDs_DMRs/specific_hypo_",group2Name, ".bed", sep=""), quote=F, sep="\t", row.names = F, col.names = F)

plotAllDMRMethylationDefaultOrder=function(dataMatrix){
  print(paste("DMR regions:", nrow(dataMatrix), sep=" "))
  plotdata_ESCCTumor=dataMatrix[,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%c("ESCC_Tumor"),])]
  col_order=hclust(dist(t(plotdata_ESCCTumor)))$order
  plotdata_ESCCTumor=plotdata_ESCCTumor[,col_order]
  
  plotdata_EACTumor=dataMatrix[,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%c("EAC/GEJ_Tumor"),])]
  col_order=hclust(dist(t(plotdata_EACTumor)))$order
  plotdata_EACTumor=plotdata_EACTumor[,col_order]
  
  plotdata_ESCCNormal=dataMatrix[,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%c("ESCC_Nonmalignant"),])]
  col_order=hclust(dist(t(plotdata_ESCCNormal)))$order
  plotdata_ESCCNormal=plotdata_ESCCNormal[,col_order]
  
  plotdata_EACNormal=dataMatrix[,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%c("GEJ_Nonmalignant"),])]
  col_order=hclust(dist(t(plotdata_EACNormal)))$order
  plotdata_EACNormal=plotdata_EACNormal[,col_order]
  
  temp=as.data.frame(matrix(rep(NA, length=nrow(plotdata_ESCCNormal)),ncol=1))
  if((length(grep("Tumor", c(sampleType1, sampleType2)))>0)&(length(grep("Nonmalignant", c(sampleType1, sampleType2)))>0)){
    plotdata=cbind(plotdata_ESCCNormal,plotdata_ESCCTumor, temp, plotdata_EACTumor,plotdata_EACNormal)
  }else{
    plotdata=cbind(plotdata_ESCCNormal,temp, plotdata_ESCCTumor, plotdata_EACTumor,temp, plotdata_EACNormal)
  }
  ann_colors2=ann_colors
  ann_colors2$Type=ann_colors2$Type[names(ann_colors2$Type)%in%as.character(unique(annotation_col$Type))]
  resultHeatmap=pheatmap(plotdata, annotation_col = annotation_col[,colnames(annotation_col)%in%"Type", drop=F], annotation_colors = ann_colors2,
                         color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(200),
                         show_rownames = F, show_colnames = T, cluster_rows = F, cluster_cols = F,na_col="white")
  return(resultHeatmap)
}
plotTargetDeltaHeatmap2=function(dataMatrix, ylim){
  plotdata=data.frame(delta=dataMatrix)
  colnames(plotdata)="delta"
  plotdata$order=nrow(plotdata):1
  p=ggplot(plotdata,aes(x=order,y=delta))+geom_bar(stat="identity",width=1)+ylim(-ylim, ylim)+theme_classic()+ylab("")+xlab("")
  p=p+ coord_flip()
  return(p)
}
plotFDRHeatmap=function(dataMatrix){
  plotdata=data.frame(delta=dataMatrix)
  resultHeatmap=pheatmap(plotdata, show_rownames = F, show_colnames = T, cluster_rows = F, cluster_cols = F,color = c("grey", "darkred"))
  return(resultHeatmap)
}
plotTumorSpecificDMRsHeatmap=function(dataMatrix, hypogroupBetaMatrix_Ttest, type1, type2, saveFile){
  type1Mean=rowMeans(dataMatrix[,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%type1,])], na.rm = T)
  type2Mean=rowMeans(dataMatrix[,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%type2,])], na.rm = T)
  deltaMean=type1Mean-type2Mean
  deltaMean=data.frame(mean=deltaMean)
  deltaMean=deltaMean[order(deltaMean$mean),,drop=F]
  rowOrder=factor(rownames(dataMatrix), levels=rownames(deltaMean))
  dataMatrix=dataMatrix[order(rowOrder),]
  
  
  FDRData=hypogroupBetaMatrix_Ttest[,c(3,4)]
  FDRData$lessFDRType=0
  FDRData[FDRData$lessFDR<0.05,]$lessFDRType=1
  FDRData$greatFDRType=0
  FDRData[FDRData$greatFDR<0.05,]$greatFDRType=1
  FDRData=FDRData[order(rowOrder),]
  
  heatmapResult=plotAllDMRMethylationDefaultOrder(dataMatrix)
  deltaHeatmap=plotTargetDeltaHeatmap2(deltaMean)
  fdrHeatmap=plotFDRHeatmap(FDRData[,3,drop=F])
  dev.off()
  plot_list=list()
  plot_list[[1]]=deltaHeatmap
  plot_list[[2]]=fdrHeatmap[[4]]
  plot_list[[3]]=heatmapResult[[4]]
  pdf(saveFile, width = 12, height=8)
  print(grid.arrange(arrangeGrob(grobs= plot_list, widths = c(2.5/20,2.5/20, 7/10), ncol=3)))
  dev.off()
  
  dataMatrix=data.frame(region=rownames(dataMatrix), dataMatrix)
  write.table(dataMatrix, gsub(".pdf", "_matrix.txt", saveFile), quote=F, row.names = F, col.names = T, sep="\t")
  FDRData=data.frame(region=rownames(FDRData), FDRData)
  write.table(FDRData, gsub(".pdf", "_fdr.txt", saveFile), quote=F, row.names = F, col.names = T, sep="\t")
  deltaMean=data.frame(region=rownames(deltaMean), deltaMean)
  write.table(deltaMean, gsub(".pdf", "_deltaMean.txt", saveFile), quote=F, row.names = F, col.names = T, sep="\t")
}
plotTumorSpecificDMRsHeatmap(hypogroup1BetaMatrix, hypogroup1BetaMatrix_Ttest, sampleType1, "GEJ_Nonmalignant", "Figure6A/Figure6A.heatmap.pdf")
plotTumorSpecificDMRsHeatmap(hypogroup2BetaMatrix, hypogroup2BetaMatrix_Ttest, sampleType2, "ESCC_Nonmalignant", "FigureS5A/FigureS5A.heatmap.pdf")

nonSpecific_hypoEAC=hypoGroup1[!hypoGroup1$DMR%in%rownames(hypogroup1BetaMatrix_Ttest[hypogroup1BetaMatrix_Ttest$lessFDR<0.05,]),1:3]
nonSpecific_hypoESCC=hypoGroup2[!hypoGroup2$DMR%in%rownames(hypogroup2BetaMatrix_Ttest[hypogroup2BetaMatrix_Ttest$lessFDR<0.05,]),1:3]
write.table(nonSpecific_hypoEAC, file=paste("Data/MaskUnionPMDs_DMRs/nonSpecific_hypo_",group1Name, ".bed", sep=""), quote=F, sep="\t", row.names = F, col.names = F)
write.table(nonSpecific_hypoESCC, file=paste("Data/MaskUnionPMDs_DMRs/nonSpecific_hypo_",group2Name, ".bed", sep=""), quote=F, sep="\t", row.names = F, col.names = F)

#########################################################FigureS6#####################
hypergroup1BetaMatrix=hypogroup2BetaMatrix
hypergroup2BetaMatrix=hypogroup1BetaMatrix
hypergroup1BetaMatrix_Ttest=calculatePvalue(hypergroup1BetaMatrix, sampleType1, "GEJ_Nonmalignant")
hypergroup2BetaMatrix_Ttest=calculatePvalue(hypergroup2BetaMatrix, sampleType2, "ESCC_Nonmalignant")

plotTumorSpecificDMRsHeatmap2=function(dataMatrix, hypergroupBetaMatrix_Ttest, type1, type2, saveFile){
  type1Mean=rowMeans(dataMatrix[,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%type1,])], na.rm = T)
  type2Mean=rowMeans(dataMatrix[,colnames(dataMatrix)%in%rownames(annotation_col[annotation_col$Type%in%type2,])], na.rm = T)
  deltaMean=type1Mean-type2Mean
  deltaMean=data.frame(mean=deltaMean)
  deltaMean=deltaMean[order(deltaMean$mean,decreasing = T),,drop=F]
  rowOrder=factor(rownames(dataMatrix), levels=rownames(deltaMean))
  dataMatrix=dataMatrix[order(rowOrder),]
  
  FDRData=hypergroupBetaMatrix_Ttest[,c(3,4)]
  FDRData$lessFDRType=0
  FDRData[FDRData$lessFDR<0.05,]$lessFDRType=1
  FDRData$greatFDRType=0
  FDRData[FDRData$greatFDR<0.05,]$greatFDRType=1
  FDRData=FDRData[order(rowOrder),]
  print(nrow(FDRData[FDRData$greatFDRType==1,]))
  
  heatmapResult=plotAllDMRMethylationDefaultOrder(dataMatrix)
  deltaHeatmap=plotTargetDeltaHeatmap2(deltaMean, 1)
  fdrHeatmap=plotFDRHeatmap(FDRData[,4,drop=F])
  dev.off()
  plot_list=list()
  plot_list[[1]]=deltaHeatmap
  plot_list[[2]]=fdrHeatmap[[4]]
  plot_list[[3]]=heatmapResult[[4]]
  pdf(saveFile, width = 12, height=8)
  print(grid.arrange(arrangeGrob(grobs= plot_list, widths = c(2.5/20,2.5/20, 7/10), ncol=3)))
  dev.off()
  
  dataMatrix=data.frame(region=rownames(dataMatrix), dataMatrix)
  write.table(dataMatrix, gsub(".pdf", "_matrix.txt", saveFile), quote=F, row.names = F, col.names = T, sep="\t")
  FDRData=data.frame(region=rownames(FDRData), FDRData)
  write.table(FDRData, gsub(".pdf", "_fdr.txt", saveFile), quote=F, row.names = F, col.names = T, sep="\t")
  deltaMean=data.frame(region=rownames(deltaMean), deltaMean)
  write.table(deltaMean, gsub(".pdf", "_deltaMean.txt", saveFile), quote=F, row.names = F, col.names = T, sep="\t")
}
plotTumorSpecificDMRsHeatmap2(hypergroup1BetaMatrix, hypergroup1BetaMatrix_Ttest, sampleType1, "GEJ_Nonmalignant", "FigureS6/FigureS6A.EACTumor_tsDMRS.hyper.heatmap.pdf")
plotTumorSpecificDMRsHeatmap2(hypergroup2BetaMatrix, hypergroup2BetaMatrix_Ttest, sampleType2, "ESCC_Nonmalignant", "FigureS6/FigureS6B.ESCCTumor_tsDMRS.hyper.heatmap.pdf")

specific_hyperEAC=hypoGroup2[hypoGroup2$DMR%in%rownames(hypergroup1BetaMatrix_Ttest[hypergroup1BetaMatrix_Ttest$greatFDR<0.05,]),1:3]
specific_hyperESCC=hypoGroup1[hypoGroup1$DMR%in%rownames(hypergroup2BetaMatrix_Ttest[hypergroup2BetaMatrix_Ttest$greatFDR<0.05,]),1:3]
write.table(specific_hyperEAC, file="FigureS6/specific_hyperEAC_Tumor.bed", quote=F, sep="\t", row.names = F, col.names = F)
write.table(specific_hyperESCC, file="FigureS6/specific_hyperESCC_Tumor.bed", quote=F, sep="\t", row.names = F, col.names = F)

#####FigureS6C
#./Figure2_annot_motifs_pathway/runAnnot.sh FigureS6/specific_hyperEAC_Tumor.bed
#./Figure2_annot_motifs_pathway/runAnnot.sh FigureS6/specific_hyperESCC_Tumor.bed

annotFileList=dir("FigureS6/", pattern=".annot.gene.txt", full=T)
for(file in annotFileList){
  annotData=read_tsv(file)
  fileName=gsub(".annot.gene.txt", "", basename(file))
  colnames(annotData)[1]="PeakID"
  annotData=as.data.frame(annotData[,colnames(annotData)%in%c("PeakID", "Annotation")])
  annotData$Annotation=gsub(" \\(.*", "", annotData$Annotation)
  annotData[annotData$Annotation%in%c("5' UTR","3' UTR", "TTS"),]$Annotation="exon"
  
  annotData$Annotation=factor(annotData$Annotation, levels=c("promoter-TSS", "exon", "intron", "Intergenic", "non-coding"))
  plotdata=as.data.frame(table(annotData$Annotation))
  colnames(plotdata)=c("Type", "Number")
  plotdata$Type=factor(plotdata$Type, levels=as.character(plotdata$Type))
  plotdata$Ratio=plotdata$Number/sum(plotdata$Number)
  plotdata=plotdata[,c(1,3)]
  colnames(plotdata)=c("Type", fileName)
  if(file==annotFileList[1]){
    result=plotdata
  }else{
    result=merge(result, plotdata, by.x="Type", by.y="Type", all=T)
  }
}
result[is.na(result)]=0
colnames(result)=c("Type","tshyperEAC", "tshyperESCC")
result=melt(result, id.vars = "Type")
colnames(result)=c("Type", "Sample", "Ratio")
p<-ggplot(data=result, aes(x=Sample, y=Ratio, fill=Type)) + geom_bar(stat="identity", width = 0.7)
p=p+ylab("Ratio")+xlab("")
p=p+theme_classic()
pdf("FigureS6/FigureS6C.pdf", width=4, height=4)
print(p)
dev.off()
write.table(result, file="FigureS6/FigureS6C_data.txt", row.names = F, col.names = T, sep="\t", quote=F)

#####FigureS6DE
##specific_hypo_EAC DMRs ("FigureS6/specific_hyperEAC_Tumor.bed" and EAC-downregulated genes in "meta/TCGA_EACvsNormal.DESeq2.txt")
##specific_hypo_ESCC DMRs ("FigureS6/specific_hyperESCC_Tumor.bed" and ESCC downregulated genes in "meta/ESCCTumorvsNorml.DESeq2.txt")
plotCistomePathway=function(cistromePathwayFile, FDRcutoff, saveFile){
  cistromePathwayResult=read_tsv(cistromePathwayFile)
  colnames(cistromePathwayResult)=c("GO_term", "detail", "enrichment", "pValue", "FDR", "Genes")
  cistromePathwayResult=cistromePathwayResult[order(cistromePathwayResult$FDR),]
  cistromePathwayResult=cistromePathwayResult[cistromePathwayResult$FDR<FDRcutoff,]
  cistromePathwayResult$FDR=-log10(cistromePathwayResult$FDR)
  if(nrow(cistromePathwayResult)>15){
    cistromePathwayResult=cistromePathwayResult[1:15,]
  }
  cistromePathwayResult$GO_term=gsub("\\(.*\\)", "", cistromePathwayResult$GO_term)
  cistromePathwayResult$GO_term=factor(cistromePathwayResult$GO_term, levels=rev(unique(cistromePathwayResult$GO_term)))
  
  p<-ggplot(data=cistromePathwayResult, aes(x=GO_term, y=FDR)) + geom_bar(stat="identity", fill="black")+ coord_flip()+theme_classic()
  p=p+ylab("-log10(FDR)")+xlab("")
  p=p+theme(axis.text = element_text(colour = "black"))
  pdf(saveFile,width=10, height=2+nrow(cistromePathwayResult)/5)
  print(p)
  dev.off()
  write.table(cistromePathwayResult, file=gsub(".pdf", "_data.txt", saveFile), row.names = F, 
              col.names = T, sep="\t", quote=F)
}
plotCistomePathway("FigureS6/specific_hyperEAC_Tumor_CistromeGO_go_bp_result.txt", 0.05, "FigureS6/FigureS6D.pdf")
plotCistomePathway("FigureS6/specific_hyperESCC_Tumor_CistromeGO_go_bp_result.txt", 0.1, "FigureS6/FigureS6E.pdf")

########################################################FigureS4C#####################################################
##get mean methylation after masking unionPMDs
###./Shell/FigureS4C/getRemovePMDsMeanBetaValues.sh
##merge results into "FigureS4C/Methylation_maskedUnionPMDs.txt"
colorInfo=read.table("FigureS4C/Methylation_maskedUnionPMDs.txt", sep="\t", stringsAsFactors = F, header=T)
plotPCA2=function(plotdata,annotation_col, colorInfo1, saveFile){
  plotdata=t(plotdata)
  pca <- prcomp(plotdata,scale = TRUE)
  xlab <- paste("PC1","(",round((summary(pca))$importance[2,1]*100,1),"%)",sep="")
  ylab <- paste("PC2","(",round((summary(pca))$importance[2,2]*100,1),"%)",sep="")
  x<-"PC1"
  y<-"PC2"
  data_x <- data.frame(varnames=rownames(plotdata), pca$x)
  data_x=data_x[,colnames(data_x)%in%c("varnames", "PC1", "PC2")]
  data_x <- merge(data_x, annotation_col, by.x="varnames", by.y="Sample")
  data_x=merge(data_x, colorInfo1[,c(1,3)], by.x="varnames", by.y="Sample")
  colnames(data_x)[5]=c("maskUnionPMDs")
  my_theme=theme_bw()+ theme(plot.background=element_blank(), panel.background=element_blank(),panel.grid.minor=element_blank(),
                             panel.grid.major=element_blank(),axis.title=element_text(color="black",size=12),
                             axis.text=element_text(size=12), legend.title = element_text(face="bold"))
  p1 <- ggplot(data_x, aes(PC1,PC2, label = varnames,color=Type))+geom_point(size=4, alpha=1)+coord_equal(ratio=1)+xlab(xlab)+ylab(ylab)
  p1=p1+scale_color_manual(values=ann_colors$Type)+my_theme
  p2 <- ggplot(data_x, aes(PC1,PC2, label = varnames,color=maskUnionPMDs))+geom_point(size=4, alpha=1)+coord_equal(ratio=1)+xlab(xlab)+ylab(ylab)
  p2= p2+scale_color_gradientn(name="global methylation", colours = c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"), limits=c(0,1))+my_theme
  
  data_x$Type2="Others"
  data_x[data_x$varnames%in%c("ESCC_Nonmalignant_1", "ESCC_Nonmalignant_2", "ESCC_Nonmalignant_3"),]$Type2="Our_data"
  data_x[data_x$varnames%in%c("ESCC_Nonmalignant_53F", "ESCC_Nonmalignant_54M"),]$Type2="ENCODE_data"
  data_x$Type2=factor(data_x$Type2, levels=c("Others", "Our_data", "ENCODE_data"))
  p3 <- ggplot(data_x, aes(PC1,PC2, color=Type2))+geom_point(size=4, alpha=1)+coord_equal(ratio=1)+xlab(xlab)+ylab(ylab)
  p3= p3+scale_color_manual(values=c("grey", "orange", "#B3672B"))+my_theme
  png(gsub(".pdf", ".png", saveFile), res=300, width = 1500, height=1500)
  print(p3)
  dev.off()
  
  pdf(saveFile,width = 18,height =5)
  print(ggarrange(p1, p2, p3,align = "hv", nrow=1))
  dev.off()
  
  write.table(data_x, gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
}
load("Data/MaskUnionPMDs_DMRs/AllSample_beta_cov7.maskPMD.RData")
colnames(betaMatrix_cov7)=gsub("_Normal", "_Nonmalignant", colnames(betaMatrix_cov7))
plotdata=betaMatrix_cov7
rowVar=rowVars(as.matrix(plotdata[,-1:-3]))
plotdata=plotdata[order(rowVar,decreasing = T),]
plotdata2=plotdata[1:8000,-1:-3]
plotPCA2(plotdata2, annotation_col, colorInfo, "FigureS4C/FigureS4C.pdf")

########################################################Figure5B, 6B, S3F, S4B##############################################
###get background
##Rscript Shell/Figure5B_6B_S4F_S5B/getBackground2.R Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed 10
##Rscript Shell/Figure5B_6B_S4F_S5B/getBackground2.R Data/MaskUnionPMDs_DMRs/hypoESCC_Tumor.bed 10
##Rscript Shell/Figure5B_6B_S4F_S5B/getBackground2.R Data/MaskUnionPMDs_DMRs/specific_hypo_EAC_Tumor.bed 10
##Rscript Shell/Figure5B_6B_S4F_S5B/getBackground2.R Data/MaskUnionPMDs_DMRs/specific_hypo_ESCC_Tumor.bed 10
##Rscript Shell/Figure5B_6B_S4F_S5B/getBackground2.R Data/MaskUnionPMDs_DMRs/nonSpecific_hypo_EAC_Tumor.bed 10
##Rscript Shell/Figure5B_6B_S4F_S5B/getBackground2.R Data/MaskUnionPMDs_DMRs/nonSpecific_hypo_ESCC_Tumor.bed 10

####annot by homer
#./Shell/Figure5B_6B_S4F_S5B/runAnnot.sh Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed
#./Shell/Figure5B_6B_S4F_S5B/runAnnot.sh Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.BackgroundMotif.10X.bed
#./Shell/Figure5B_6B_S4F_S5B/runAnnot.sh Data/MaskUnionPMDs_DMRs/hypoESCC_Tumor.ned
#./Shell/Figure5B_6B_S4F_S5B/runAnnot.sh Data/MaskUnionPMDs_DMRs/hypoESCC_Tumor.BackgroundMotif.10X.bed
#./Shell/Figure5B_6B_S4F_S5B/runAnnot.sh Data/MaskUnionPMDs_DMRs/specific_hypo_EAC_Tumor.bed
#./Shell/Figure5B_6B_S4F_S5B/runAnnot.sh Data/MaskUnionPMDs_DMRs/specific_hypo_EAC_Tumor.BackgroundMotif.10X.bed
#./Shell/Figure5B_6B_S4F_S5B/runAnnot.sh Data/MaskUnionPMDs_DMRs/specific_hypo_ESCC_Tumor.bed
#./Shell/Figure5B_6B_S4F_S5B/runAnnot.sh Data/MaskUnionPMDs_DMRs/specific_hypo_ESCC_Tumor.BackgroundMotif.10X.bed
#./Shell/Figure5B_6B_S4F_S5B/runAnnot.sh Data/MaskUnionPMDs_DMRs/nonSpecific_hypo_EAC_Tumor.bed
#./Shell/Figure5B_6B_S4F_S5B/runAnnot.sh Data/MaskUnionPMDs_DMRs/nonSpecific_hypo_ESCC_Tumor.bed

annotFileList=dir("Data/Figure5B_6B_S4F_S5B/", pattern=".annot.gene.txt", full=T)
for(file in annotFileList){
  annotData=read_tsv(file)
  fileName=gsub(".annot.gene.txt", "", basename(file))
  colnames(annotData)[1]="PeakID"
  annotData=as.data.frame(annotData[,colnames(annotData)%in%c("PeakID", "Annotation")])
  annotData$Annotation=gsub(" \\(.*", "", annotData$Annotation)
  annotData[annotData$Annotation%in%c("5' UTR","3' UTR", "TTS"),]$Annotation="exon"
  
  annotData$Annotation=factor(annotData$Annotation, levels=c("promoter-TSS", "exon", "intron", "Intergenic", "non-coding"))
  plotdata=as.data.frame(table(annotData$Annotation))
  colnames(plotdata)=c("Type", "Number")
  plotdata$Type=factor(plotdata$Type, levels=as.character(plotdata$Type))
  plotdata$Ratio=plotdata$Number/sum(plotdata$Number)
  plotdata=plotdata[,c(1,3)]
  colnames(plotdata)=c("Type", fileName)
  if(file==annotFileList[1]){
    result=plotdata
  }else{
    result=merge(result, plotdata, by.x="Type", by.y="Type", all=T)
  }
}
result[is.na(result)]=0
colnames(result)=c("Type","hypoEAC", "hypoEAC.Back.Nicole", 
                   "hypoESCC", "hypoESCC.Back.Nicole", "ntshypoEAC", "ntshypoESCC",
                   "tshypoEAC", "tshypoEAC.Back.Nicole",  
                   "tshypoESCC", "tshypoESCC.Back.Nicole")
result=melt(result, id.vars = "Type")
colnames(result)=c("Type", "Sample", "Ratio")
result$Group=""
result[grep("^hypoEAC", result$Sample),]$Group="hypoEAC"
result[grep("^hypoESCC", result$Sample),]$Group="hypoESCC"
result[grep("^ntshypoEAC", result$Sample),]$Group="ntshypoEAC"
result[grep("^ntshypoESCC", result$Sample),]$Group="ntshypoESCC"
result[grep("^tshypoEAC", result$Sample),]$Group="tshypoEAC"
result[grep("^tshypoESCC", result$Sample),]$Group="tshypoESCC"
result$SampleGroup="DMRs"
result[grep("Back.Nicole", result$Sample),]$SampleGroup="Nicole_NC"
result$Group=factor(result$Group, levels=c("hypoEAC", "hypoESCC", "tshypoEAC", "tshypoESCC", "ntshypoEAC", "ntshypoESCC"))
result$SampleGroup=factor(result$SampleGroup, levels=c("DMRs", "Nicole_NC"))

###Figure 5B
result1=result[result$Sample%in%c("hypoEAC","hypoESCC"),]
p<-ggplot(data=result1, aes(x=Sample, y=Ratio, fill=Type)) + geom_bar(stat="identity", width = 0.7)
p=p+ylab("Ratio")+xlab("")+scale_fill_manual(values=c("#E41A1C", "#984EA3", "#4DAF4A", "#377EB8", "#E76BF3"))
p=p+theme_classic()
pdf("Figure5B_6B_S4F_S5B/Figure5B_hypoDMR.pdf", height=4, width=4)
print(p)
dev.off()
write.table(result1, "Figure5B_6B_S4F_S5B/Figure5B_hypoDMR.txt", row.names = F, col.names = T, sep="\t", quote=F)

###Figure 6B
result2=result[result$Sample%in%c("tshypoEAC","tshypoESCC"),]
p<-ggplot(data=result2, aes(x=Sample, y=Ratio, fill=Type)) + geom_bar(stat="identity", width = 0.7)
p=p+ylab("Ratio")+xlab("")+scale_fill_manual(values=c("#E41A1C", "#984EA3", "#4DAF4A", "#377EB8", "#E76BF3"))
p=p+theme_classic()
pdf("Figure5B_6B_S4F_S5B/Figure6B_tshypoDMR.pdf", height=4, width=4)
print(p)
dev.off()
write.table(result1, "Figure5B_6B_S4F_S5B/Figure6B_tshypoDMR.txt", row.names = F, col.names = T, sep="\t", quote=F)

###FigureS4F
result3=result[result$Sample%in%c("hypoEAC.Back.Nicole","hypoESCC.Back.Nicole"),]
p<-ggplot(data=result3, aes(x=Sample, y=Ratio, fill=Type)) + geom_bar(stat="identity", width = 0.7)
p=p+ylab("Ratio")+xlab("")+scale_fill_manual(values=c("#E41A1C", "#984EA3", "#4DAF4A", "#377EB8", "#E76BF3"))
p=p+theme_classic()
pdf("Figure5B_6B_S4F_S5B/FigureS4F_hypoDMR_background.pdf", height=4, width=4)
print(p)
dev.off()
write.table(result1, "Figure5B_6B_S4F_S5B/FigureS4F_hypoDMR_background.txt", row.names = F, col.names = T, sep="\t", quote=F)

###FigureS5B
result4=result[result$Sample%in%c("tshypoEAC.Back.Nicole","tshypoESCC.Back.Nicole"),]
p<-ggplot(data=result4, aes(x=Sample, y=Ratio, fill=Type)) + geom_bar(stat="identity", width = 0.7)
p=p+ylab("Ratio")+xlab("")+scale_fill_manual(values=c("#E41A1C", "#984EA3", "#4DAF4A", "#377EB8", "#E76BF3"))
p=p+theme_classic()
pdf("Figure5B_6B_S4F_S5B/FigureS6B_tshypoDMR_background.pdf", height=4, width=4)
print(p)
dev.off()
write.table(result4, "Figure5B_6B_S4F_S5B/FigureS5B_tshypoDMR_background.txt", row.names = F, col.names = T, sep="\t", quote=F)

###Figure5B length####
getDMRLength=function(file){
  peakData=read.table(file)
  peakData$width=peakData$V3-peakData$V2
  peakData[peakData$width>10000,]$width=10000
  return(peakData)
}
hypoEAC_length=getDMRLength("Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed")
hypoEAC_length$Sample="hypoEAC_DMRs"
hypoESCC_length=getDMRLength("Data/MaskUnionPMDs_DMRs/hypoESCC_Tumor.bed")
hypoESCC_length$Sample="hypoESCC_DMRs"
plotdata=rbind(hypoEAC_length[,c(5,4)], hypoESCC_length[,c(5,4)])
plotdata$Sample=factor(plotdata$Sample, levels=c("hypoESCC_DMRs", "hypoEAC_DMRs"))
p=ggplot(plotdata, aes(x=width, color=Sample)) + geom_density(size=1)
p=p+scale_color_manual(values=c("#EA3323", "#0000F5"))
p = p+ theme_classic() +  ggtitle(paste("DMR length distribution" ,sep="")) + xlab("Base pair")
p = p + theme( plot.title = element_text(size = 16, face = "bold",hjust = 0.5),
               axis.title = element_text(face="bold",size=14),
               axis.text.x=element_text(colour="black", size = 12),
               axis.text.y=element_text(colour="black", size = 12),
               axis.line = element_line(size=0.5, colour = "black"),
               legend.title = element_text(face="bold",size=14),
               legend.text = element_text(colour="black", size = 12))
pdf("Figure5B_6B_S4F_S5B/Figure5B_length.pdf",height=4, width=5)
print(p)
dev.off()
write.table(plotdata, file="Figure5B_6B_S4F_S5B/Figure5B_length.txt", row.names = F, col.names = T, sep="\t", quote=F)

########################################################Figure5CD and S4GH##############################################
#./Shell/Figure5CD_S4GH/run.plotheatmapTCGA_ATAC.sh tumor_hypoDMRs Data/Figure5CD_S4GH/
#./Shell/Figure5CD_S4GH/run.plotheatmapTCGA_ATAC.sh tumor_hypoDMRs_background10X Data/Figure5CD_S4GH/

#./Shell/Figure5CD_S4GH/run.plotheatmapESCACellsH3K27ac.sh tumor_hypoDMRs Data/Figure5CD_S4GH/
#./Shell/Figure5CD_S4GH/run.plotheatmapESCACellsH3K27ac.sh tumor_hypoDMRs_background10X Data/Figure5CD_S4GH/

plotProfile2=function(fileName, title, ysize, saveFile){
  deeptoolsHeatmapHeatmap=read.table(fileName, nrows=1, comment.char ="",stringsAsFactors = F)
  num=c(as.numeric(strsplit(deeptoolsHeatmapHeatmap[1,1],":")[[1]][2]), as.numeric(strsplit(deeptoolsHeatmapHeatmap[1,2],":")[[1]][2]))
  startList=c(1, num[1]+1)
  endList=c(num[1], num[1]+num[2])
  deeptoolsHeatmapHeatmapHeader=read.table(fileName, skip =2, nrows=1, comment.char ="",stringsAsFactors = F)
  deeptoolsHeatmapHeatmapHeader=deeptoolsHeatmapHeatmapHeader[,-1:-(length(num))]
  deeptoolsHeatmapHeatmapHeader=data.frame(type=t(deeptoolsHeatmapHeatmapHeader),stringsAsFactors = F)
  sampleNames=unique(deeptoolsHeatmapHeatmapHeader$type)
  deeptoolsHeatmapHeatmap=read.table(fileName, skip =3, comment.char ="",stringsAsFactors = F)
  rowNum=ncol(deeptoolsHeatmapHeatmap)/length(sampleNames)
  
  if(title%in%"ATACseq"){
    sampleLength=6
  }else if(title%in%"Cell_line"){
    sampleLength=7
  }else if(title%in%"EAC_H3K27ac"){
    sampleLength=10
  }else if(title%in%"ESCA_tumor_H3K27ac"){
    sampleLength=10
  }else if(title%in%"ESCC_tumor_H3K27ac"){
    sampleLength=10
  }
  par(mfrow=c(1, 2))
  plotList=list()
  for(j in 1:length(num)){
    plotResult=as.data.frame(matrix(numeric(0),nrow=rowNum))
    for(i in 1:length(sampleNames)){
      print(i)
      colStart=rowNum*(i-1)+1
      colEnd=rowNum*i
      tempData=deeptoolsHeatmapHeatmap[startList[j]:endList[j], colStart:colEnd]
      tempData=data.frame(colMeans(tempData))
      rownames(tempData)=paste("V", 1:rowNum, sep="")
      plotResult=cbind(plotResult, tempData)
    }
    colorList=c("#0000F5", "#EA3323")
    legendList=c("EAC Tumor", "ESCC Tumor")
    
    par(mar = c(3, 5.5, 2, 2))
    if(length(grep("ATAC",fileName))>0){
      plot(plotResult[,1], type="l", col = colorList[1], lwd = 3,axes=F,  xlab="", ylab="Signal Confidence Scores", xaxt = "n", ylim=c(0,ysize),cex.lab=1.4, cex.axis=1.4)
    }else{
      plot(plotResult[,1], type="l", col = colorList[1], lwd = 3,axes=F,  xlab="", ylab="CPM", xaxt = "n", ylim=c(0,ysize),cex.lab=1.4, cex.axis=1.4)
    }
    for(i in 2:sampleLength){
      lines(plotResult[,i], type="l", col = colorList[1], lwd = 3)
    }
    if(sampleLength<ncol(plotResult)){
      for(k in (sampleLength+1):ncol(plotResult)){
        lines(plotResult[,k], type="l", col = colorList[2], lwd = 3)
      }
    }
    axis(1, c(0,100), c("Start","end"), cex.axis=1.3)
    axis(2, cex.axis=1.5)
    
    plotResult=data.frame(t(plotResult))
    plotResult$type=c(rep(legendList[1], sampleLength), rep(legendList[2], nrow(plotResult)-sampleLength))
    
    if(j==1){
      title("EAC/GEJ tumor hypoDMRs",cex.main=1.4)
      legend(0,ysize,legend=legendList,pch=16, col=colorList, bty="n", lwd = 4,cex =1)
      write.table(plotResult, gsub(".pdf", ".EAC.txt", saveFile), sep="\t", row.names = F, col.names = T, quote=F)
    }else{
      title("ESCC tumor hypoDMRs",cex.main=1.4)
      write.table(plotResult, gsub(".pdf", ".ESCC.txt", saveFile), sep="\t", row.names = F, col.names = T, quote=F)
    }
  }
  p=recordPlot()
  pdf(saveFile,width = 8, height=3.5)
  print(p)
  dev.off()
}
plotProfile2("Data/Figure5CD_S4GH/tumor_hypoDMRs.TCGA_ATAC.pValue.tab",  "ATACseq", 22, "Figure5CD_S4GH/Figure5CD.DMR_ATAC.pdf")
plotProfile2("Data/Figure5CD_S4GH/tumor_hypoDMRs_background10X.TCGA_ATAC.pValue.tab", "ATACseq", 22, "Figure5CD_S4GH/Figure5CD.DMR_background_ATAC.pdf")
plotProfile2("Data/Figure5CD_S4GH/tumor_hypoDMRs.ESCACellsH3K27ac.SubtractCPM.tab", "Cell_line", 0.12, "Figure5CD_S4GH/FigureS4GH_cellH3K27ac.pdf")
plotProfile2("Data/Figure5CD_S4GH/tumor_hypoDMRs_background10X.ESCACellsH3K27ac.SubtractCPM.tab", "Cell_line", 0.12, "Figure5CD_S4GH/FigureS4GH_cellH3K27ac_background.pdf")

########################################################Figure5EF and Figure6CD##############################################
##using cistrom-GO website
##hypoEAC DMRs("Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed" and EAC up-regulated genes in "meta/TCGA_ESCCvsEAC.DESeq2.txt")
##hypoESCC DMRs("Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed" and ESCC up-regulated genes in "meta/TCGA_ESCCvsEAC.DESeq2.txt")
##specific_hypo_EAC DMRs ("Data/MaskUnionPMDs_DMRs/specific_hypo_EAC_Tumor.bed" and EAC-upregulated genes in "meta/TCGA_EACvsNormal.DESeq2.txt")
##specific_hypo_ESCC DMRs ("Data/MaskUnionPMDs_DMRs/specific_hypo_ESCC_Tumor.bed" and ESCC up-regulated genes in "meta/ESCCTumorvsNorml.DESeq2.txt")

plotCistomePathway=function(cistromePathwayFile, saveFile){
  cistromePathwayResult=read_tsv(cistromePathwayFile)
  colnames(cistromePathwayResult)=c("GO_term", "detail", "enrichment", "pValue", "FDR", "Genes")
  cistromePathwayResult=cistromePathwayResult[order(cistromePathwayResult$FDR),]
  cistromePathwayResult=cistromePathwayResult[cistromePathwayResult$FDR<0.05,]
  cistromePathwayResult$FDR=-log10(cistromePathwayResult$FDR)
  if(nrow(cistromePathwayResult)>15){
    cistromePathwayResult=cistromePathwayResult[1:15,]
  }
  cistromePathwayResult$GO_term=gsub("\\(.*\\)", "", cistromePathwayResult$GO_term)
  cistromePathwayResult$GO_term=factor(cistromePathwayResult$GO_term, levels=rev(unique(cistromePathwayResult$GO_term)))
  
  p<-ggplot(data=cistromePathwayResult, aes(x=GO_term, y=FDR)) + geom_bar(stat="identity", fill="black")+ coord_flip()+theme_classic()
  p=p+ylab("-log10(FDR)")+xlab("")
  p=p+theme(axis.text = element_text(colour = "black"))
  pdf(saveFile,width=10, height=2+nrow(cistromePathwayResult)/5)
  print(p)
  dev.off()
}
plotCistomePathway("Figure5EF_6CD/hypoEACTumor_CistromeGO_go_bp_result.txt", "Figure5EF_6CD/Figure5E.pdf")
plotCistomePathway("Figure5EF_6CD/hypoESCCTumor_CistromeGO_go_bp_result.txt", "Figure5EF_6CD/Figure5F.pdf")
plotCistomePathway("Figure5EF_6CD/tshypoEACTumor_CistromeGO_go_bp_result.txt", "Figure5EF_6CD/Figure6C.pdf")
plotCistomePathway("Figure5EF_6CD/tshypoESCCTumor_CistromeGO_go_bp_result.txt", "Figure5EF_6CD/Figure6D.pdf")


########################################################Figure5GH##############################################
###motif analysis####
runSplitBackground=function(fileIndex, inputPath, outputPath){
  load(paste(inputPath, "/", fileIndex, ".RData",sep=""))
  data=data.frame(chr=background.set@seqnames, start=background.set@ranges@start, end=(background.set@ranges@start+background.set@ranges@width-1),stringsAsFactors =F)
  num=ceiling(nrow(data)/7000)
  if(num>1){
    for(i in 1:(num-1)){
      temp=data[(7000*(i-1)+1):(7000*i),]
      write.table(temp, file=paste(outputPath, "/", fileIndex, "_", i,".bed",sep=""), sep="\t", row.names=F, col.names=F, quote=F)
    }
    temp=data[(7000*(num-1)+1):nrow(data),]
    write.table(temp, file=paste(outputPath, "/", fileIndex, "_", num,".bed",sep=""), sep="\t", row.names=F, col.names=F, quote=F)
  }else{
    write.table(data, file=paste(outputPath, "/", fileIndex, "_1.bed",sep=""), sep="\t", row.names=F, col.names=F, quote=F)
  }
}
runSplitBackground("hypoEAC_Tumor.10L.BackgroundMotif", "Data/MaskUnionPMDs_DMRs/ESCCTumorvsEACTumor/DMRPlot_20201121_500_0.1", "Data/Figure5GH/")
runSplitBackground("hypoESCC_Tumor.10L.BackgroundMotif", "MaskUnionPMDs_DMRs/ESCCTumorvsEACTumor/DMRPlot_20201121_500_0.1", "Data/Figure5GH/")
runSplitBackground("specific_hypo_EAC_Tumor.10L.BackgroundMotif", "MaskUnionPMDs_DMRs/ESCCTumorvsEACTumor/DMRPlot_20201121_500_0.1", "Data/Figure5GH/")
runSplitBackground("specific_hypo_ESCC_Tumor.10L.BackgroundMotif", "MaskUnionPMDs_DMRs/ESCCTumorvsEACTumor/DMRPlot_20201121_500_0.1", "Data/Figure5GH/")

###bash Shell/Figure5GH/runMotifAnnotBatch.sh
library(ELMER)
runELMER=function(foregroundFile, backgroundFile, outputFile){
  foreground <- ELMER:::getMatrix(filename = foregroundFile)
  background <- ELMER:::getMatrix(filename = backgroundFile)
  or <- calculateEnrichement(foreground,background)
  save(or, file=paste(outputFile, "_or.RData",sep=""))
}
runELMER("Data/Figure5GH/hypoEAC_Tumor.annot.txt", "Data/Figure5GH//hypoEAC_Tumor.BackgroundMotif.10X.annot.txt" ,"Data/Figure5GH//hypoEAC_Tumor")
runELMER("Data/Figure5GH/hypoESCC_Tumor.annot.txt", "Data/Figure5GH//hypoESCC_Tumor.BackgroundMotif.10X.annot.txt" ,"Data/Figure5GH//hypoESCC_Tumor")
runELMER("Data/Figure5GH/specific_hypo_EAC_Tumor.annot.txt", "Data/Figure5GH//specific_hypo_EAC_Tumor.BackgroundMotif.annot.txt" ,"Data/Figure5GH//specific_hypo_EAC_Tumor")
runELMER("Data/Figure5GH/specific_hypo_ESCC_Tumor.annot.txt", "Data/Figure5GH//specific_hypo_ESCC_Tumor.BackgroundMotif.10X.annot.txt" ,"Data/Figure5GH//specific_hypo_ESCC_Tumor")
runELMER("Data/Figure5GH/nonSpecific_hypo_EAC_Tumor.annot.txt", "Data/Figure5GH//nonSpecific_hypo_EAC_Tumor.BackgroundMotif.10X.annot.txt" ,"Data/Figure5GH//nonSpecific_hypo_EAC_Tumor")
runELMER("Data/Figure5GH/nonSpecific_hypo_ESCC_Tumor.annot.txt", "Data/Figure5GH//nonSpecific_hypo_ESCC_Tumor.BackgroundMotif.10X.annot.txt" ,"Data/Figure5GH//nonSpecific_hypo_ESCC_Tumor")

###plot figures
tfs=read.table("meta/HOCOMOCOv11_core_annotation_HUMAN_mono.tsv", sep="\t", header=T, stringsAsFactors = F)
tfs=tfs[,c(1:2,14,15)]
colnames(tfs)=c("motif", "gene", "family","subfamily")
tfs$family=gsub("\\{.*", "", tfs$family)
tfs$subfamily=gsub("\\{.*", "", tfs$subfamily)
geneInformation=read.table("meta/gencode.v22.all.20190716.txt", sep="\t", header=T, stringsAsFactors = F)
geneInformation=unique(geneInformation[,colnames(geneInformation)%in%c("GeneName","GeneId","Genebiotype")])
geneInformation$GeneId=gsub("\\.[0-9]*", "", geneInformation$GeneId)
load("meta/TCGA-ESCA.clinc.RData")
load("meta/TCGA-ESCA_Expression_hg38.v2.RData")
EACTumor_samples=sampleInformation[sampleInformation$Label%in%"EAC_Tumor",]$SampleName
ESCCTumor_samples=sampleInformation[sampleInformation$Label%in%"ESCC_Tumor",]$SampleName
EACNormal_samples=sampleInformation[sampleInformation$Label%in%"EAC_Normal",]$SampleName

ESCAMeanExpression=data.frame(ESCCTumorMean=rowMeans(expression[,colnames(expression)%in%ESCCTumor_samples]),
                              EACTumorMean=rowMeans(expression[,colnames(expression)%in%EACTumor_samples]),
                              EACNormalMean=rowMeans(expression[,colnames(expression)%in%EACNormal_samples]))
ESCAMeanExpression$GeneId=rownames(ESCAMeanExpression)
ESCAMeanExpression=merge(ESCAMeanExpression,geneInformation,by.x="GeneId", by.y="GeneId")

plotMotif=function(ESCAMeanExpression, fileIndex, saveFile){
  load(paste(fileIndex, "_or.RData",sep=""))
  or=or[or$FDR<0.05,]
  or=merge(or, tfs, by.x="motif", by.y="motif")
  or=merge(or, ESCAMeanExpression, by.x="gene", by.y="GeneName")
  if(length(grep("ESCC_Tumor", fileIndex))>0){
    or=or[or$ESCCTumorMean>=5,]
  }else if(length(grep("EAC_Tumor", fileIndex))>0){
    or=or[or$EACTumorMean>=5,]
  }
  or=or[order(or$OR, decreasing = T),]
  motif.enrichment=or[1:30,]
  motif.enrichment$gene <- factor(motif.enrichment$gene,
                                  levels = as.character(motif.enrichment$gene))
  limits <- aes_string(ymax = "upperOR", ymin = "lowerOR")
  
  P <- ggplot(motif.enrichment, aes_string(x = "gene", y = "OR")) + geom_point() +
    geom_errorbar(limits, width = 0.3) + geom_abline(intercept = 1, slope = 0, linetype = "3313") + theme_bw() + theme(panel.grid.major = element_blank()) +
    xlab("Motifs") + ylab("Odds Ratio")+theme(axis.text.x=element_text(face="bold", angle = 90, hjust = 1, vjust = 1),axis.text.y=element_text(face="bold"), axis.title = element_text(face="bold"))+ylim(1, ceiling(max(or$upperOR)))
  
  pdf(saveFile, width=8,heigh=3.5)
  print(P)
  dev.off()
}
plotMotif(ESCAMeanExpression, "Data/Figure5GH/hypoESCC_Tumor", "Figure5GH/Figure5G_EAChypoDMRs.pdf")
plotMotif(ESCAMeanExpression, "Data/Figure5GH/hypoEAC_Tumor", "Figure5GH/Figure5H_ESCChypoDMRs.pdf")

########################################################Figure5IJ##############################################
#./Shell/Figure5IJ/runGATA4.sh
#./Shell/Figure5IJ/runTP63.sh
#./Shell/Figure5IJ/runAnnot.sh Eso26_Gata4flagvsInput_summits.chr.bed GATA4.motif GATA4
#./Shell/Figure5IJ/runAnnot.sh TE5_TP63vsInput_summits.chr.bed P63.motif P63

TFValidation=function(fileName, TF, outputPath, motifAnnotFile, tfPeakFile){
  deeptoolsHeatmapHeatmap=read.table(fileName, nrows=1, comment.char ="",stringsAsFactors = F)
  num=as.numeric(strsplit(deeptoolsHeatmapHeatmap[1,1],":")[[1]][2])
  deeptoolsHeatmapHeatmapHeader=read.table(fileName, skip =2, nrows=1, comment.char ="",stringsAsFactors = F)
  deeptoolsHeatmapHeatmapHeader=deeptoolsHeatmapHeatmapHeader[,-1]
  deeptoolsHeatmapHeatmapHeader=data.frame(type=t(deeptoolsHeatmapHeatmapHeader),stringsAsFactors = F)
  sampleNames=unique(deeptoolsHeatmapHeatmapHeader$type)
  deeptoolsHeatmapHeatmap=read.table(fileName, skip =3, comment.char ="",stringsAsFactors = F)
  test1=deeptoolsHeatmapHeatmap[,1:100]
  breakList=seq(0,2, by=0.01)
  png(paste0(outputPath, TF, "_signal.png"),res=300, width=500, height=1000)
  pheatmap(test1, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, breaks = breakList, 
           color = colorRampPalette(c("#FFF5F0", "red"))(n = length(breakList)))
  dev.off()
  write.table(test1, paste0(outputPath, TF, "_signal.txt"), row.names = F, col.names = F, sep="\t", quote=F)
  test1=colMeans(test1, na.rm=T)
  
  test2=deeptoolsHeatmapHeatmap[,101:200]
  breakList=seq(0,2, by=0.01)
  png(paste0(outputPath, "_H3K27ac_signal.png"),res=300, width=500, height=1000)
  pheatmap(test2, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, breaks = breakList,
           color = colorRampPalette(c("#FFF5F0", "red"))(n = length(breakList)))
  dev.off()
  write.table(test2, paste0(outputPath, "_H3K27ac_signal.txt"), row.names = F, col.names = F, sep="\t", quote=F)
  test2=colMeans(test2, na.rm=T)
  
  test3=deeptoolsHeatmapHeatmap[,201:300]
  breakList=seq(0,1, by=0.01)
  png(paste0(outputPath, "_WGBS.png"),res=300, width=500, height=1000)
  pheatmap(test3, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F, breaks = breakList,
           color = colorRampPalette(c("blue", "red"))(n = length(breakList)), na_col="white")
  dev.off()
  write.table(test3, paste0(outputPath, "_WGBS.txt"), row.names = F, col.names = F, sep="\t", quote=F)
  test3=colMeans(test3, na.rm=T)
  
  ###motif
  sortBedFile=gsub(".tab", ".sort.bed", fileName)
  data=read_tsv(motifAnnotFile)
  data=data[,c(2:4,ncol(data))]
  matrixData=matrix(rep(0, 3000*nrow(data)),ncol=3000, nrow=nrow(data))
  rownames(matrixData)=paste0(data$Chr,":",data$Start+1499,"-",data$End-1500)
  for(i in 1:nrow(data)){
    motifs=as.character(data[i,4])
    if(!is.na(motifs)){
      motifsList=strsplit(motifs, "),")[[1]]
      for(motifInfo in motifsList){
        motifInfoStart=as.numeric(strsplit(motifInfo, "\\(")[[1]][1])
        motifInfoSize=str_count(strsplit(strsplit(motifInfo, "\\(")[[1]][2],",")[[1]][1])
        motifInfoStrand=strsplit(strsplit(motifInfo, "\\(")[[1]][2],",")[[1]][2]
        if((motifInfoStart+motifInfoSize)>3000){
          print(motifInfo)
          matrixData[i, motifInfoStart:3000]=1
        }else{
          matrixData[i, motifInfoStart:(motifInfoStart+motifInfoSize)]=1
        }
      }
    }
  }
  
  ###review
  tmp=matrixData[,1401:1600]
  tmp=rowSums(tmp)
  tmpChoose=tmp[tmp>10]
  tmpChoose=data.frame(number=tmpChoose)
  tmpChoose$chrom=gsub(":.*", "", rownames(tmpChoose))
  tmpChoose$start=as.numeric(gsub("-.*", "", gsub(".*:", "", rownames(tmpChoose))))
  tmpChoose$end=as.numeric(gsub(".*-", "", rownames(tmpChoose)))
  tmpChoose$start=tmpChoose$start-100
  tmpChoose$end=tmpChoose$end+100
  
  matrixData2=matrix(rep(0, (3000/30)*nrow(matrixData)),ncol=100, nrow=nrow(matrixData))
  for(i in 1:nrow(matrixData)){
    for(j in 1:100){
      if(sum(matrixData[i,(j*30-30+1):(j*30)]==1)>=0.1){
        matrixData2[i,j]=1
      }
    }
  }
  rownames(matrixData2)=rownames(matrixData)
  sortBed=read.table(sortBedFile, sep="\t", stringsAsFactors = F)
  rowOrder=factor(rownames(matrixData2), levels=sortBed$V4)
  matrixData2=matrixData2[order(rowOrder),]
  png(paste0(outputPath, TF, "_motif.png"),res=300, width=500, height=1000)
  print(pheatmap(matrixData2, cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F,
           color = c("#FFF5F0", "red"), na_col="white"))
  dev.off()
  write.table(matrixData2, paste0(outputPath, TF, "_motif.txt"), row.names = F, col.names = F, sep="\t", quote=F)
  test4=colMeans(matrixData2)
  
  gsea.layout <- layout(matrix(1:4,nrow = 1), heights = c(2,2,2,2))
  layout.show(gsea.layout)
  test=data.frame(A=unlist(test1), B=unlist(test2), C=unlist(test3), D=unlist(test4), stringsAsFactors = F)
  colnames(test)=c(TF, "H3K27ac", "WGBS", "TF motif")
  for(i in 1:ncol(test)){
    par(mar = c(3, 4, 2, 0))
    if(colnames(test)[i]=="WGBS"){
      name="Beta value"
      ylim=c(0,1)
    }else if(colnames(test)[i]=="TF motif"){
      name="Density"
      ylim=c(0,ceiling(max(test[,i])*10)/10)
    }else{
      name="CPM"
      ylim=c(0,ceiling(max(test[,i])*10)/10)
    }
    plot(test[,i], type="l", col = "darkblue", lwd = 3,axes=F,  xlab="", ylab=name, xaxt = "n", 
         ylim=ylim,cex.lab=1.4, cex.axis=1.4)
    axis(1, c(0,50,100), c("-1.5 kb","center","1.5 kb"), cex.axis=1.3)
    axis(2, cex.axis=1.5)
    title(colnames(test)[i],cex.main=1.4)
  }
  p=recordPlot()
  pdf(paste(outputPath, TF, ".profile.pdf",sep=""),width = 12, height=3)
  print(p)
  dev.off()
  dev.off()

  ######  
  sortBedFile=gsub(".tab", ".sort.bed", fileName)
  sortbedData=read.table(sortBedFile, sep="\t", stringsAsFactors = F)
  sorbedData2=tibble(chrom=sortbedData$V1, start=sortbedData$V2, end=sortbedData$V3, name=sortbedData$V4)
  peakData=read.table(tfPeakFile, sep="\t", stringsAsFactors = F)
  peakData=peakData[,1:4]
  peakData2=tibble(chrom=peakData$V1, start=peakData$V2, end=peakData$V3, name=peakData$V4)
  merge_data=bed_intersect(sorbedData2, peakData2)
  merge_data=merge_data[,c(1,5,6,4)]
  colnames(merge_data)=c("chrom", "start", "end", "name")
  EAC_hypoDMRs=read_bed("Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed",n_fields = 3)
  EAC_hypoDMRs_result=as.data.frame(bed_intersect(merge_data, EAC_hypoDMRs))
  ESCC_hypoDMRs=read_bed("Data/MaskUnionPMDs_DMRs/hypoESCC_Tumor.bed", n_fields = 3)
  ESCC_hypoDMRs_result=as.data.frame(bed_intersect(merge_data, ESCC_hypoDMRs))
  
  motifResult=as.data.frame(bed_intersect(merge_data, tmpChoose[,c("chrom", "start", "end", "number")]))
  result=matrix(rep(0,3*nrow(merge_data)), ncol=3)
  result=as.data.frame(result)
  colnames(result)=c("EAC_hypoDMRs", "ESCC_hypoDMRs", "motif")
  rownames(result)=sortbedData$V4
  result[rownames(result)%in%unique(EAC_hypoDMRs_result$name.x),]$EAC_hypoDMRs=1
  result[rownames(result)%in%unique(ESCC_hypoDMRs_result$name.x),]$ESCC_hypoDMRs=1
  result[rownames(result)%in%unique(motifResult$name.x),]$motif=1
  result2=data.frame(region=rownames(result), result)
  write.table(result2, paste(outputPath, TF, ".PerPeakOverlapDMR.txt",sep=""), row.names = F, col.names = T, sep="\t", quote=F)
  
  plotBarplotData=data.frame(Type=c("EAC_hypoDMRs","EAC_hypoDMRs", "ESCC_hypoDMRs", "ESCC_hypoDMRs"),
                             Motif=c("with", "without", "with", "without"),
                             Ratio=c(nrow(result[result$EAC_hypoDMRs==1&result$motif==1,])/nrow(result),
                                     nrow(result[result$EAC_hypoDMRs==1&result$motif==0,])/nrow(result),
                                     nrow(result[result$ESCC_hypoDMRs==1&result$motif==1,])/nrow(result),
                                     nrow(result[result$ESCC_hypoDMRs==1&result$motif==0,])/nrow(result)),
                             number=c(nrow(result[result$EAC_hypoDMRs==1&result$motif==1,]),
                                      nrow(result[result$EAC_hypoDMRs==1&result$motif==0,]),
                                      nrow(result[result$ESCC_hypoDMRs==1&result$motif==1,]),
                                      nrow(result[result$ESCC_hypoDMRs==1&result$motif==0,])))
  plotBarplotData$Type=factor(plotBarplotData$Type, levels=c("EAC_hypoDMRs","ESCC_hypoDMRs"))
  plotBarplotData$Motif=factor(plotBarplotData$Motif, levels=c("with","without"))
  p=ggplot(data=plotBarplotData, aes(x=Type, y=Ratio, fill=Motif)) + geom_bar(stat="identity", width = 0.7)
  p=p+ylab("Fraction of peaks overlapped with DMRs")+xlab("")
  p=p+scale_fill_manual(values=c("darkred", "darkgrey"))
  p=p+theme_classic()+theme(axis.text.x = element_text(angle=45, size=10, color = "black", hjust = 1, vjust = 1))
  pdf(paste(outputPath, TF, ".PeaksPercontainsDMR.pdf",sep=""), width=3, height=4)
  print(p)
  dev.off()
  
  pie1=ggplot(plotBarplotData[plotBarplotData$Type%in%"EAC_hypoDMRs",], aes(x="", y=number, fill=Motif)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0)+theme_void()+scale_fill_manual(values=c("#8B0000", "#000000"))
  pie2=ggplot(plotBarplotData[plotBarplotData$Type%in%"ESCC_hypoDMRs",], aes(x="", y=number, fill=Motif)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0)+theme_void()+scale_fill_manual(values=c("#8B0000", "#000000"))
  pdf(paste(outputPath, TF, ".PeaksPercontainsDMR.MotifPie.pdf",sep=""), width=8, height=3)
  print(ggarrange(pie1,pie2))
  dev.off()
  
  result[result$EAC_hypoDMRs==1&result$motif==1,]$EAC_hypoDMRs=2
  result[result$ESCC_hypoDMRs==1&result$motif==1,]$ESCC_hypoDMRs=2
  pdf(paste(outputPath, TF, ".PerPeakOverlapDMR.pdf",sep=""), width=2, height=7)
  print(pheatmap(result[,1:2], show_rownames = F, show_colnames = F, cluster_rows = F, cluster_cols = F, color = c("#FAFAFA", "black", "darkred")))
  dev.off()

}
TFValidation("Data/Figure5IJ/ESO26_GATA4_H3K27ac_WGBS.tab", "GATA4", "Figure5IJ/Eso26_", "Data/Figure5IJ/Eso26_Gata4flagvsInput_summits.chr.3kb.GATA4.txt", "Data/Figure5IJ/Eso26_Gata4flagvsInput_peaks.narrowPeak")
TFValidation("Data/Figure5IJ/TE5_TP63_H3K27ac_WGBS.tab", "TP63", "Figure5IJ/P63_", "Data/Figure5IJ/TE5_TP63vsInput_summits.chr.3kb.P63.txt", "Data/Figure5IJ/TE5_TP63vsInput_peaks.narrowPeak")

########################################################Figure6E and FigureS4E##############################################
combinedMotifs=function(nonTsFile, tsFile, type){
  load(paste(nonTsFile, "_or.RData",sep=""))
  # nonTs_or=or[or$FDR<0.05,]
  nonTs_or=or
  nonTs_or=merge(nonTs_or, tfs, by.x="motif", by.y="motif")
  nonTs_or=merge(nonTs_or, ESCAMeanExpression, by.x="gene", by.y="GeneName")
  if(type%in%"ESCCTumor"){
    nonTs_or=nonTs_or[nonTs_or$ESCCTumorMean>=5,]
  }else if(type%in%"EACTumor"){
    nonTs_or=nonTs_or[nonTs_or$EACTumorMean>=5,]
  }
  nonTs_or=nonTs_or[order(nonTs_or$OR, decreasing = T),]
  nonTs_or_top30=nonTs_or[1:30,]
  
  load(paste(tsFile, "_or.RData",sep=""))
  # ts_or=or[or$FDR<0.05,]
  ts_or=or
  ts_or=merge(ts_or, tfs, by.x="motif", by.y="motif")
  ts_or=merge(ts_or, ESCAMeanExpression, by.x="gene", by.y="GeneName")
  if(type%in%"ESCCTumor"){
    ts_or=ts_or[ts_or$ESCCTumorMean>=5,]
  }else if(type%in%"EACTumor"){
    ts_or=ts_or[ts_or$EACTumorMean>=5,]
  }
  ts_or=ts_or[order(ts_or$OR, decreasing = T),]
  ts_or_top30=ts_or[1:30,]
  
  unionTop30=unique(c(nonTs_or_top30$motif, ts_or_top30$motif))
  unionTop30_nonTs=nonTs_or[nonTs_or$motif%in%unionTop30,c(2,1,6,8)]
  unionTop30_ts=ts_or[ts_or$motif%in%unionTop30,c(2,1,6,8)]
  unionTop30_merge=merge(unionTop30_nonTs, unionTop30_ts, by.x=c("motif", "gene"), by.y=c("motif", "gene"),all=T)
  unionTop30_merge=unionTop30_merge[,c(1,2,3,5,4,6)]
  colnames(unionTop30_merge)=c("TF", "gene","non_ts", "ts", "non_ts_FDR","ts_FDR")
  return(unionTop30_merge)
}
EACTumor_mergeTFs=combinedMotifs("Data/Figure5GH/nonSpecific_hypo_EAC_Tumor", "Data/Figure5GH/specific_hypo_EAC_Tumor", "EACTumor")
EACTumor_mergeTFs$TF=gsub("_.*", "", EACTumor_mergeTFs$TF)
ESCCTumor_mergeTFs=combinedMotifs("Data/Figure5GH/nonSpecific_hypo_ESCC_Tumor", "Data/Figure5GH/specific_hypo_ESCC_Tumor", "ESCCTumor")
ESCCTumor_mergeTFs$TF=gsub("_.*", "", ESCCTumor_mergeTFs$TF)

load("meta/TCGA_ESCCvsEAC.DESeq2.RData")
load("meta/TCGA_EACvsNormal.DESeq2.RData")
load("meta/ESCC_Paired_result.RData")
ESCC_EAC_targetResults=ESCC_EAC_result[,c(2,6,7,8)]
ESCC_EAC_targetResults$GeneID=rownames(ESCC_EAC_targetResults)
ESCC_EAC_targetResults=merge(ESCC_EAC_targetResults, geneInformation[,c(2,1)], by.x="GeneID", by.y="GeneId")
EAC_TN_targetResults=EAC_TN_result[,c(2,6,7,8)]
EAC_TN_targetResults$GeneID=rownames(EAC_TN_targetResults)
EAC_TN_targetResults=merge(EAC_TN_targetResults, geneInformation[,c(2,1)], by.x="GeneID", by.y="GeneId")
ESCC_TN_targetResults=ESCC_Paired_resuts[,c(1,26,30,23:24,2)]

###Figure6E
EACTumor_mergeTFs$gene=factor(EACTumor_mergeTFs$gene, levels=as.character(EACTumor_mergeTFs$gene))
ESCCTumor_mergeTFs$gene=factor(ESCCTumor_mergeTFs$gene, levels=as.character(ESCCTumor_mergeTFs$gene))
EACTumor_mergeTFs_result=merge(EACTumor_mergeTFs, ESCC_EAC_targetResults[,c(6,2,4:5)], by.x="gene", by.y="GeneName")
EACTumor_mergeTFs_result$log2FoldChange=-(EACTumor_mergeTFs_result$log2FoldChange)
EACTumor_mergeTFs_result=merge(EACTumor_mergeTFs_result, EAC_TN_targetResults[,c(6,2,4:5)], by.x="gene", by.y="GeneName")
EACTumor_mergeTFs_result=EACTumor_mergeTFs_result[,c(2:6,8,9,12,7,10)]
colnames(EACTumor_mergeTFs_result)=c("TF", "non_tsDMRs", "tsDMRs", "non_tsDMRs_FDR", "tsDMRs_FDR",
                                     "ESCCMean", "EACMean", "EACNonmamlignant", "EACvsESCC","EACTvsEACN")
EACTumor_mergeTFs_result$deltaOR=EACTumor_mergeTFs_result$tsDMRs-EACTumor_mergeTFs_result$non_tsDMRs
EACTumor_mergeTFs_result=EACTumor_mergeTFs_result[order(EACTumor_mergeTFs_result$deltaOR,decreasing = T),]
p=ggplot(EACTumor_mergeTFs_result, aes(x=EACTvsEACN, y=deltaOR, label=TF)) + geom_point(pch=15, size=3)+theme_classic()+geom_text(size=3, angle=15, hjust=1)
p=p+xlab("log2 fold change(EAC tumor v.s. nonmalignant EAC)")+ylab("delta OR (ts hypoDMRs v.s. non-ts hypoDMRs)")
p=p+geom_hline(yintercept=0, linetype="dotted", size=1, color="grey")+geom_vline(xintercept=log2(1.5), linetype="dotted", size=1, color="grey")
pdf("Figure6E_S5E/Figure6E_EACTumor_tsvsnonts.pdf", height=5, width=7)
print(p)
dev.off()
write.table(EACTumor_mergeTFs_result, "Figure6E_S5E/Figure6E_EACTumor_tsvsnonts.txt", row.names = F, col.names = T, sep="\t", quote=F)

###FigureS5E
ESCCTumor_mergeTFs_result=merge(ESCCTumor_mergeTFs, ESCC_EAC_targetResults[,c(6,2,4:5)], by.x="gene", by.y="GeneName")
ESCCTumor_mergeTFs_result=merge(ESCCTumor_mergeTFs_result, ESCC_TN_targetResults[,c(6,2,4,5)], by.x="gene", by.y="GeneName")
ESCCTumor_mergeTFs_result=ESCCTumor_mergeTFs_result[,c(2:6,8:9,12,11,7,10)]
colnames(ESCCTumor_mergeTFs_result)=c("TF", "non_tsDMRs", "tsDMRs", "non_tsDMRs_FDR", "tsDMRs_FDR",
                                      "ESCCMean", "EACMean", "ESCCTMean","ESCCNonmamlignant", "ESCCvsEAC","ESCCTvsESCCN")
ESCCTumor_mergeTFs_result$deltaOR=ESCCTumor_mergeTFs_result$tsDMRs-ESCCTumor_mergeTFs_result$non_tsDMRs
ESCCTumor_mergeTFs_result=ESCCTumor_mergeTFs_result[order(ESCCTumor_mergeTFs_result$deltaOR, decreasing = T),]
p=ggplot(ESCCTumor_mergeTFs_result, aes(x=ESCCTvsESCCN, y=deltaOR, label=TF)) + geom_point(pch=15, size=3)+theme_classic()+geom_text(size=3, angle=15, hjust=1)
p=p+xlab("log2 fold change(ESCC tumor v.s. nonmalignant ES)")+ylab("delta OR (ts hypoDMRs v.s. non-ts hypoDMRs)")
p=p+geom_hline(yintercept=0, linetype="dotted", size=1, color="grey")+geom_vline(xintercept=log2(1.5), linetype="dotted", size=1, color="grey")
pdf("Figure6E_S5E/FigureS5E_ESCCTumor_tsvsnonts.pdf", height=5, width=7)
print(p)
dev.off()
write.table(ESCCTumor_mergeTFs_result, "Figure6E_S5E/FigureS5E_ESCCTumor_tsvsnonts.txt", row.names = F, col.names = T, sep="\t", quote=F)

#########################################FigureS4KL_S6CD########################################################
getDMRMethylationAndExpressionCor=function(annotFile, DEExpMatrix, type, saveFile){
  annotResult=read_tsv(annotFile)
  annotResult=annotResult[,c(2:4,15)]
  targetGenes=unique(annotResult$`Nearest Ensembl`)
  targetGenes=targetGenes[!is.na(targetGenes)]
  targetGenesExp=DEExpMatrix[DEExpMatrix$GeneID%in%targetGenes,]
  if(type%in%"hypoEAC_Tumor"){
    targetGenesExp=targetGenesExp[targetGenesExp$EAC_Mean>1&targetGenesExp$log2FoldChange<(-log2(1.5))&targetGenesExp$padj<0.05,]
  }else if(type%in%"hypoESCC_Tumor"){
    targetGenesExp=targetGenesExp[targetGenesExp$ESCC_Mean>1&targetGenesExp$log2FoldChange>log2(1.5)&targetGenesExp$padj<0.05,]
  }else if(type%in%"tshypoEAC_Tumor"){
    targetGenesExp=targetGenesExp[targetGenesExp$EAC_Mean>1&targetGenesExp$log2FoldChange>(log2(1.5))&targetGenesExp$padj<0.05,]
  }else if(type%in%"tshypoESCC_Tumor"){
    targetGenesExp=targetGenesExp[targetGenesExp$ESCC_Tumor_Mean>1&targetGenesExp$log2FoldChange>(log2(1.5))&targetGenesExp$padj<0.05,]
  }
  print(length(targetGenesExp$GeneID))
  write.table(targetGenesExp$GeneID, file=saveFile, row.names = F, col.names = F, sep="\t", quote=F)
}
getDMRMethylationAndExpressionCor("Data/Figure5GH/hypoEAC_Tumor.annot.txt", ESCC_EAC_targetResults, "hypoEAC_Tumor", "FigureS4KL_S5CD/hypoEAC_Tumor.EACup.gene.txt")
getDMRMethylationAndExpressionCor("Data/Figure5GH/hypoESCC_Tumor.annot.txt", ESCC_EAC_targetResults, "hypoESCC_Tumor", "FigureS4KL_S5CD/hypoESCC_Tumor.ESCCup.gene.txt")
getDMRMethylationAndExpressionCor("Data/Figure5GH/specific_hypo_EAC_Tumor.annot.txt", EAC_TN_targetResults, "tshypoEAC_Tumor", "FigureS4KL_S5CD/specific_hypo_EAC_Tumor.EACup.gene.txt")
getDMRMethylationAndExpressionCor("Data/Figure5GH/specific_hypo_ESCC_Tumor.annot.txt", ESCC_TN_targetResults, "tshypoESCC_Tumor", "FigureS4KL_S5CD/specific_hypo_ESCC_Tumor.ESCCup.gene.txt")

########################################################Figure6F##############################################
getTargetTFDMRsByELMER=function(tsAnnotFile, nontsAnnotFile, TF, saveFile){
  ts_hypoEACDMRs=read_tsv(tsAnnotFile)
  nonts_hypoEACDMRs=read_tsv(nontsAnnotFile)
  targetColnames=colnames(ts_hypoEACDMRs)[grep(TF, colnames(ts_hypoEACDMRs))]
  ts_hypoEACDMRs_TF=ts_hypoEACDMRs[rowSums(!is.na(ts_hypoEACDMRs[,colnames(ts_hypoEACDMRs)%in%targetColnames,drop=F]))>0,]
  colnames(ts_hypoEACDMRs_TF)[1]="PeakID"
  nonts_hypoEACDMRs_TF=nonts_hypoEACDMRs[rowSums(!is.na(nonts_hypoEACDMRs[,colnames(nonts_hypoEACDMRs)%in%targetColnames,drop=F]))>0,]
  colnames(nonts_hypoEACDMRs_TF)[1]="PeakID"
  if(length(grep("nonSpecific", tsAnnotFile))>0){
    plotdata=data.frame(Sample=c("nts hypDMRs", "ts hypoDMRs"), 
                        Ratio=c(nrow(ts_hypoEACDMRs_TF)/nrow(ts_hypoEACDMRs),
                                nrow(nonts_hypoEACDMRs_TF)/nrow(nonts_hypoEACDMRs)))
  }else{
    plotdata=data.frame(Sample=c("ts hypDMRs", "nts hypoDMRs"), 
                        Ratio=c(nrow(ts_hypoEACDMRs_TF)/nrow(ts_hypoEACDMRs),
                                nrow(nonts_hypoEACDMRs_TF)/nrow(nonts_hypoEACDMRs)))
  }
  print(plotdata)
  
  p<-ggplot(data=plotdata, aes(x=Sample, y=Ratio, fill=Sample)) + geom_bar(stat="identity", width = 0.7)
  p=p+ylab("Ratio")+xlab("")+scale_fill_manual(values=c("#BEBEBE", "#A020F0"))+ylim(0, ceiling(max(plotdata$Ratio*10))/10)
  p=p+theme_classic()
  pdf(saveFile, height=4, width=4)
  print(p)
  dev.off()
  num1=nrow(ts_hypoEACDMRs_TF)
  num2=nrow(nonts_hypoEACDMRs_TF)
  num3=nrow(ts_hypoEACDMRs)-num1
  num4=nrow(nonts_hypoEACDMRs)-num2
  if(length(grep("nonSpecific", tsAnnotFile))>0){
    resultData <-matrix(c(num1, num2, num3, num4),nrow = 2,
                        dimnames = list(catalog = c("nts hypoDMRs", "ts hypoDMRs"), Type = c("Target", "NonTarget")))
  }else{
    resultData <-matrix(c(num1, num2, num3, num4),nrow = 2,
                        dimnames = list(catalog = c("ts hypoDMRs", "nts hypoDMRs"), Type = c("Target", "NonTarget")))
  }
  plotdata=rbind(plotdata, data.frame(Sample="pvalue", Ratio=fisher.test(resultData)$p.value))
  write.table(plotdata, gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
}
getTargetTFDMRsByELMER("Data/Figure5GH/specific_hypo_EAC_Tumor.annot.txt", "Data/Figure5GH/nonSpecific_hypo_EAC_Tumor.annot.txt", "HNF4A", "Figure6F/DMR_motif_HNF4A.pdf")

########################################################Figure6G##############################################
##Shell
#wc -l Data/MaskUnionPMDs_DMRs/specific_hypo_EAC_Tumor.bed
##1972
#bedtools intersect -a Data/MaskUnionPMDs_DMRs//specific_hypo_EAC_Tumor.bed -b Data/Figure6G/ESO26_HNF4AvsInput_peaks.narrowPeak -wa | sort | uniq | wc -l
##288
#bedtools intersect -a Data/MaskUnionPMDs_DMRs//specific_hypo_EAC_Tumor.bed -b Data/Figure6G/OE19_HNF4Av2_intersect.bwa.bed -wa | sort | uniq | wc -l
#274

#wc -l Data/MaskUnionPMDs_DMRs//nonSpecific_hypo_EAC_Tumor.bed
##5762
#bedtools intersect -a Data/MaskUnionPMDs_DMRs//nonSpecific_hypo_EAC_Tumor.bed -b Data/Figure6G/ESO26_HNF4AvsInput_peaks.narrowPeak -wa | sort | uniq | wc -l
#418
#bedtools intersect -a Data/MaskUnionPMDs_DMRs//nonSpecific_hypo_EAC_Tumor.bed -b Data/Figure6G/OE19_HNF4Av2_intersect.bwa.bed -wa | sort | uniq | wc -l
#455

plotdata=data.frame(type=c("ts hypoDMRs", "non-ts hypoDMRs","ts hypoDMRs", "non-ts hypoDMRs"), 
                    ratio=c(288/1972, 418/5762, 274/1972, 455/5762),
                    cells=c("ESO26", "ESO26", "OE19", "OE19"))
p<-ggplot(data=plotdata, aes(x=cells, y=ratio, fill=type)) + geom_bar(stat="identity", position=position_dodge(), color="white")+theme_classic()
p=p+xlab("")+ylab("Fraction of DMRs \noverlapping with HNF4A peaks")+scale_fill_manual(values=c("grey", "purple"))
p=p+theme(axis.text =element_text(size=12), axis.title=element_text(size=14))
pdf("Figure6G/Figure6G_HNF4A_peaks_validation.pdf", height=4, width=6)
print(p)
dev.off()

resultData <-matrix(c(288, 418, 1972-288, 5762-418),nrow = 2,
                    dimnames = list(catalog = c("ts hypoDMRs", "non-ts hypoDMRs"), Type = c("Target", "NonTarget")))
p=fisher.test(resultData)$p.value
plotdata=rbind(plotdata, data.frame(type="pvalue", ratio=p, cells="ESO26"))
resultData <-matrix(c(274, 455, 1972-274, 5762-455),nrow = 2,
                    dimnames = list(catalog = c("ts hypoDMRs", "non-ts hypoDMRs"), Type = c("Target", "NonTarget")))
p=fisher.test(resultData)$p.value
plotdata=rbind(plotdata, data.frame(type="pvalue", ratio=p, cells="OE19"))
write.table(plotdata, "Figure6G/Figure6G_HNF4A_peaks_validation.txt", row.names = F, col.names = T, sep="\t", quote=F)

########################################################Figure6H##############################################
###run motif analyis
getMatrix=function(motifsData){
  matrix <- Matrix::Matrix(0, nrow = nrow(motifsData), ncol = ncol(motifsData) - 
                             21, sparse = TRUE)
  colnames(matrix) <- gsub(" Distance From Peak\\(sequence,strand,conservation\\)", 
                           "", colnames(motifsData)[-c(1:21)])
  rownames(matrix) <- motifsData$PeakID
  matrix[!is.na(motifsData[, -c(1:21)])] <- 1
  matrix <- as(matrix, "nsparseMatrix")
  return(matrix)
}
getTargetTFDMRsByELMER=function(foregroundFile, backgroundFile, type, TF, lab, saveFile){
  ts_hypoEACDMRs=read_tsv(foregroundFile)
  nonts_hypoEACDMRs=read_tsv(backgroundFile)
  targetColnames=colnames(ts_hypoEACDMRs)[grep(TF, colnames(ts_hypoEACDMRs))]
  ts_hypoEACDMRs_TF=ts_hypoEACDMRs[rowSums(!is.na(ts_hypoEACDMRs[,colnames(ts_hypoEACDMRs)%in%targetColnames,drop=F]))>0,]
  colnames(ts_hypoEACDMRs_TF)[1]="PeakID"
  nonts_hypoEACDMRs_TF=nonts_hypoEACDMRs[rowSums(!is.na(nonts_hypoEACDMRs[,colnames(nonts_hypoEACDMRs)%in%targetColnames,drop=F]))>0,]
  colnames(nonts_hypoEACDMRs_TF)[1]="PeakID"

  foreground <- getMatrix(ts_hypoEACDMRs_TF)
  background <- getMatrix(nonts_hypoEACDMRs_TF)
  or <- calculateEnrichement(foreground,background)
  or=or[order(or$FDR),]
  or=or[or$FDR<0.05,]
  or=merge(or, tfs, by.x="motif", by.y="motif")
  or=merge(or, ESCAMeanExpression, by.x="gene", by.y="GeneName")
  if(type%in%"EAC"){
    or=or[or$EACTumorMean>=5,]
  }else if(type%in%"ESCC"){
    or=or[or$ESCCTumorMean>=5,]
  }else if(type%in%"EACNormal"){
    or=or[or$EACNormalMean>=5,]
  }
  or=or[order(or$OR, decreasing = T),]
  if(nrow(or)>30){
    motif.enrichment=or[1:30,]
  }else{
    motif.enrichment=or
  }
  motif.enrichment$gene <- factor(motif.enrichment$gene,
                                  levels = as.character(motif.enrichment$gene[nrow(motif.enrichment):1]))
  limits <- aes_string(ymax = "upperOR", ymin = "lowerOR")
  P <- ggplot(motif.enrichment, aes_string(x = "gene", y = "OR")) + geom_point() +
    geom_errorbar(limits, width = 0.3) + coord_flip() + geom_abline(intercept = 1, slope = 0, linetype = "3313") + theme_bw() + theme(panel.grid.major = element_blank()) +
    xlab("Motifs") + ylab(paste0("Odds Ratio (", lab, ")"))+theme(axis.text=element_text(face="bold"),axis.title = element_text(face="bold"))+ylim(0, ceiling(max(or$upperOR)))
  P=P+ggtitle(TF)+theme(plot.title = element_text(hjust = 0.5, size=16))
  pdf(saveFile, width=4, height=4)
  print(P)
  dev.off()
  write.table(motif.enrichment, gsub(".pdf", ".txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
}
getTargetTFDMRsByELMER("Data/Figure5GH/specific_hypo_EAC_Tumor.annot.txt", "Data/Figure5GH/nonSpecific_hypo_EAC_Tumor.annot.txt", 
                       "EAC", "HNF4A", "ts/nts hypoDMR", "Figure6H/Figure6H_HNF4A_cobinding_TF.pdf")


########################################################FigureS4DE and Figure7C#############################################
## bash Shell/FigureS4DE_7C/getDMRMeanBetaValues.sh Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed
## bash Shell/FigureS4DE_7C/getDMRMeanBetaValues.sh Data/MaskUnionPMDs_DMRs/hypoESCC_Tumor.bed
## merge the result into "FigureS4DE_7C/hypoDMR_methylation.txt" 
data=read.table("FigureS4DE_7C/hypoDMR_methylation.txt", sep="\t", stringsAsFactors = F, header=T, check.names = F)
EAC_DMR_WGBS=data[data$Type%in%c("EAC/GEJ_Tumor", "GEJ_Nonmalignant"),]
EAC_DMR_WGBS$Source="WGBS"
colnames(EAC_DMR_WGBS)=c("Sample", "hypoEAC", "hypoESCC", "Type", "Source")
data=melt(data, id.vars = c("Sample", "Type"))
colnames(data)=c("Sample", "SampleType", "RegionType", "Methylation")
plotRegionLineplot2=function(plotdata, type, ylab){
  plotdata$RegionType=factor(plotdata$RegionType, levels=c("EAC/GEJ_tumor_hypoDMRs", "ESCC_tumor_hypoDMRs"))
  p=ggplot(data=plotdata, aes(x=RegionType, y=Methylation, group=Sample)) + geom_line(alpha=0.2, color="darkblue")+ geom_point(alpha=0.5, color="darkblue", size=0.8)+theme_classic()
  p=p+xlab("")+ylab(ylab)+ggtitle(type)+ylim(0, 1)
  p=p+theme(plot.title = element_text(hjust = 0.5, size=14, color="black", face="bold"),
            axis.text.x = element_text(size=11, color="black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size=11, color="black", angle = 0),
            axis.title = element_text(size=13, color="black", face="bold"),
            legend.title =element_text(size=12, color="black", face="bold"),
            legend.text =element_text(size=12, color="black"),legend.position = "none")
  return(p)
}
p1=plotRegionLineplot2(data[data$SampleType%in%"EAC/GEJ_Tumor",], "EAC/GEJ_Tumor", "DMR methylation")
p2=plotRegionLineplot2(data[data$SampleType%in%"ESCC_Tumor",], "ESCC_Tumor", "DMR methylation")
p3=plotRegionLineplot2(data[data$SampleType%in%"GEJ_Nonmalignant",], "GEJ_Nonmalignant", "DMR methylation")
p4=plotRegionLineplot2(data[data$SampleType%in%"ESCC_Nonmalignant",], "ESCC_Nonmalignant", "DMR methylation")
pdf("FigureS4DE_7C/FigureS4DE_hypoDMR_methylation.pdf", width=6, height=4.5)
print(ggarrange(p1,p2, nrow=1))
dev.off()

pdf("FigureS4DE_7C/Figure7C_hypoDMR_methylation.pdf", width=6, height=4.5)
print(ggarrange(p3,p4, nrow=1))
dev.off()

########################################################Figure 7B#############################################
#./Shell/Figure7B/getPerRegionMethylation.sh EAC_specificPMDs.bed
#./Shell/Figure7B/getPerRegionMethylation.sh ESCA_sharedHMDs.bed
#./Shell/Figure7B/getPerRegionMethylation.sh ESCA_sharedPMDs.bed
#./Shell/Figure7B/getPerRegionMethylation.sh ESCC_specificPMDs.bed

calculatePvalue=function(dataMatrix, sampleList1, sampleList2, sampleList1Name, sampleList2Name){
  TtestResult=lapply(1:nrow(dataMatrix), function(x){
    t.test(dataMatrix[x,colnames(dataMatrix)%in%sampleList1],
           dataMatrix[x,colnames(dataMatrix)%in%sampleList2])$p.value
  })
  TtestResult3=do.call(c, TtestResult)
  
  sampleList1Mean=rowMeans(dataMatrix[,colnames(dataMatrix)%in%sampleList1], na.rm=T)
  sampleList2Mean=rowMeans(dataMatrix[,colnames(dataMatrix)%in%sampleList2], na.rm=T)
  
  TtestResult=data.frame(sampleList1Mean= sampleList1Mean, sampleList2Mean= sampleList2Mean, Pvalue=TtestResult3)
  colnames(TtestResult)[1:2]=c(sampleList1Name, sampleList2Name)
  TtestResult$FDR=p.adjust(TtestResult$Pvalue, method="BH")
  rownames(TtestResult)=rownames(dataMatrix)
  return(TtestResult)
}
getDomainMethylation=function(domainPath, pattern){
  for(i in 1:length(annotation_col$Sample)){
    sample=annotation_col$Sample[i]
    domainData=read.table(paste0(domainPath, "/", pattern, ".", sample, ".txt"))
    domainData=domainData[,1:2]
    colnames(domainData)=c("PMD", sample)
    if(i==1){
      domainMeth=domainData
    }else{
      domainMeth=merge(domainMeth, domainData, by.x="PMD", by.y="PMD")
    }
  }
  rownames(domainMeth)=domainMeth$PMD
  domainMeth=domainMeth[,-1]
  
  ESCC_Tumor_sampleList=paste0("ESCC_", c(1:17,19:22))
  EAC_Tumor_sampleList=c(paste0("EAC_", c(1:4,6:7)), paste0("GEJ_", c(1:7)))
  ESCC_Normal_sampleList=paste0("ESCC_Nonmalignant_", c(1:3, "53F", "54M"))
  GEJ_Normal_sampleList=paste0("GEJ_Nonmalignant_", c(1:7))
  targetMethylation_TumorPvalue=calculatePvalue(domainMeth, ESCC_Tumor_sampleList, EAC_Tumor_sampleList, "ESCC_Tumor", "EAC_Tumor")
  targetMethylation_TumorPvalue$deltaMethylation=targetMethylation_TumorPvalue$ESCC_Tumor-targetMethylation_TumorPvalue$EAC_Tumor
  targetGenesMethylation_NormalPvalue=calculatePvalue(domainMeth, ESCC_Normal_sampleList, GEJ_Normal_sampleList, "ESCC_Nonmalignant", "EAC_Nonmalignant")
  targetGenesMethylation_NormalPvalue$deltaMethylation= targetGenesMethylation_NormalPvalue$ESCC_Nonmalignant-targetGenesMethylation_NormalPvalue$EAC_Nonmalignant
  result=list()
  result$tumor=targetMethylation_TumorPvalue
  result$normal=targetGenesMethylation_NormalPvalue
  return(result)
}
EAC_regionLevelPMD=getDomainMethylation("Data/Figure7B/PMD_methylation/", "EAC_specificPMDs")
ESCC_regionLevelPMD=getDomainMethylation("Data/Figure7B/PMD_methylation/", "ESCC_specificPMDs")
ESCA_regionLevelSharedPMD=getDomainMethylation("Data/Figure7B/PMD_methylation/", "ESCA_sharedPMDs")
ESCA_regionLevelSharedHMD=getDomainMethylation("Data/Figure7B/PMD_methylation/", "ESCA_sharedHMDs")

Volcanoplot=function(dataMatrix, title, cutoff, color1, color2, saveFile){
  plotdata=dataMatrix
  xlab=paste0(colnames(plotdata)[1:2], collapse = "-")
  plotdata=plotdata[,colnames(plotdata)%in%c("FDR", "deltaMethylation")]
  plotdata$Type="non-change"
  if(nrow(plotdata[plotdata$FDR<cutoff&plotdata$deltaMethylation>0,])>0){
    plotdata[plotdata$FDR<cutoff&plotdata$deltaMethylation>0,]$Type="Up"
  }
  if(nrow(plotdata[plotdata$FDR<cutoff&plotdata$deltaMethylation<0,])>0){
    plotdata[plotdata$FDR<cutoff&plotdata$deltaMethylation<0,]$Type="Down"
  }
  plotdata$FDR=-log10(plotdata$FDR)
  typeColors=data.frame(type=c("Up", "Down", "non-change"), colors=c(color1, color2, "black"), stringsAsFactors = F)
  plotdata$Type=factor(plotdata$Type, levels=c("Up", "Down", "non-change"))
  
  write.table(plotdata, saveFile, row.names = F, col.names = T, sep="\t", quote=F)
  p <- ggplot(data=plotdata,aes(x=deltaMethylation,y=FDR,colour=Type))+geom_point(pch=15)
  p <- p+labs(x = paste0("Delta methlation (", xlab, ")"), y = "-log10 FDR")+ggtitle(title)+scale_color_manual(values=typeColors[typeColors$type%in%unique(plotdata$Type),]$colors)
  xlim1=min(plotdata$deltaMethylation)+(max(plotdata$deltaMethylation)-min(plotdata$deltaMethylation))/10
  xlim2=max(plotdata$deltaMethylation)-(max(plotdata$deltaMethylation)-min(plotdata$deltaMethylation))/10
  ylim=max(plotdata$FDR)-(max(plotdata$FDR)-min(plotdata$FDR))/10
  if(nrow(plotdata[plotdata$Type%in%"Down",])>0){
    p=p+annotate(geom="text", x=xlim1, y=ylim, label=paste0(round(nrow(plotdata[plotdata$Type%in%"Down",])/nrow(plotdata), 4)*100, "%\n(", comma(nrow(plotdata[plotdata$Type%in%"Down",]),0),"/",comma(nrow(plotdata),0),")"),color="black", size=5)
  }
  if(nrow(plotdata[plotdata$Type%in%"Up",])>0){
    p=p+annotate(geom="text", x=xlim2, y=ylim, label=paste0(round(nrow(plotdata[plotdata$Type%in%"Up",])/nrow(plotdata), 4)*100, "%\n(", comma(nrow(plotdata[plotdata$Type%in%"Up",]),0),"/",comma(nrow(plotdata),0),")"),color="black", size=5)
  }
  p <- p +theme_classic() + theme(axis.text = element_text(colour = "black",size=10), 
                                  axis.title=element_text(colour="black",size=12),
                                  plot.title = element_text(hjust = 0.5,size=16))
  # legend.position = "none")
}
getSignificant=function(dataMatrix, type, targetCutoff, title, color1, color2, saveIndex){
  if(type%in%"EAC"){
    tumorTemp=dataMatrix$tumor[dataMatrix$tumor$FDR<targetCutoff&dataMatrix$tumor$deltaMethylation>0,]
    normalTemp=dataMatrix$normal[rownames(dataMatrix$normal)%in%rownames(tumorTemp),]
    print(paste0(targetCutoff, " ", nrow(normalTemp[normalTemp$FDR<targetCutoff&normalTemp$deltaMethylation>0,]), 
                 " (", round(nrow(normalTemp[normalTemp$FDR<targetCutoff&normalTemp$deltaMethylation>0,])/nrow(normalTemp),4)*100, "%)"))
    p=Volcanoplot(normalTemp, title, targetCutoff,  color1, color2, paste0(saveIndex, "_EAC_specificPMDs.normal.txt"))
  }else if(type%in%"ESCC"){
    tumorTemp=dataMatrix$tumor[dataMatrix$tumor$FDR<targetCutoff&dataMatrix$tumor$deltaMethylation<0,]
    normalTemp=dataMatrix$normal[rownames(dataMatrix$normal)%in%rownames(tumorTemp),]
    print(paste0(targetCutoff, " ", nrow(normalTemp[normalTemp$FDR<targetCutoff&normalTemp$deltaMethylation<0,]), 
                 " (", round(nrow(normalTemp[normalTemp$FDR<targetCutoff&normalTemp$deltaMethylation<0,])/nrow(normalTemp),4)*100, "%)"))
    p=Volcanoplot(normalTemp, title, targetCutoff,  color1, color2, paste0(saveIndex, "_ESCC_specificPMDs.normal.txt"))
  }
  return(p)
}

# p1=Volcanoplot(EAC_regionLevelPMD$tumor, "EAC-specific PMDs", 0.1, "#0000F5", "#EA3323", "Figure7B/EAC_specificPMDs.txt")
# p2=getSignificant(EAC_regionLevelPMD, "EAC", 0.1, "EAC-specific PMDs", "#75FBFD", "#F7CE46", "Figure7B/Figure7B")
# p3=Volcanoplot(ESCC_regionLevelPMD$tumor, "ESCC-specific PMDs", 0.1, "#0000F5", "#EA3323", "Figure7B/ESCC_specificPMDs.txt")
# p4=getSignificant(ESCC_regionLevelPMD, "ESCC", 0.1, "ESCC-specific PMDs", "#75FBFD", "#F7CE46", "Figure7B/Figure7B")
# p5=Volcanoplot(ESCA_regionLevelSharedPMD$tumor, "shared PMD genes", 0.1, "#0000F5", "#EA3323", "Figure7B/ESCA_sharedPMDs.txt")
# p6=Volcanoplot(ESCA_regionLevelSharedHMD$tumor, "shared HMD genes", 0.1, "#0000F5", "#EA3323", "Figure7B/ESCA_sharedHMDs.txt")

p2=getSignificant(EAC_regionLevelPMD, "EAC", 0.1, "EAC-specific PMDs", "#75FBFD", "#F7CE46", "Figure7B/Figure7B")
p4=getSignificant(ESCC_regionLevelPMD, "ESCC", 0.1, "ESCC-specific PMDs", "#75FBFD", "#F7CE46", "Figure7B/Figure7B")
pdf("Figure7B/Figure7B_specificPMD.volcanoplot.pdf", width=9.5, height=4)
print(ggarrange(p2, p4))
dev.off()

##############################################Figure7D#############################################################
calculatePvalue=function(dataMatrix, sampleList1, sampleList2, sampleList1Name, sampleList2Name){
  TtestResult=lapply(1:nrow(dataMatrix), function(x){
    t.test(dataMatrix[x,colnames(dataMatrix)%in%sampleList1],
           dataMatrix[x,colnames(dataMatrix)%in%sampleList2])$p.value
  })
  TtestResult3=do.call(c, TtestResult)
  
  sampleList1Mean=rowMeans(dataMatrix[,colnames(dataMatrix)%in%sampleList1], na.rm=T)
  sampleList2Mean=rowMeans(dataMatrix[,colnames(dataMatrix)%in%sampleList2], na.rm=T)
  
  TtestResult=data.frame(sampleList1Mean= sampleList1Mean, sampleList2Mean= sampleList2Mean, Pvalue=TtestResult3)
  colnames(TtestResult)[1:2]=c(sampleList1Name, sampleList2Name)
  TtestResult$FDR=p.adjust(TtestResult$Pvalue, method="BH")
  rownames(TtestResult)=rownames(dataMatrix)
  return(TtestResult)
}
getDMRsMethylation=function(DMRMethyMatrixData, dmrData){
  targetDMRsMethylation=DMRMethyMatrixData[rownames(DMRMethyMatrixData)%in%dmrData$DMR,]
  print(nrow(targetDMRsMethylation))
  targetDMRsMethylation=targetDMRsMethylation[rowSums(is.na(targetDMRsMethylation))==0,]
  print(nrow(targetDMRsMethylation))
  ESCC_Tumor_sampleList=paste0("ESCC_", c(1:17,19:22))
  EAC_Tumor_sampleList=c(paste0("EAC_", c(1:4,6:7)), paste0("GEJ_", c(1:7)))
  ESCC_Normal_sampleList=paste0("ESCC_Nonmalignant_", c(1:3, "53F", "54M"))
  GEJ_Normal_sampleList=paste0("GEJ_Nonmalignant_", c(1:7))
  targetDMRsMethylation_TumorPvalue=calculatePvalue(targetDMRsMethylation, ESCC_Tumor_sampleList, EAC_Tumor_sampleList, "ESCCTumor", "EACTumor")
  targetDMRsMethylation_TumorPvalue$deltaMethylation=targetDMRsMethylation_TumorPvalue$ESCCTumor-targetDMRsMethylation_TumorPvalue$EACTumor
  targetDMRsMethylation_NormalPvalue=calculatePvalue(targetDMRsMethylation, ESCC_Normal_sampleList, GEJ_Normal_sampleList, "ESCCNormal", "EACNormal")
  targetDMRsMethylation_NormalPvalue$deltaMethylation=targetDMRsMethylation_NormalPvalue$ESCCNormal-targetDMRsMethylation_NormalPvalue$EACNormal
  result=list()
  result$tumor=targetDMRsMethylation_TumorPvalue
  result$normal=targetDMRsMethylation_NormalPvalue
  return(result)
}
EAC_regionLevelDMRs=getDMRsMethylation(hypogroup1BetaMatrix, hypoGroup1)
ESCC_regionLevelDMRs=getDMRsMethylation(hypogroup2BetaMatrix, hypoGroup2)
Volcanoplot=function(dataMatrix, title, cutoff, color1, color2, saveFile){
  plotdata=dataMatrix
  xlab=paste0(colnames(plotdata)[1:2], collapse = "-")
  plotdata=plotdata[,colnames(plotdata)%in%c("FDR", "deltaMethylation")]
  plotdata$Type="non-change"
  if(nrow(plotdata[plotdata$FDR<cutoff&plotdata$deltaMethylation>0,])>0){
    plotdata[plotdata$FDR<cutoff&plotdata$deltaMethylation>0,]$Type="Up"
  }
  if(nrow(plotdata[plotdata$FDR<cutoff&plotdata$deltaMethylation<0,])>0){
    plotdata[plotdata$FDR<cutoff&plotdata$deltaMethylation<0,]$Type="Down"
  }
  plotdata$FDR=-log10(plotdata$FDR)
  typeColors=data.frame(type=c("Up", "Down", "non-change"), colors=c(color1, color2, "black"), stringsAsFactors = F)
  plotdata$Type=factor(plotdata$Type, levels=c("Up", "Down", "non-change"))
  write.table(plotdata, saveFile, row.names = F, col.names = T, sep="\t", quote=F)
  
  p <- ggplot(data=plotdata,aes(x=deltaMethylation,y=FDR,colour=Type))+geom_point(pch=15)
  p <- p+labs(x = paste0("Delta methlation (", xlab, ")"), y = "-log10 FDR")+ggtitle(title)+scale_color_manual(values=typeColors[typeColors$type%in%unique(plotdata$Type),]$colors)
  xlim1=min(plotdata$deltaMethylation)+(max(plotdata$deltaMethylation)-min(plotdata$deltaMethylation))/10
  xlim2=max(plotdata$deltaMethylation)-(max(plotdata$deltaMethylation)-min(plotdata$deltaMethylation))/10
  ylim=max(plotdata$FDR)-(max(plotdata$FDR)-min(plotdata$FDR))/10
  if(nrow(plotdata[plotdata$Type%in%"Down",])>0){
    p=p+annotate(geom="text", x=xlim1, y=ylim, label=paste0(round(nrow(plotdata[plotdata$Type%in%"Down",])/nrow(plotdata), 4)*100, "%\n(", comma(nrow(plotdata[plotdata$Type%in%"Down",]),0),"/",comma(nrow(plotdata),0),")"),color="black", size=5)
  }
  if(nrow(plotdata[plotdata$Type%in%"Up",])>0){
    p=p+annotate(geom="text", x=xlim2, y=ylim, label=paste0(round(nrow(plotdata[plotdata$Type%in%"Up",])/nrow(plotdata), 4)*100, "%\n(", comma(nrow(plotdata[plotdata$Type%in%"Up",]),0),"/",comma(nrow(plotdata),0),")"),color="black", size=5)
  }
  p <- p +theme_classic() + theme(axis.text = element_text(colour = "black",size=10), 
                                  axis.title=element_text(colour="black",size=12),
                                  plot.title = element_text(hjust = 0.5,size=16))
  # legend.position = "none")
}
getDMRSignificant=function(dataMatrix, type, targetCutoff, title, color1, color2, saveFile){
  if(type%in%"EAC"){
    tumorTemp=dataMatrix$tumor[dataMatrix$tumor$FDR<targetCutoff&dataMatrix$tumor$deltaMethylation>0,]
    normalTemp=dataMatrix$normal[rownames(dataMatrix$normal)%in%rownames(tumorTemp),]
    print(paste0(targetCutoff, " ", nrow(normalTemp[normalTemp$FDR<targetCutoff&normalTemp$deltaMethylation>0,]), 
                 " (", round(nrow(normalTemp[normalTemp$FDR<targetCutoff&normalTemp$deltaMethylation>0,])/nrow(normalTemp),4)*100, "%)"))
    p=Volcanoplot(normalTemp, title, targetCutoff,  color1, color2, saveFile)
  }else if(type%in%"ESCC"){
    tumorTemp=dataMatrix$tumor[dataMatrix$tumor$FDR<targetCutoff&dataMatrix$tumor$deltaMethylation<0,]
    normalTemp=dataMatrix$normal[rownames(dataMatrix$normal)%in%rownames(tumorTemp),]
    print(paste0(targetCutoff, " ", nrow(normalTemp[normalTemp$FDR<targetCutoff&normalTemp$deltaMethylation<0,]), 
                 " (", round(nrow(normalTemp[normalTemp$FDR<targetCutoff&normalTemp$deltaMethylation<0,])/nrow(normalTemp),4)*100, "%)"))
    p=Volcanoplot(normalTemp, title, targetCutoff,  color1, color2, saveFile)
  }
  return(p)
}
p1=getDMRSignificant(EAC_regionLevelDMRs, "EAC", 0.1, "EAC/GEJ tumor hypoDMRs", "#75FBFD", "#F7CE46", "Figure7D/Figure7D_EAC.txt")
p2=getDMRSignificant(ESCC_regionLevelDMRs, "ESCC", 0.1, "ESCC tumor hypoDMRs", "#75FBFD", "#F7CE46", "Figure7D/Figure7D_ESCC.txt")
pdf("Figure7D/Figure7D.volcanoplot.pdf", width=9.5, height=4)
print(ggarrange(p1,p2))
dev.off()


##############################################Figure7EFG_S7A#############################################################
###Figure7E_left
###download bed file from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184462
targetTissueType=read.table("Data/Figure7EFG_S7A/GSE184462_tissueType.txt", header=T, sep="\t")
targetEpithelialCells=targetTissueType[targetTissueType$Cells%in%"Epithelial cells",]$Tissue_contribution
metaInfo=read.table("Data/Figure7EFG_S7A/GSE184462_metadata.tsv", sep="\t", stringsAsFactors = F, header=T)
targetMetaInfo=metaInfo[metaInfo$cell.type%in%targetEpithelialCells,]
targetMetaInfo$barcodeID=gsub(".*\\+", "", targetMetaInfo$cellID)

for(PMDfile in dir("Data/Figure7EFG_S7A/PMDresult/", pattern = ".txt", full.names = T)){
  tmpData=read.table(PMDfile, header=T, sep="\t", stringsAsFactors = F)
  if(PMDfile%in%dir("Data/Figure7EFG_S7A/PMDresult/", pattern = ".txt", full.names = T)[1]){
    PMD_result=tmpData
  }else{
    PMD_result=rbind(PMD_result, tmpData)
  }
}
write.table(PMD_result, "Data/Figure7EFG_S7A/PMD_result.txt", row.names = F, col.names = T, sep="\t", quote=F)
PMD_result=merge(PMD_result, targetMetaInfo[,c("cellID", "tissue", "barcodeID")], by=c("barcodeID", "tissue"))
umapData=read.table("Data/Figure7EFG_S7A/adult.tsv", header=T, sep="\t", stringsAsFactors = F)
colnames(umapData)=c("cellID", "UMAP1", "UMAP2")
umapData=merge(umapData, PMD_result, by="cellID")
umapData=umapData[,c(1:3,6,8:11)]
umapData=merge(umapData, targetTissueType[,c(1,3,4)], by.x="cell.type", by.y="Tissue_contribution")
umapData$ESCC_EAC_specificPMDs=umapData$ESCC_specificPMDs-umapData$EAC_specificPMDs
umapData$HMD_PMD=umapData$ESCA_sharedHMDs-umapData$ESCA_sharedPMDs
umapData$Source=factor(umapData$Source, levels=c("GI", "squamous", "other"))
write.table(umapData, "Figure7EFG_S7A/Figure7EF_S6C.data.txt", row.names = F, col.names = T, sep="\t", quote=F)

umapData2=umapData
umapData2[umapData2$ESCC_EAC_specificPMDs>1,]$ESCC_EAC_specificPMDs=1
umapData2[umapData2$ESCC_EAC_specificPMDs<(-1),]$ESCC_EAC_specificPMDs=(-1)
if(nrow(umapData2[umapData2$HMD_PMD>2,])>0){
  umapData2[umapData2$HMD_PMD>2,]$HMD_PMD=2
}
if(nrow(umapData2[umapData2$HMD_PMD<(-2),])>0){
  umapData2[umapData2$HMD_PMD<(-2),]$HMD_PMD=(-2)
}
p1=ggplot(umapData2, aes(x=UMAP1, y=UMAP2, color=Source)) + geom_point(alpha=0.7,size=0.1)
p1=p1+theme_classic()+scale_color_manual(values=c("blue", "red", "grey"))
p1=p1+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p1=p1+theme(legend.position = "none")
png("Figure7EFG_S7A/Figure7E.cluster.png", res=300, height=1600, width=2000)
print(p1)
dev.off()

###Figure7E_middle
p2=ggplot(umapData2, aes(x=UMAP1, y=UMAP2, color=ESCC_EAC_specificPMDs)) + geom_point(alpha=0.7,size=0.1)
p2=p2+theme_classic()+scale_color_gradient2(low ="red", mid = "white", high ="blue", limits=c(-1, 1), 
                                            midpoint = 0, name="delta methylation\n(ESCC_tumor_PMDs \nV.S. EAC_tumor_PMDs)")
p2=p2+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
pdf("Figure7EFG_S7A/Figure7E.middle.pdf", height=5, width=7)
print(p2)
dev.off()
p2=p2+theme(legend.position = "none")
png("Figure7EFG_S7A/Figure7E.middle.png", res=300, height=1700, width=2000)
print(p2)
dev.off()

###FigureS7A
p3=ggplot(umapData2, aes(x=UMAP1, y=UMAP2, color=HMD_PMD)) + geom_point(alpha=0.7,size=0.1)
p3=p3+theme_classic()+scale_color_gradient2(low ="red", mid = "white", high ="blue", limits=c(-2, 2), 
                                            midpoint = 0, name="delta methylation\n(shared_HMDs \nV.S. shared_PMDs)")
p3=p3+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
pdf("Figure7EFG_S7A/FigureS7A.pdf", height=5, width=7)
print(p3)
dev.off()
p3=p3+theme(legend.position = "none")
png("Figure7EFG_S7A/FigureS7A.png", res=300, height=1600, width=2000)
print(p3)
dev.off()

###Figure7F
meanData= do.call(c, lapply(unique(umapData$Subtype), function(x){
  mean(umapData[umapData$Subtype%in%x,]$ESCC_EAC_specificPMDs)
}))
meanData=data.frame(Subtype=unique(umapData$Subtype), mean=meanData)
meanData=meanData[order(meanData$mean),]
umapData$Subtype=factor(umapData$Subtype, levels=meanData$Subtype)
plotScATACPointPlot=function(plotdata, ylabel, saveFile){
  plotdata=plotdata[order(plotdata$cluster, plotdata$delta),]
  name<-as.character(unique(plotdata$cluster))
  sum=0
  for(i in 1:length(name)){
    tumorRow = nrow(plotdata[plotdata$cluster==name[i],])
    for(j in 1:tumorRow){
      plotdata[sum+j,"Postion"]<- (1000*(i-1)+350+250*j/tumorRow)
    }
    sum=sum+tumorRow
  } 
  for(i in 1:length(name)){
    plotdata[plotdata$cluster==name[i],"Tick"]<-(1000*(i-1)+350)
  }
  meanData=plotdata[1,]
  sum2=0
  for(i in 1:length(name)){
    tumorRow = nrow(plotdata[plotdata$cluster==name[i],])
    temp=plotdata[sum2+ceiling(tumorRow/2),]
    temp$delta=mean(plotdata[plotdata$cluster==name[i],]$delta)
    meanData <-rbind(meanData,temp)
    sum2=sum2+tumorRow
  }
  meanData$xmin=meanData$Postion-125
  meanData$xmax=meanData$Postion+125
  meanData$ymin=meanData$delta
  meanData$ymax=meanData$delta
  meanData=meanData[-1,]
  p = ggplot(plotdata, aes(x=Postion, y=delta, color=source))+geom_point(size=1,stat="identity")+
    geom_rect(data=meanData,aes(xmin = xmin, xmax = xmax, ymin = delta, ymax = delta),
              fill="white",size=1,color="orange")
  p=p+theme_classic()+ylab(ylabel)+ggtitle("scATACseq")
  p=p+scale_color_manual(name="cell types", values=c("blue", "red", "grey"))
  p=p+scale_x_continuous(name = "", breaks = unique(plotdata$Tick), labels = unique(plotdata$cluster))
  p=p+geom_hline(yintercept=0, linetype="dashed", color = "grey")
  p=p+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text = element_text(color="black", size=12),
            plot.title = element_text(hjust = 0.5, size=16, face="bold"), legend.position = "bottom",
            legend.text = element_text(size=12))
  png(saveFile, width=2000, height=1500,res=300)
  print(p)
  dev.off()
}
plotdata=umapData[,c("Subtype", "ESCC_EAC_specificPMDs", "Source")]
colnames(plotdata)=c("cluster", "delta", "source")
plotScATACPointPlot(plotdata, "delta accessiblity (ESCC-specific \nPMDs V.S. EAC-specific PMDs)", "Figure7EFG_S7A/Figure7F.png")

######Figure7E_right
for(DMRfile in dir("Data/Figure7EFG_S7A/DMRresult//", pattern = ".txt", full.names = T)){
  tmpData=read.table(DMRfile, header=T, sep="\t", stringsAsFactors = F)
  if(DMRfile%in%dir("Data/Figure7EFG_S7A/DMRresult//", pattern = ".txt", full.names = T)[1]){
    DMR_result=tmpData
  }else{
    DMR_result=rbind(DMR_result, tmpData)
  }
}
write.table(DMR_result, "Data/Figure7EFG_S7A/DMR_result.txt", row.names = F, col.names = T, sep="\t", quote=F)
DMR_result=merge(DMR_result, targetMetaInfo[,c("cellID", "tissue", "barcodeID")], by=c("barcodeID", "tissue"))
umapData=read.table("Data/Figure7EFG_S7A/adult.tsv", header=T, sep="\t", stringsAsFactors = F)
colnames(umapData)=c("cellID", "UMAP1", "UMAP2")
umapData=merge(umapData, DMR_result, by="cellID")
umapData=umapData[,c(1:3,6,8:9)]
umapData=merge(umapData, targetTissueType[,c(1,3,4)], by.x="cell.type", by.y="Tissue_contribution")
umapData$ESCC_EAC_DMRs=umapData$ESCC_DMRs-umapData$EAC_DMRs
umapData$Source=factor(umapData$Source, levels=c("GI", "squamous", "other"))
write.table(umapData, "Data/Figure7EFG_S7A/Figure7EG.data.txt", row.names = F, col.names = T, sep="\t", quote=F)
umapData2=umapData
umapData2[umapData2$ESCC_EAC_DMRs>10,]$ESCC_EAC_DMRs=10
umapData2[umapData2$ESCC_EAC_DMRs<(-10),]$ESCC_EAC_DMRs=(-10)
p2=ggplot(umapData2, aes(x=UMAP1, y=UMAP2, color=ESCC_EAC_DMRs)) + geom_point(alpha=0.7,size=0.1)
p2=p2+theme_classic()+scale_color_gradient2(low ="blue", mid = "white", high ="red", limits=c(-10, 10), 
                                            midpoint = 0, name="delta methylation\n(ESCC_tumor_PMDs \nV.S. EAC_tumor_PMDs)")
p2=p2+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
pdf("Figure7EFG_S7A/Figure7E.right.pdf", height=5, width=7)
print(p2)
dev.off()
p2=p2+theme(legend.position = "none")
png("Figure7EFG_S7A/Figure7E.right.png", res=300, height=2000, width=2000)
print(p2)
dev.off()

######Figure7G
meanData= do.call(c, lapply(unique(umapData$Subtype), function(x){
  mean(umapData[umapData$Subtype%in%x,]$ESCC_EAC_DMRs)
}))
meanData=data.frame(Subtype=unique(umapData$Subtype), mean=meanData)
meanData=meanData[order(meanData$mean),]
umapData$Subtype=factor(umapData$Subtype, levels=meanData$Subtype)
plotdata=umapData[,c("Subtype", "ESCC_EAC_DMRs", "Source")]
colnames(plotdata)=c("cluster", "delta", "source")
plotScATACPointPlot(plotdata, "delta accessiblity (ESCC-specific \nDMRs V.S. EAC-specific DMRs)", "Figure7EFG_S7A/Figure7G.png")

##############################################Figure8ABCDE_S7B#############################################################
######Figure8DE
cancerSubtypeInfo=read.table("Data/Figure8ABCDE_S7B/TCGA_cancerSample.txt", sep="\t", stringsAsFactors = F, header=T)
bladderSubtypeInfo=read.table("Data/Figure8ABCDE_S7B/TCGA_bladderSample.txt", sep="\t", stringsAsFactors = F, header=T)
bladderSubtypeInfo$PATIENT_BARCODE=do.call(c, lapply(bladderSubtypeInfo$ID, function(x){
  paste0(strsplit(x, "-")[[1]][1:3], collapse = "-")
}))
bladderSubtypeInfo$SAMPLE_BARCODE=bladderSubtypeInfo$ID
bladderSubtypeInfo$DISEASE="BLCA"
bladderSubtypeInfo$SUBTYPE=bladderSubtypeInfo$Consensus_class
bladderSubtypeInfo=bladderSubtypeInfo[,c(4:7)]
cancerSubtypeInfo=cancerSubtypeInfo[!cancerSubtypeInfo$DISEASE%in%"BLCA",]
cancerSubtypeInfo=rbind(cancerSubtypeInfo, bladderSubtypeInfo)
cancerSubtypeInfo$Label="Others"
cancerSubtypeInfo[cancerSubtypeInfo$DISEASE%in%"COAD",]$Label="COAD"
cancerSubtypeInfo[cancerSubtypeInfo$DISEASE%in%"READ",]$Label="READ"
cancerSubtypeInfo[cancerSubtypeInfo$DISEASE%in%"HNSC",]$Label="HNSC"
cancerSubtypeInfo[cancerSubtypeInfo$DISEASE%in%"LUSC",]$Label="LUSC"
cancerSubtypeInfo[cancerSubtypeInfo$DISEASE%in%"STAD",]$Label="STAD"
cancerSubtypeInfo[cancerSubtypeInfo$DISEASE%in%"BLCA"&cancerSubtypeInfo$SUBTYPE%in%"Ba/Sq",]$Label="BLCA (Ba/Sq)"
cancerSubtypeInfo[cancerSubtypeInfo$DISEASE%in%"CESC"&cancerSubtypeInfo$SUBTYPE%in%"SquamousCarcinoma",]$Label="CESC (Sq)"
load("meta/TCGA-ESCA.clinc.RData")
EAC_id=sampleInformation[sampleInformation$Label%in%"EAC_Tumor",1]
ESCC_id=sampleInformation[sampleInformation$Label%in%"ESCC_Tumor",1]

######Figure8D
plotMutationPointPlot=function(plotdata2, ylabel, saveFile){
  colnames(plotdata2)[4]=c("delta")
  plotdata2=plotdata2[order(plotdata2$cancerType, plotdata2$delta),]
  name<-as.character(unique(plotdata2$cancerType))
  sum=0
  for(i in 1:length(name)){
    tumorRow = nrow(plotdata2[plotdata2$cancerType==name[i],])
    for(j in 1:tumorRow){
      plotdata2[sum+j,"Postion"]<- (1000*(i-1)+350+250*j/tumorRow)
    }
    sum=sum+tumorRow
  } 
  for(i in 1:length(name)){
    plotdata2[plotdata2$cancerType==name[i],"Tick"]<-(1000*(i-1)+350)
  }
  meanData=plotdata2[1,]
  sum2=0
  for(i in 1:length(name)){
    tumorRow = nrow(plotdata2[plotdata2$cancerType==name[i],])
    temp=plotdata2[sum2+ceiling(tumorRow/2),]
    temp$delta=mean(plotdata2[plotdata2$cancerType==name[i],]$delta)
    meanData <-rbind(meanData,temp)
    sum2=sum2+tumorRow
  }
  meanData$xmin=meanData$Postion-125
  meanData$xmax=meanData$Postion+125
  meanData$ymin=meanData$delta
  meanData$ymax=meanData$delta
  meanData=meanData[-1,]
  p = ggplot(plotdata2, aes(x=Postion, y=delta, color=Label))+geom_point(size=0.1,stat="identity")+
    geom_rect(data=meanData,aes(xmin = xmin, xmax = xmax, ymin = delta, ymax = delta),
              fill="white",size=0.4,color="orange")
  p=p+theme_classic()+ylab(ylabel)+ggtitle("TCGA pan-cancer methylation")
  p=p+scale_color_manual(name="Cancer groups", values=c("#0000F5", "black", "#EA3323"))
  p=p+scale_x_continuous(name = "", breaks = unique(plotdata2$Tick), labels = unique(plotdata2$cancerType))
  p=p+geom_hline(yintercept=0, linetype="dashed", color = "grey")
  p=p+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text = element_text(color="black", size=12),
            plot.title = element_text(hjust = 0.5, size=16, face="bold"), legend.position = "bottom",
            legend.text = element_text(size=12))
  pdf(saveFile, width=9, height=4)
  print(p)
  dev.off()
}
getMethylationHT450k=function(regionFile, probesFile, type, methylationMatrix){
  tempData=read_bed(regionFile,n_fields = 5)
  probesMatrix=read_bed(probesFile, n_fields = 4)
  tempData_probes=bed_intersect(probesMatrix,tempData)
  tempData_probes=as.data.frame(tempData_probes)
  tempData_probes=unique(tempData_probes[tempData_probes$.overlap==2,]$name.x)
  print(length(tempData_probes))
  tempData=methylationMatrix[rownames(methylationMatrix)%in%tempData_probes,, drop=F]
  tempData=data.frame(colMeans(tempData, na.rm = T))
  colnames(tempData)=type
  return(tempData)
}
getAllCancerMethylation=function(cancerSampleInfo, saveFile){
  # probesFile="meta/HT450k.probe.rmCGI.rmblackList.solo.bed"
  # sampleNum=data.frame(table(cancerSampleInfo$DISEASE),stringsAsFactors = F)
  # resultMatrix= as.data.frame(matrix(numeric(0),ncol=5))
  # sampleInfoMatrix=as.data.frame(matrix(numeric(0),ncol=3))
  # for(cancerType in as.character(sampleNum$Var1)){
  #   print(cancerType)
  #   load(paste0("/Volumes/Yuan_2T/TCGA/Methylation_array/Methylation/TCGA-", cancerType, "_methylation_hg38.v2.RData"))
  #   TCGAprobeMethylationMatrix=result
  #   TCGAprobeSampleInfo=match.file.cases
  #   if(cancerType%in%"CHOL"){
  #     TCGAprobeSampleInfo=TCGAprobeSampleInfo[,c(1,4,3)]
  #   }else{
  #     TCGAprobeSampleInfo=TCGAprobeSampleInfo[,c(1,4,2)]
  #   }
  #   colnames(TCGAprobeSampleInfo)=c("Sample","SampleName", "Type")
  #   print(paste0(cancerType, " ", nrow(TCGAprobeSampleInfo)))
  #   TCGAprobeMethylationMatrix=TCGAprobeMethylationMatrix[,colnames(TCGAprobeMethylationMatrix)%in%TCGAprobeSampleInfo$SampleName]
  #   temp1=suppressWarnings(getMethylationHT450k("Data/MMSeekR_PMDs/ESCC_specificPMDs.bed", probesFile, "ESCC_specificPMDs",  TCGAprobeMethylationMatrix))
  #   temp2=suppressWarnings(getMethylationHT450k("Data/MMSeekR_PMDs/EAC_specificPMDs.bed", probesFile, "EAC_specificPMDs", TCGAprobeMethylationMatrix))
  #   temp3=suppressWarnings(getMethylationHT450k("Data/MMSeekR_PMDs/ESCA_sharedPMDs.bed", probesFile, "SharedPMDs", TCGAprobeMethylationMatrix))
  #   temp4=suppressWarnings(getMethylationHT450k("Data/MMSeekR_PMDs/ESCA_sharedHMDs.bed", probesFile, "SharedHMDs", TCGAprobeMethylationMatrix))
  #   temp=cbind(temp1, temp2, temp3, temp4)
  #   temp$cancerType=cancerType
  #   resultMatrix=rbind(resultMatrix, temp)
  #   sampleInfoMatrix=rbind(sampleInfoMatrix, TCGAprobeSampleInfo)
  # }
  # sampleInfoMatrix2=sampleInfoMatrix
  # resultMatrix2=resultMatrix
  # resultMatrix$SampleName=rownames(resultMatrix)
  # resultMatrix=merge(resultMatrix, sampleInfoMatrix, by.x="SampleName", by.y="SampleName", all.x=T)
  # resultMatrix=merge(resultMatrix, cancerSubtypeInfo[,c(1,5)], by.x="Sample", by.y="PATIENT_BARCODE", all.x=T)
  # resultMatrix[is.na(resultMatrix$Label),]$Label="Others"
  # resultMatrix$delta_ESCCspecificVSEACspecific=resultMatrix$ESCC_specificPMDs-resultMatrix$EAC_specificPMDs
  # resultMatrix$delta_sharedHMDVSsharedPMD= resultMatrix$SharedHMDs-resultMatrix$SharedPMDs
  # 
  # tumorResultMatrix=resultMatrix[resultMatrix$Type%in%"Tumor",]
  # tumorResultMatrix[tumorResultMatrix$Sample%in%EAC_id,]$cancerType="EAC"
  # tumorResultMatrix[tumorResultMatrix$Sample%in%ESCC_id,]$cancerType="ESCC"
  # tumorResultMatrix=tumorResultMatrix[!tumorResultMatrix$cancerType%in%"ESCA",]
  # tumorResultMatrix[tumorResultMatrix$Label%in%"BLCA (Ba/Sq)",]$cancerType="BLCA (Ba/Sq)"
  # tumorResultMatrix[tumorResultMatrix$cancerType%in%"BLCA",]$cancerType="BLCA (Others)"
  # tumorResultMatrix[tumorResultMatrix$Label%in%"CESC (Sq)",]$cancerType="CESC (Sq)"
  # tumorResultMatrix[tumorResultMatrix$cancerType%in%"CESC",]$cancerType="CESC (Others)"
  # tumorResultMatrix=tumorResultMatrix[,c(2,7,10,11)]
  # save(tumorResultMatrix, file="Data/Figure8ABCDE_S7B/TCGA_methylation.Rdata")
  
  load("Data/Figure8ABCDE_S7B/TCGA_methylation.Rdata")
  medianData= do.call(c, lapply(unique(tumorResultMatrix$cancerType), function(x){
    mean(tumorResultMatrix[tumorResultMatrix$cancerType%in%x,]$delta_ESCCspecificVSEACspecific)
  }))
  medianData=data.frame(cancerType=unique(tumorResultMatrix$cancerType), median=unique(medianData))
  medianData=medianData[order(medianData$median),]
  tumorResultMatrix$cancerType=factor(tumorResultMatrix$cancerType, levels=as.character(medianData$cancerType))
  tumorResultMatrix$Label="Other cancers"
  tumorResultMatrix[tumorResultMatrix$cancerType%in%c("ESCC", "HNSC", "LUSC", "CESC (Sq)", "BLCA (Ba/Sq)"),]$Label="Squamous cancers"
  tumorResultMatrix[tumorResultMatrix$cancerType%in%c("EAC", "COAD", "STAD", "READ"),]$Label="GI cancers"
  tumorResultMatrix$Label=factor(tumorResultMatrix$Label, levels=c("GI cancers", "Other cancers", "Squamous cancers"))
  tumorResultMatrix[,4]=tumorResultMatrix[,3]
  write.table(tumorResultMatrix, gsub(".pdf", ".data.txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
  plotMutationPointPlot(tumorResultMatrix, "Methyation(ESCC only PMDs v.s. EAC only PMDs)", saveFile)
}
getAllCancerMethylation(cancerSubtypeInfo, "Figure8ABCDE_S7B/Figure8D.pdf")

######Figure8E
getAllCancerMethylationDMRs=function(cancerSampleInfo, saveFile){
  # probesFile="meta/HT450k.probe.rmblackList.bed"
  # sampleNum=data.frame(table(cancerSampleInfo$DISEASE),stringsAsFactors = F)
  # resultMatrix= as.data.frame(matrix(numeric(0),ncol=3))
  # sampleInfoMatrix=as.data.frame(matrix(numeric(0),ncol=3))
  # for(cancerType in as.character(sampleNum$Var1)){
  #   print(cancerType)
  #   load(paste0("/Volumes/Yuan_2T/TCGA/Methylation_array/Methylation/TCGA-", cancerType, "_methylation_hg38.v2.RData"))
  #   TCGAprobeMethylationMatrix=result
  #   TCGAprobeSampleInfo=match.file.cases
  #   if(cancerType%in%"CHOL"){
  #     TCGAprobeSampleInfo=TCGAprobeSampleInfo[,c(1,4,3)]
  #   }else{
  #     TCGAprobeSampleInfo=TCGAprobeSampleInfo[,c(1,4,2)]
  #   }
  #   colnames(TCGAprobeSampleInfo)=c("Sample","SampleName", "Type")
  #   print(paste0(cancerType, " ", nrow(TCGAprobeSampleInfo)))
  #   TCGAprobeMethylationMatrix=TCGAprobeMethylationMatrix[,colnames(TCGAprobeMethylationMatrix)%in%TCGAprobeSampleInfo$SampleName]
  #   temp1=suppressWarnings(getMethylationHT450k("Data/MaskUnionPMDs_DMRs/hypoESCC_Tumor.bed", probesFile, "ESCC_DMRs",  TCGAprobeMethylationMatrix))
  #   temp2=suppressWarnings(getMethylationHT450k("Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed", probesFile, "EAC_DMRs", TCGAprobeMethylationMatrix))
  #   temp=cbind(temp1, temp2)
  #   temp$cancerType=cancerType
  #   resultMatrix=rbind(resultMatrix, temp)
  #   sampleInfoMatrix=rbind(sampleInfoMatrix, TCGAprobeSampleInfo)
  # }
  # sampleInfoMatrix2=sampleInfoMatrix
  # resultMatrix2=resultMatrix
  # resultMatrix$SampleName=rownames(resultMatrix)
  # resultMatrix=merge(resultMatrix, sampleInfoMatrix, by.x="SampleName", by.y="SampleName", all.x=T)
  # resultMatrix=merge(resultMatrix, cancerSubtypeInfo[,c(1,5)], by.x="Sample", by.y="PATIENT_BARCODE", all.x=T)
  # resultMatrix[is.na(resultMatrix$Label),]$Label="Others"
  # resultMatrix$delta_ESCCDMRsVSEACDMRs=resultMatrix$ESCC_DMRs-resultMatrix$EAC_DMRs
  # tumorResultMatrix=resultMatrix[resultMatrix$Type%in%"Tumor",]
  # tumorResultMatrix[tumorResultMatrix$Sample%in%EAC_id,]$cancerType="EAC"
  # tumorResultMatrix[tumorResultMatrix$Sample%in%ESCC_id,]$cancerType="ESCC"
  # tumorResultMatrix=tumorResultMatrix[!tumorResultMatrix$cancerType%in%"ESCA",]
  # tumorResultMatrix[tumorResultMatrix$Label%in%"BLCA (Ba/Sq)",]$cancerType="BLCA (Ba/Sq)"
  # tumorResultMatrix[tumorResultMatrix$cancerType%in%"BLCA",]$cancerType="BLCA (Others)"
  # tumorResultMatrix[tumorResultMatrix$Label%in%"CESC (Sq)",]$cancerType="CESC (Sq)"
  # tumorResultMatrix[tumorResultMatrix$cancerType%in%"CESC",]$cancerType="CESC (Others)"
  # tumorResultMatrix=tumorResultMatrix[,c(2,5,7,8)]
  # save(tumorResultMatrix, file="Data/Figure8ABCDE_S7B/TCGA_methylation_AllDMR.Rdata")
  load("Data/Figure8ABCDE_S7B/TCGA_methylation_AllDMR.Rdata")
  medianData= do.call(c, lapply(unique(tumorResultMatrix$cancerType), function(x){
    mean(tumorResultMatrix[tumorResultMatrix$cancerType%in%x,]$delta_ESCCDMRsVSEACDMRs)
  }))
  medianData=data.frame(cancerType=unique(tumorResultMatrix$cancerType), median=unique(medianData))
  medianData=medianData[order(medianData$median),]
  tumorResultMatrix$cancerType=factor(tumorResultMatrix$cancerType, levels=as.character(medianData$cancerType))
  tumorResultMatrix$Label="Other cancers"
  tumorResultMatrix[tumorResultMatrix$cancerType%in%c("ESCC", "HNSC", "LUSC", "CESC (Sq)", "BLCA (Ba/Sq)"),]$Label="Squamous cancers"
  tumorResultMatrix[tumorResultMatrix$cancerType%in%c("EAC", "COAD", "STAD", "READ"),]$Label="GI cancers"
  tumorResultMatrix$Label=factor(tumorResultMatrix$Label, levels=c("GI cancers", "Other cancers", "Squamous cancers"))
  write.table(tumorResultMatrix, gsub(".pdf", ".data.txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
  plotMutationPointPlot(tumorResultMatrix, "Methyation(ESCC tumor hypoDMRs v.s. EAC tumor hypoDMRs)",saveFile)
}
getAllCancerMethylationDMRs(cancerSubtypeInfo, "Figure8ABCDE_S7B/Figure8E.pdf")

#####Figure8A
tumorColorMap=read.table("Data/Figure8ABCDE_S7B/PancanAtlas_tumor_colors.tab", sep="\t", stringsAsFactors = F, header=T, comment.char = "@")
pancanAtlas_iCluster=read.table("Data/Figure8ABCDE_S7B/PancanAtlas_euclideaniCluster_hexagonPos.tab", sep="\t", header=F, stringsAsFactors = F)
colnames(pancanAtlas_iCluster)=c("samples", "x", "y")
pancanAtlas_iCluster_samples=read.table("Data/Figure8ABCDE_S7B/PancanAtlas_euclideaniCluster_sample.txt", sep="\t", stringsAsFactors = F, header=T)
pancanAtlas_iCluster=merge(pancanAtlas_iCluster, pancanAtlas_iCluster_samples, by.x="samples", by.y="samples")
pancanAtlas_iCluster$type=do.call(c, lapply(pancanAtlas_iCluster$samples, function(x){as.numeric(strsplit(x, "-")[[1]][4])}))
pancanAtlas_iCluster$barcode=do.call(c, lapply(pancanAtlas_iCluster$samples, function(x){paste0(strsplit(x, "-")[[1]][1:3], collapse = "-")}))
pancanAtlas_iCluster=merge(pancanAtlas_iCluster, tumorColorMap, by.x="disease", by.y="disease")
pancanAtlas_iCluster=pancanAtlas_iCluster[!duplicated(pancanAtlas_iCluster$barcode),]

load("Data/Figure8ABCDE_S7B/TCGA_methylation.Rdata")
tumorResultMatrix$SampleName=gsub("_T$", "", tumorResultMatrix$SampleName)
pancanAtlas_iCluster_methylation=merge(pancanAtlas_iCluster, tumorResultMatrix, by.x="barcode", by.y="SampleName")
pancanAtlas_iCluster_methylation[pancanAtlas_iCluster_methylation$cancerType%in%"BLCA (Ba/Sq)",]$color="green"
pancanAtlas_iCluster_methylation[pancanAtlas_iCluster_methylation$cancerType%in%"CESC (Sq)",]$color="#98FF98"
pancanAtlas_iCluster_methylation[pancanAtlas_iCluster_methylation$cancerType%in%"ESCC",]$color="#254117"
pancanAtlas_iCluster_methylation=pancanAtlas_iCluster_methylation[order(pancanAtlas_iCluster_methylation$cancerType),]
pancanAtlas_iCluster_methylation$PanCan="Others"
pancanAtlas_iCluster_methylation[pancanAtlas_iCluster_methylation$disease%in%c("BLCA (Ba/Sq)", "CESC (Sq)", "ESCC", "HNSC", "LUSC"),]$PanCan="Pan-Squamous"
pancanAtlas_iCluster_methylation[pancanAtlas_iCluster_methylation$disease%in%c("EAC", "COAD", "STAD", "READ"),]$PanCan="Pan-GI"
pancanAtlas_iCluster_methylation$PanCan=factor(pancanAtlas_iCluster_methylation$PanCan, levels=c("Others","Pan-Squamous","Pan-GI"))
write.table(pancanAtlas_iCluster_methylation, file="Figure8ABCDE_S7B/Figure8AB_S7B.data.txt", row.names = F, col.names = T, sep="\t", quote=F)

p1=ggplot(pancanAtlas_iCluster_methylation, aes(x=x, y=y, color=cancerType)) + geom_point(alpha=0.7,size=0.1)
p1=p1+theme_classic()+scale_color_manual(values=unique(pancanAtlas_iCluster_methylation$color))
p1=p1+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p1=p1+theme(legend.position = "none")
p1=ggplot(pancanAtlas_iCluster_methylation, aes(x=x, y=y, color=PanCan)) + geom_point(alpha=0.7,size=0.8)
p1=p1+theme_classic()+scale_color_manual(values=c("grey", "#EA3323", "#0000F5"))
p1=p1+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p1=p1+theme(legend.position = "none")
png("Figure8ABCDE_S7B/Figure8A.png", res=300, width=1800, height=1800)
print(p1)
dev.off()

######Figure8B
pancanAtlas_iCluster_methylation2=pancanAtlas_iCluster_methylation[order(pancanAtlas_iCluster_methylation$delta_ESCCspecificVSEACspecific),]
pancanAtlas_iCluster_methylation2_1=pancanAtlas_iCluster_methylation2[pancanAtlas_iCluster_methylation2$delta_ESCCspecificVSEACspecific<(-0.2)|pancanAtlas_iCluster_methylation2$delta_ESCCspecificVSEACspecific>0.2,]
pancanAtlas_iCluster_methylation2_2=pancanAtlas_iCluster_methylation2[(pancanAtlas_iCluster_methylation2$delta_ESCCspecificVSEACspecific<(-0.1)&pancanAtlas_iCluster_methylation2$delta_ESCCspecificVSEACspecific>(-0.2))
                                                                      |(pancanAtlas_iCluster_methylation2$delta_ESCCspecificVSEACspecific>0.1&pancanAtlas_iCluster_methylation2$delta_ESCCspecificVSEACspecific<0.2),]
pancanAtlas_iCluster_methylation2_3=pancanAtlas_iCluster_methylation2[(pancanAtlas_iCluster_methylation2$delta_ESCCspecificVSEACspecific>(-0.1)&pancanAtlas_iCluster_methylation2$delta_ESCCspecificVSEACspecific<0.1),]
p2=ggplot(pancanAtlas_iCluster_methylation2_3, aes(x=x, y=y, color=delta_ESCCspecificVSEACspecific)) + geom_point(alpha=0.7,size=0.8)
p2=p2+geom_point(data=pancanAtlas_iCluster_methylation2_2, aes(x=x, y=y, color=delta_ESCCspecificVSEACspecific),alpha=0.7,size=0.8)
p2=p2+geom_point(data=pancanAtlas_iCluster_methylation2_1, aes(x=x, y=y, color=delta_ESCCspecificVSEACspecific),alpha=0.7,size=0.8)
p2=p2+scale_colour_gradientn(colours=colorRampPalette(c("#FF0000", "#F22626","#D87272", "grey", "grey", "grey", "#7171D8", "#2626F2", "#0000FF"))(50),
                             limits=c(-0.33, 0.33), name="delta methylation\n(ESCC_tumor_PMDs \nV.S. EAC_tumor_PMDs)")+theme_classic()
p2=p2+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p3=p2+theme(legend.position = "none")
png("Figure8ABCDE_S7B/Figure8B.png", res=300, width=1800, height=1800)
print(p3)
dev.off()
pdf("Figure8ABCDE_S7B/Figure8B.legend.pdf", width=6, height=6)
print(p2)
dev.off()

######FigureS7B
p2=ggplot(pancanAtlas_iCluster_methylation, aes(x=x, y=y, color=delta_sharedHMDVSsharedPMD)) + geom_point(alpha=0.7,size=0.1)
p2=p2+theme_classic()+scale_color_gradient2(low ="red", mid = "white", high ="blue", limits=c(-0.44, 0.44),
                                            midpoint = 0, name="delta methylation\n(shared HMDs \nV.S. shared PMDs)")
p2=p2+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p3=p2+theme(legend.position = "none")
png("Figure8ABCDE_S7B/FigureS7B.png", res=300, width=1800, height=1800)
print(p3)
dev.off()
pdf("Figure8ABCDE_S7B/FigureS7B.legend.pdf", width=6, height=6)
print(p2)
dev.off()

######Figure8C
load("Data/Figure8ABCDE_S7B/TCGA_methylation_AllDMR.Rdata")
tumorResultMatrix$SampleName=gsub("_T$", "", tumorResultMatrix$SampleName)
pancanAtlas_iCluster_methylation=merge(pancanAtlas_iCluster, tumorResultMatrix, by.x="barcode", by.y="SampleName")
pancanAtlas_iCluster_methylation[pancanAtlas_iCluster_methylation$cancerType%in%"BLCA (Ba/Sq)",]$color="green"
pancanAtlas_iCluster_methylation[pancanAtlas_iCluster_methylation$cancerType%in%"CESC (Sq)",]$color="#98FF98"
pancanAtlas_iCluster_methylation[pancanAtlas_iCluster_methylation$cancerType%in%"ESCC",]$color="#254117"
pancanAtlas_iCluster_methylation=pancanAtlas_iCluster_methylation[order(pancanAtlas_iCluster_methylation$cancerType),]
pancanAtlas_iCluster_methylation$cancerType=factor(pancanAtlas_iCluster_methylation$cancerType, levels=unique(pancanAtlas_iCluster_methylation$cancerType))
pancanAtlas_iCluster_methylation$color=factor(pancanAtlas_iCluster_methylation$color, levels=unique(pancanAtlas_iCluster_methylation$color))
write.table(pancanAtlas_iCluster_methylation, "Figure8ABCDE_S7B/Figure8C_data.txt", row.names = F, col.names = T, sep="\t", quote=F)

pancanAtlas_iCluster_methylation2=pancanAtlas_iCluster_methylation[order(pancanAtlas_iCluster_methylation$delta_ESCCDMRsVSEACDMRs),]
pancanAtlas_iCluster_methylation2_1=pancanAtlas_iCluster_methylation2[pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs<(-0.4)|pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs>0.4,]
pancanAtlas_iCluster_methylation2_2=pancanAtlas_iCluster_methylation2[(pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs<(-0.3)&pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs>(-0.4))|
                                                                        (pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs>0.3&pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs<0.4),]
pancanAtlas_iCluster_methylation2_3=pancanAtlas_iCluster_methylation2[(pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs<(-0.2)&pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs>(-0.3))|
                                                                        (pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs>0.2&pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs<0.3),]
pancanAtlas_iCluster_methylation2_4=pancanAtlas_iCluster_methylation2[(pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs>(-0.2)&pancanAtlas_iCluster_methylation2$delta_ESCCDMRsVSEACDMRs<0.2),]

p3=ggplot(pancanAtlas_iCluster_methylation2_4, aes(x=x, y=y, color=delta_ESCCDMRsVSEACDMRs)) + geom_point(alpha=0.7,size=0.8)
p3=p3+geom_point(data=pancanAtlas_iCluster_methylation2_3, aes(x=x, y=y, color=delta_ESCCDMRsVSEACDMRs), alpha=0.7,size=0.8)
p3=p3+geom_point(data=pancanAtlas_iCluster_methylation2_2, aes(x=x, y=y, color=delta_ESCCDMRsVSEACDMRs), alpha=0.7,size=0.8)
p3=p3+geom_point(data=pancanAtlas_iCluster_methylation2_1, aes(x=x, y=y, color=delta_ESCCDMRsVSEACDMRs), alpha=0.7,size=0.8)
p3=p3+theme_classic()+scale_colour_gradientn(colours=colorRampPalette(c("#FF0000", "#F12828", "#DC6464", "#C1B3B3", "grey", "grey", "grey", "#A0A0C8", "#6464DC", "#2828F1", "#0000FF"))(50),
                                             limits=c(-0.5, 0.5), name="delta methylation\n(ESCC_tumor_PMDs \nV.S. EAC_tumor_PMDs)")+theme_classic()
p3=p3+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p4=p3+theme(legend.position = "none")
png("Figure8ABCDE_S7B/Figure8C.png", res=300, width=1800, height=1800)
print(p4)
dev.off()
pdf("Figure8ABCDE_S7B/Figure8C.legend.pdf", width=6, height=6)
print(p3)
dev.off()


##############################################Figure8F-J##############################################
######Figure8F
cancerTypeInfo=read.table("Data/Figure8FHIJ_S7C/CancerType_info2.txt", header=T, sep="\t", stringsAsFactors = F, comment.char = "@")
domainResultFile="Data/Figure8FHIJ_S7C/ATAC_domain_result.txt"
domainResult=read.table(domainResultFile, sep="\t", stringsAsFactors = F, header=T)
domainResult[domainResult$cancerType%in%c("Basal","LumA", "LumB", "Her2"),]$cancerType="BRCA"
domainResult$cluster=gsub("_T[0-9]*", "", domainResult$cluster)
for(i in 2:5){
  domainResult[,i]=domainResult[,i]/domainResult$whole_genome
}
domainResult$delta_ESCCspecificVSEACspecific=domainResult$ESCC_specificPMDs-domainResult$EAC_specificPMDs
domainResult$delta_sharedHMD_PMD=domainResult$ESCA_sharedHMDs-domainResult$ESCA_sharedPMDs
domainResult=domainResult[,c(1,7:9)]
pancanAtlas_iCluster_ATAC=merge(pancanAtlas_iCluster, domainResult, by.x="barcode", by.y="cluster", all.x=T)
pancanAtlas_iCluster_ATAC1=pancanAtlas_iCluster_ATAC[is.na(pancanAtlas_iCluster_ATAC$cancerType),]
pancanAtlas_iCluster_ATAC2=pancanAtlas_iCluster_ATAC[!is.na(pancanAtlas_iCluster_ATAC$cancerType),]
pancanAtlas_iCluster_ATAC2[pancanAtlas_iCluster_ATAC2$cancerType%in%"BLCA (Ba/Sq)",]$color="green"
pancanAtlas_iCluster_ATAC2[pancanAtlas_iCluster_ATAC2$cancerType%in%"CESC (Sq)",]$color="#98FF98"
pancanAtlas_iCluster_ATAC2[pancanAtlas_iCluster_ATAC2$cancerType%in%"ESCC",]$color="#254117"
pancanAtlas_iCluster_ATAC2=pancanAtlas_iCluster_ATAC2[order(pancanAtlas_iCluster_ATAC2$cancerType),]
pancanAtlas_iCluster_ATAC2$cancerType=factor(pancanAtlas_iCluster_ATAC2$cancerType, levels=unique(pancanAtlas_iCluster_ATAC2$cancerType))
pancanAtlas_iCluster_ATAC2$color=factor(pancanAtlas_iCluster_ATAC2$color, levels=unique(pancanAtlas_iCluster_ATAC2$color))
pancanAtlas_iCluster_ATAC2$Type2="Others"
pancanAtlas_iCluster_ATAC2[pancanAtlas_iCluster_ATAC2$cancerType%in%c("BLCA (Ba/Sq)", "CESC (Sq)", "ESCC", "HNSC", "LUSC"),]$Type2="Pan-Squamous"
pancanAtlas_iCluster_ATAC2[pancanAtlas_iCluster_ATAC2$cancerType%in%c("EAC", "COAD", "STAD", "READ"),]$Type2="Pan-GI"
pancanAtlas_iCluster_ATAC2$Type2=factor(pancanAtlas_iCluster_ATAC2$Type2, levels=c("Others","Pan-Squamous","Pan-GI"))

write.table(pancanAtlas_iCluster_ATAC1, "Figure8FHIJ_S7C/Figure8FG_S7C.data.txt", row.names = F, col.names = T, sep="\t", quote=F)
write.table(pancanAtlas_iCluster_ATAC2, "Figure8FHIJ_S7C/Figure8FG_S7C.data.txt", row.names = F, col.names = F, sep="\t", quote=F, append = T)
p1=ggplot(pancanAtlas_iCluster_ATAC1, aes(x=x, y=y)) + geom_point(alpha=0.7,size=0.1, color="grey")
p1=p1+geom_point(data=pancanAtlas_iCluster_ATAC2, aes(x=x, y=y, color=Type2),size=2)
p1=p1+theme_classic()+scale_color_manual(values=c("grey", "#EA3323", "#0000F5"))
p1=p1+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p1=p1+theme(legend.position = "none")
png("Figure8FHIJ_S7C/Figure8F.png", res=300, width=2200, height=1800)
print(p1)
dev.off()

#####Figure8G
pancanAtlas_iCluster_ATAC3=pancanAtlas_iCluster_ATAC2
pancanAtlas_iCluster_ATAC3[pancanAtlas_iCluster_ATAC3$delta_ESCCspecificVSEACspecific>0.7,]$delta_ESCCspecificVSEACspecific=0.7
pancanAtlas_iCluster_ATAC3[pancanAtlas_iCluster_ATAC3$delta_ESCCspecificVSEACspecific<(-0.35),]$delta_ESCCspecificVSEACspecific=(-0.35)
p1=ggplot(pancanAtlas_iCluster_ATAC1, aes(x=x, y=y)) + geom_point(alpha=0.7,size=0.1, color="grey")
p1=p1+geom_point(data=pancanAtlas_iCluster_ATAC3, aes(x=x, y=y, color=delta_ESCCspecificVSEACspecific),size=2)
p1=p1+theme_classic()
p1=p1+scale_color_gradientn(name="delta accessiblity\n(ESCC only PMDs \nV.S. EAC only PMDs)",limits=c(-0.35,0.7),
                            colours = c("red", "white", "#AC82FF", "blue"))
p1=p1+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p2=p1+theme(legend.position = "none")
png("Figure8FHIJ_S7C/Figure8G.png", res=300, width=2200, height=1800)
print(p2)
dev.off()
pdf("Figure8FHIJ_S7C/Figure8G_withlegend.pdf", width=8, height=8)
print(p1)
dev.off()

#####FigureS7C
p1=ggplot(pancanAtlas_iCluster_ATAC1, aes(x=x, y=y)) + geom_point(alpha=0.7,size=0.1, color="grey")
p1=p1+geom_point(data=pancanAtlas_iCluster_ATAC3, aes(x=x, y=y, color=delta_sharedHMD_PMD),size=2)
p1=p1+theme_classic()
p1=p1+scale_color_gradientn(name="delta accessiblity\n(shared HMDs \nV.S. shared PMDs)",limits=c(-1.7,1.7),
                            colours = c("red", "white", "#AC82FF", "blue"))
p1=p1+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p2=p1+theme(legend.position = "none")
png("Figure8FHIJ_S7C/FigureS7C.png", res=300, width=2200, height=1800)
print(p2)
dev.off()
pdf("Figure8FHIJ_S7C/FigureS7C_withlegend.pdf", width=8, height=8)
print(p1)
dev.off()

######Figure8I
plotTumorATACPointPlot=function(plotdata2, ylabel, saveFile){
  plotdata2=plotdata2[order(plotdata2$cancerType, plotdata2$delta),]
  name<-as.character(unique(plotdata2$cancerType))
  sum=0
  for(i in 1:length(name)){
    tumorRow = nrow(plotdata2[plotdata2$cancerType==name[i],])
    for(j in 1:tumorRow){
      plotdata2[sum+j,"Postion"]<- (1000*(i-1)+350+250*j/tumorRow)
    }
    sum=sum+tumorRow
  } 
  for(i in 1:length(name)){
    plotdata2[plotdata2$cancerType==name[i],"Tick"]<-(1000*(i-1)+350)
  }
  meanData=plotdata2[1,]
  sum2=0
  for(i in 1:length(name)){
    tumorRow = nrow(plotdata2[plotdata2$cancerType==name[i],])
    temp=plotdata2[sum2+ceiling(tumorRow/2),]
    temp$delta=mean(plotdata2[plotdata2$cancerType==name[i],]$delta)
    meanData <-rbind(meanData,temp)
    sum2=sum2+tumorRow
  }
  meanData$xmin=meanData$Postion-125
  meanData$xmax=meanData$Postion+125
  meanData$ymin=meanData$delta
  meanData$ymax=meanData$delta
  meanData=meanData[-1,]
  write.table(plotdata2, gsub(".pdf", ".data.txt", saveFile), row.names = F, col.names = T, sep="\t", quote=F)
  p = ggplot(plotdata2, aes(x=Postion, y=delta, color=label))+geom_point(size=0.6,stat="identity")+
    geom_rect(data=meanData,aes(xmin = xmin, xmax = xmax, ymin = delta, ymax = delta),
              fill="white",size=0.4,color="orange")
  p=p+theme_classic()+ylab(ylabel)+ggtitle("tumor ATACseq")
  p=p+scale_color_manual(name="Cancer groups", values=c("#EA3323","black", "#0000F5"))
  p=p+scale_x_continuous(name = "", breaks = unique(plotdata2$Tick), labels = unique(plotdata2$cancerType))
  p=p+geom_hline(yintercept=0, linetype="dashed", color = "grey")
  p=p+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1), axis.text = element_text(color="black", size=12),
            plot.title = element_text(hjust = 0.5, size=16, face="bold"), legend.position = "bottom",
            legend.text = element_text(size=12))
  pdf(saveFile, width=9, height=4)
  print(p)
  dev.off()
}
plotdata=pancanAtlas_iCluster_ATAC2[,c(1,8,9,11)]
colnames(plotdata)=c("sample", "cancerType", "delta", "label")
medianData= do.call(c, lapply(unique(plotdata$cancerType), function(x){
  mean(plotdata[plotdata$cancerType%in%x,]$delta)
}))
medianData=data.frame(cancerType=unique(plotdata$cancerType), median=unique(medianData))
medianData=medianData[order(medianData$median),]
plotdata$cancerType=factor(plotdata$cancerType, levels=as.character(medianData$cancerType))
plotdata$label=factor(plotdata$label, levels=c("Pan-Squamous","Others","Pan-GI"))
plotTumorATACPointPlot(plotdata, "Detal accessibility (ESCC only PMDs V.S. EAC only PMDs)", "Figure8FHIJ_S7C/Figure8I.pdf")

######Figure8H
cancerTypeInfo=read.table("Data/Figure8FHIJ_S7C/CancerType_info2.txt", header=T, sep="\t", stringsAsFactors = F, comment.char = "@")
domainResultFile="Data/Figure8FHIJ_S7C/ATAC_domain_result.txt"
domainResult=read.table(domainResultFile, sep="\t", stringsAsFactors = F, header=T)

DMRResultFile="Data/Figure8FHIJ_S7C/ATAC_hypoDMR_result.txt"
DMRResult=read.table(DMRResultFile, sep="\t", stringsAsFactors = F, header=T)
DMRResult[grep("Basal|LumA|LumB|Her2", DMRResult$cancerType),]$cancerType="BRCA"

DMRResult=merge(DMRResult, domainResult[,colnames(domainResult)%in%c("cluster", "whole_genome")], by.x="cluster", by.y="cluster")
DMRResult$cluster= gsub("_.*", "", DMRResult$cluster)
for(i in 2:3){
  DMRResult[,i]=DMRResult[,i]/DMRResult$whole_genome
}
DMRResult=DMRResult[,-ncol(DMRResult)]
color_info=list(cancerType=as.character(cancerTypeInfo$Color))
names(color_info$cancerType)=as.character(cancerTypeInfo$CancerType)
DMRResult$delta=DMRResult$hypoESCC_DMRs-DMRResult$hypoEAC_DMRs
DMRResult=DMRResult[,c(1,4,5)]

pancanAtlas_iCluster_ATAC=merge(pancanAtlas_iCluster, DMRResult, by.x="barcode", by.y="cluster", all.x=T)
pancanAtlas_iCluster_ATAC1=pancanAtlas_iCluster_ATAC[is.na(pancanAtlas_iCluster_ATAC$cancerType),]
pancanAtlas_iCluster_ATAC2=pancanAtlas_iCluster_ATAC[!is.na(pancanAtlas_iCluster_ATAC$cancerType),]
pancanAtlas_iCluster_ATAC2[pancanAtlas_iCluster_ATAC2$cancerType%in%"BLCA (Ba/Sq)",]$color="green"
pancanAtlas_iCluster_ATAC2[pancanAtlas_iCluster_ATAC2$cancerType%in%"CESC (Sq)",]$color="#98FF98"
pancanAtlas_iCluster_ATAC2[pancanAtlas_iCluster_ATAC2$cancerType%in%"ESCC",]$color="#254117"
pancanAtlas_iCluster_ATAC2=pancanAtlas_iCluster_ATAC2[order(pancanAtlas_iCluster_ATAC2$cancerType),]
pancanAtlas_iCluster_ATAC2$cancerType=factor(pancanAtlas_iCluster_ATAC2$cancerType, levels=unique(pancanAtlas_iCluster_ATAC2$cancerType))
pancanAtlas_iCluster_ATAC2$color=factor(pancanAtlas_iCluster_ATAC2$color, levels=unique(pancanAtlas_iCluster_ATAC2$color))
pancanAtlas_iCluster_ATAC2$Type2="Others"
pancanAtlas_iCluster_ATAC2[pancanAtlas_iCluster_ATAC2$cancerType%in%c("BLCA (Ba/Sq)", "CESC (Sq)", "ESCC", "HNSC", "LUSC"),]$Type2="Pan-Squamous"
pancanAtlas_iCluster_ATAC2[pancanAtlas_iCluster_ATAC2$cancerType%in%c("EAC", "COAD", "STAD", "READ"),]$Type2="Pan-GI"
pancanAtlas_iCluster_ATAC3=pancanAtlas_iCluster_ATAC2
pancanAtlas_iCluster_ATAC3[pancanAtlas_iCluster_ATAC3$delta>7,]$delta=7
pancanAtlas_iCluster_ATAC3[pancanAtlas_iCluster_ATAC3$delta<(-7),]$delta=(-7)
write.table(pancanAtlas_iCluster_ATAC1, "Figure8FHIJ_S7C/Figure8H.data.txt", row.names = F, col.names = T, sep="\t", quote=F)
write.table(pancanAtlas_iCluster_ATAC2, "Figure8FHIJ_S7C/Figure8H.data.txt", row.names = F, col.names = F, sep="\t", quote=F, append = T)
p1=ggplot(pancanAtlas_iCluster_ATAC1, aes(x=x, y=y)) + geom_point(alpha=0.7,size=0.1, color="grey")
p1=p1+geom_point(data=pancanAtlas_iCluster_ATAC3, aes(x=x, y=y, color=delta),size=2)
p1=p1+theme_classic()
p1=p1+scale_color_gradient2(name="delta accessiblity\n(ESCC tumor hypoDMRs \nV.S. EAC tumor hypoDMRs)",limits=c(-7,7),
                            low ="blue", mid = "white", high ="red", guide = "colourbar",midpoint = 0)
p1=p1+theme(axis.line=element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
p2=p1+theme(legend.position = "none")
png("Figure8FHIJ_S7C/Figure8H.png", res=300, width=2200, height=1800)
print(p2)
dev.off()
pdf("Figure8FHIJ_S7C/Figure8H.withlegend.pdf", width=8, height=7)
print(p1)
dev.off()

######Figure8J
plotdata=pancanAtlas_iCluster_ATAC2[,c(1,8,9,10)]
colnames(plotdata)=c("sample", "cancerType", "delta", "label")
medianData= do.call(c, lapply(unique(plotdata$cancerType), function(x){
  mean(plotdata[plotdata$cancerType%in%x,]$delta)
}))
medianData=data.frame(cancerType=unique(plotdata$cancerType), median=unique(medianData))
medianData=medianData[order(medianData$median,decreasing = T),]
plotdata$cancerType=factor(plotdata$cancerType, levels=as.character(medianData$cancerType))
plotdata$label=factor(plotdata$label, levels=c("Pan-Squamous","Others","Pan-GI"))
plotTumorATACPointPlot(plotdata, "Detal accessibility (ESCC tumor hypoDMRs V.S. EAC/GEJ tumor hypoDMRs)", "Figure8FHIJ_S7C/Figure8J.pdf")

##############################################FigureS7D and S7E##############################################
#######FigureS7D
load("Data/FigureS7D/GSE72874_HM450k.RData")
normal_samples=sampleInfo[sampleInfo$Type%in%"Normal",]
GERD_samples=sampleInfo[sampleInfo$Type%in%"GERD",]
BE_samples=sampleInfo[sampleInfo$Type%in%"BE",]
EAC_samples=sampleInfo[sampleInfo$Type%in%"Tumour",]
targetSamples=rbind(normal_samples, GERD_samples, BE_samples, EAC_samples)

getMethylationHT450k=function(regionFile, probesFile, type, methylationMatrix){
  tempData=read_bed(regionFile,n_fields = 5)
  probesMatrix=read_bed(probesFile, n_fields = 4)
  tempData_probes=bed_intersect(probesMatrix,tempData)
  tempData_probes=as.data.frame(tempData_probes)
  tempData_probes=unique(tempData_probes[tempData_probes$.overlap==2,]$name.x)
  print(length(tempData_probes))
  tempData=methylationMatrix[rownames(methylationMatrix)%in%tempData_probes,, drop=F]
  tempData=data.frame(colMeans(tempData, na.rm = T))
  colnames(tempData)=type
  return(tempData)
}
plotHM450kPMDMethylationLinePlot=function(probesFile, methylationMatrix, sampleInfo){
  temp1=getMethylationHT450k("Data/MMSeekR_PMDs/ESCC_specificPMDs.bed", probesFile, "ESCC_specificPMDs",  methylationMatrix)
  temp2=getMethylationHT450k("Data/MMSeekR_PMDs/EAC_specificPMDs.bed", probesFile, "EAC_specificPMDs", methylationMatrix)
  temp3=getMethylationHT450k("Data/MMSeekR_PMDs/ESCA_sharedPMDs.bed", probesFile, "SharedPMDs", methylationMatrix)
  temp4=getMethylationHT450k("Data/MMSeekR_PMDs/ESCA_sharedHMDs.bed", probesFile, "SharedHMDs", methylationMatrix)
  plotdata=cbind(temp1, temp2, temp3, temp4)
  plotdata=plotdata[rownames(plotdata)%in%as.character(sampleInfo$Sample),]
  plotdata=data.frame(Sample=rownames(plotdata), plotdata)
  plotdata=merge(plotdata, sampleInfo, by.x="Sample", by.y="Sample")
  return(plotdata)
}
plotHM450kDMRMethylationLinePlot=function(probesFile, methylationMatrix, sampleInfo){
  temp1=getMethylationHT450k("Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed", probesFile, "hypoEAC",  methylationMatrix)
  temp2=getMethylationHT450k("Data/MaskUnionPMDs_DMRs/hypoESCC_Tumor.bed", probesFile, "hypoESCC", methylationMatrix)
  plotdata=cbind(temp1, temp2)
  plotdata=plotdata[rownames(plotdata)%in%as.character(sampleInfo$Sample),]
  plotdata=data.frame(Sample=rownames(plotdata), plotdata)
  plotdata=merge(plotdata, sampleInfo, by.x="Sample", by.y="Sample")
  return(plotdata)
}

GSE72874_PMD_methylation=plotHM450kPMDMethylationLinePlot("meta/HT450k.probe.rmCGI.rmblackList.solo.bed", methylationMatrix, targetSamples)
GSE72874_PMD_methylation$Source="GSE72874"
write.table(GSE72874_PMD_methylation, "FigureS7D/FigureS7D.PMD.txt", row.names = F, col.names = T, sep="\t", quote=F)
GSE72874_DMR_methylation=plotHM450kDMRMethylationLinePlot("meta/HT450k.probe.rmblackList.bed", methylationMatrix, targetSamples)
GSE72874_DMR_methylation$Source="GSE72874" 
write.table(GSE72874_DMR_methylation, "FigureS7D/FigureS7D.DMR.txt", row.names = F, col.names = T, sep="\t", quote=F)

plotPMDMethylationHeatmapAndLinePlot=function(plotdata, targetTypes, saveFile, height, width){
  my_theme=theme_classic()+ theme(plot.title = element_text(hjust = 0.5, size=14, color="black", face="bold"),
                                  axis.text.x = element_text(size=11, color="black", angle = 90, hjust = 1, vjust = 0.5),
                                  axis.text.y = element_text(size=11, color="black", angle = 0),
                                  axis.title = element_text(size=13, color="black", face="bold"),
                                  legend.title =element_text(size=12, color="black", face="bold"),
                                  legend.text =element_text(size=12, color="black"),legend.position = "bottom")
  plotList=list()
  plotdata2=do.call(rbind, lapply(unique(plotdata$Type), function(x){
    a=data.frame(mean=colMeans(plotdata[plotdata$Type%in%x,2:3]), sd=colSds(as.matrix(plotdata[plotdata$Type%in%x,2:3])),
                 sampleType=x)
    a$type=rownames(a)
    return(a)
  }))
  plotdata2$type=factor(plotdata2$type, levels=c("EAC_specificPMDs", "ESCC_specificPMDs"))
  plotdata2=plotdata2[plotdata2$sampleType%in%targetTypes,]
  plotdata2$sampleType=factor(plotdata2$sampleType,levels=targetTypes)
  ymin=0
  ymax=0.8
  p1<- ggplot(plotdata2, aes(x=type, y=mean, group=sampleType,color=sampleType)) + geom_point(shape=15) + geom_line() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position=position_dodge(0))+ ylim(ymin, ymax)+my_theme
  p1=p1+scale_color_manual(values=c("orange", "blue"))
  plotList[["trend"]]=p1
  
  
  for(type in targetTypes){
    data=plotdata[plotdata$Type%in%type,]
    rownames(data)=data$Sample
    data=data[,c(3,2)]
    data=data[order(rowSums(data)),]
    temp=data.frame(type=rep(type, nrow(data)))
    rownames(temp)=rownames(data)
    if(type%in%targetTypes[1]){
      plotdata3=data
      sampleList=temp
      annotate_col=c("orange")
    }else{
      plotdata3=rbind(plotdata3, data)
      sampleList=rbind(sampleList, temp)
      annotate_col=c(annotate_col, "blue")
    }
  }
  names(annotate_col)=targetTypes
  annotate_col=list(type=annotate_col)
  breaks=seq(0,0.8,by=0.001)
  
  p1=pheatmap(as.matrix(plotdata3), show_rownames = F, cluster_rows = F, cluster_cols = F, annotation_row = sampleList,
              breaks=breaks, border_color = NA, annotation_colors = annotate_col, 
              color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(length(breaks)))
  plotList[["heatmap"]]=p1[[4]]
  
  pdf(saveFile, height = height, width=width)
  print(ggarrange(plotlist = plotList))
  dev.off()
}
plotPMDMethylationHeatmapAndLinePlot(GSE72874_PMD_methylation, c("BE", "Tumour"), "FigureS7D/FigureS7D_PMD_heatmapLine.pdf", 5,10)

plotDMRMethylationHeatmapAndLinePlot=function(plotdata, targetTypes, saveFile, height, width){
  my_theme=theme_classic()+ theme(plot.title = element_text(hjust = 0.5, size=14, color="black", face="bold"),
                                  axis.text.x = element_text(size=11, color="black", angle = 90, hjust = 1, vjust = 0.5),
                                  axis.text.y = element_text(size=11, color="black", angle = 0),
                                  axis.title = element_text(size=13, color="black", face="bold"),
                                  legend.title =element_text(size=12, color="black", face="bold"),
                                  legend.text =element_text(size=12, color="black"),legend.position = "bottom")
  plotList=list()
  plotdata2=do.call(rbind, lapply(unique(plotdata$Type), function(x){
    a=data.frame(mean=colMeans(plotdata[plotdata$Type%in%x,2:3]), sd=colSds(as.matrix(plotdata[plotdata$Type%in%x,2:3])),
                 sampleType=x)
    a$type=rownames(a)
    return(a)
  }))
  plotdata2$type=factor(plotdata2$type, levels=c("hypoEAC", "hypoESCC"))
  plotdata2=plotdata2[plotdata2$sampleType%in%targetTypes,]
  plotdata2$sampleType=factor(plotdata2$sampleType,levels=targetTypes)
  ymin=0
  ymax=0.8
  p1<- ggplot(plotdata2, aes(x=type, y=mean, group=sampleType,color=sampleType)) + geom_point(shape=15) + geom_line() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, position=position_dodge(0))+ ylim(ymin, ymax)+my_theme
  p1=p1+scale_color_manual(values=c("orange", "blue"))
  plotList[["trend"]]=p1
  
  
  for(type in targetTypes){
    data=plotdata[plotdata$Type%in%type,]
    rownames(data)=data$Sample
    data=data[,c(2:3)]
    data=data[order(rowSums(data)),]
    temp=data.frame(type=rep(type, nrow(data)))
    rownames(temp)=rownames(data)
    if(type%in%targetTypes[1]){
      plotdata3=data
      sampleList=temp
      annotate_col=c("orange")
    }else{
      plotdata3=rbind(plotdata3, data)
      sampleList=rbind(sampleList, temp)
      annotate_col=c(annotate_col, "blue")
    }
  }
  names(annotate_col)=targetTypes
  annotate_col=list(type=annotate_col)
  breaks=seq(0,0.8,by=0.001)
  
  p1=pheatmap(as.matrix(plotdata3), show_rownames = F, cluster_rows = F, cluster_cols = F, annotation_row = sampleList,
              breaks=breaks, border_color = NA, annotation_colors = annotate_col, 
              color = colorRampPalette(c("#000436","#021EA9","#1632FB","#6E34FC","#C732D5","#FD619D","#FF9965","#FFD32B","#FFFC5A"))(length(breaks)))
  plotList[["heatmap"]]=p1[[4]]
  
  pdf(saveFile, height = height, width=width)
  print(ggarrange(plotlist = plotList))
  dev.off()
}
plotDMRMethylationHeatmapAndLinePlot(GSE72874_DMR_methylation, c("BE", "Tumour"), "FigureS7D/Figure8K_DMR_heatmapLine.pdf", 5,10)

#######Figure S7E
load("Data/FigureS7E/GSE81334_HM450k.RData")
normal_samples=sampleInfo[sampleInfo$Type%in%"squamous",]
BE_samples=sampleInfo[sampleInfo$Type%in%c("BE", "EAC.BE"),]
EAC_samples=sampleInfo[sampleInfo$Type%in%"EAC",]
targetSamples=rbind(normal_samples, BE_samples, EAC_samples)
HM450k=read_bed("meta/HT450k.probe.bed",n_fields = 4)
GSE81334_PMD_methylation=plotHM450kPMDMethylationLinePlot("meta/HT450k.probe.rmCGI.rmblackList.solo.bed", methyMatrix, targetSamples)
GSE81334_PMD_methylation$Source="GSE81334"
GSE81334_PMD_methylation[GSE81334_PMD_methylation$Type%in%"EAC.BE",]$Type="BE"
write.table(GSE81334_PMD_methylation, "FigureS7E/FigureS7E.PMD.txt", row.names = F, col.names = T, sep="\t", quote=F)
GSE81334_DMR_methylation=plotHM450kDMRMethylationLinePlot("meta/HT450k.probe.rmblackList.bed", methyMatrix, targetSamples)
GSE81334_DMR_methylation$Source="GSE81334" 
GSE81334_DMR_methylation[GSE81334_DMR_methylation$Type%in%"EAC.BE",]$Type="BE"
write.table(GSE81334_DMR_methylation, "FigureS7E/FigureS7E.DMR.txt", row.names = F, col.names = T, sep="\t", quote=F)
plotPMDMethylationHeatmapAndLinePlot(GSE81334_PMD_methylation, c("BE", "EAC"), "FigureS7E/FigureS7E_PMD_heatmapLine.pdf", 5,10)
plotDMRMethylationHeatmapAndLinePlot(GSE81334_DMR_methylation, c("BE", "EAC"), "FigureS7E/FigureS7E_DMR_heatmapLine.pdf", 5,10)


##############################################Figure8KLM_S7F##############################################
load("Data/Figure8KLM_S7F//TCGA_methylation.full.Rdata")
tumorResultMatrix2$Label="Other cancers"
tumorResultMatrix2[tumorResultMatrix2$cancerType%in%c("ESCC", "HNSC", "LUSC", "CESC (Sq)", "BLCA (Ba/Sq)"),]$Label="Squamous cancers"
tumorResultMatrix2[tumorResultMatrix2$cancerType%in%c("EAC", "COAD", "STAD", "READ"),]$Label="GI cancers"
tumorResultMatrix2$Label=factor(tumorResultMatrix2$Label, levels=c("GI cancers", "Other cancers", "Squamous cancers"))
PMD_tumorResultMatrix=tumorResultMatrix2[,c(1:3,7,6,9)]

load("Data/Figure8KLM_S7F//TCGA_methylation_AllDMR.all.Rdata")
tumorResultMatrix2$Label="Other cancers"
tumorResultMatrix2[tumorResultMatrix2$cancerType%in%c("ESCC", "HNSC", "LUSC", "CESC (Sq)", "BLCA (Ba/Sq)"),]$Label="Squamous cancers"
tumorResultMatrix2[tumorResultMatrix2$cancerType%in%c("EAC", "COAD", "STAD", "READ"),]$Label="GI cancers"
tumorResultMatrix2$Label=factor(tumorResultMatrix2$Label, levels=c("GI cancers", "Other cancers", "Squamous cancers"))
DMR_tumorResultMatrix=tumorResultMatrix2[,c(1:3,6,4,5)]

plotdata=merge(PMD_tumorResultMatrix, DMR_tumorResultMatrix, by=c("SampleName", "cancerType", "Label"))
save(plotdata, file="Figure8KLM_S7F/TCGA_specific_PMDsDMRs_forPredict.v2.RData")

library(pROC)
library(ROCR)
library(multiROC)
load("Figure8KLM_S7F/TCGA_specific_PMDsDMRs_forPredict.v2.RData")
predictData=plotdata[,c(3:9)]
rownames(predictData)=plotdata$SampleName
predictData$Label=factor(predictData$Label, levels=c("Other cancers", "Squamous cancers", "GI cancers"))

calculateMultiROC_change=function(predictData, saveIndex, features){
  # myMultiPredict_round_list=list()
  # for(i in 1:100){
  #   print(paste0("Round", i))
  #   myMultiPredict_round_list[[i]]=myMultiPredict(predictData, features, paste0(gsub(basename(saveIndex), "", saveIndex), "/Temp/", basename(saveIndex), "_downsampling.round", i, ".RData"))
  # }
  # save(myMultiPredict_round_list, file=paste0(saveIndex, "_downsampling_result.RData"))
  load(paste0(saveIndex, "_downsampling_result.RData"))
  for(i in 1:100){
    tmp=myMultiPredict_round_list[[i]]
    rownames(tmp)=paste0(rownames(tmp), "_", i)
    if(i==1){
      my_prediction=tmp
    }else{
      my_prediction=rbind(my_prediction, tmp)
    }
  }
  
  AUC_list=list()
  res <- multi_roc(my_prediction, force_diag=T)
  plot_roc_df <- plot_roc_data(res)
  plot_roc_df=plot_roc_df[plot_roc_df$Group%in%c("Squamous", "GI", "Other"),]
  plot_roc_df$Group=factor(plot_roc_df$Group, levels=c("Other", "Squamous", "GI"))
  p=ggplot(plot_roc_df, aes(x = 1-Specificity, y=Sensitivity)) +
    geom_path(aes(color = Group), size=1) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour='grey', linetype = 'dotdash') +
    scale_color_manual(values=c("darkgrey", "red", "blue"))+
    theme_classic()+theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
                          axis.ticks.length = unit(.18, "cm"), 
                          axis.ticks = element_line(colour = "black"),
                          plot.title = element_text(hjust = 0.5), 
                          legend.position = "none")
  
  png(paste0(saveIndex, "_combined.png"), width=1200, height=1000, res = 300)
  print(p)
  dev.off()
  AUC_list[["roc"]]=res$AUC
  
  pr_res <- multi_pr(my_prediction, force_diag=T)
  plot_pr_df <- plot_pr_data(pr_res)
  plot_pr_df=plot_pr_df[plot_pr_df$Group%in%c("Squamous", "GI", "Other"),]
  plot_pr_df$Group=factor(plot_pr_df$Group, levels=c("Other", "Squamous", "GI"))
  
  p1=ggplot(plot_pr_df, aes(x=Recall, y=Precision)) + 
    geom_path(aes(color = Group), size=1) +
    geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), colour='grey', linetype = 'dotdash') +
    theme_classic()+theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
                          axis.ticks.length = unit(.18, "cm"), 
                          axis.ticks = element_line(colour = "black"),
                          plot.title = element_text(hjust = 0.5), 
                          legend.position = "none")+ 
    scale_color_manual(values=c("darkgrey", "red", "blue"))
  png(paste0(saveIndex, "_combined.PR.png"), width=1200, height=1000, res = 300)
  print(p1)
  dev.off()
  AUC_list[["pr"]]=pr_res$AUC
  return(AUC_list)
}
allList=calculateMultiROC_change(predictData, "Figure8KLM_S7F/ROC-multiLogic_PMDDMR", c("ESCC_specificPMDs", "EAC_specificPMDs", "ESCC_DMRs", "EAC_DMRs"))
PMDList=calculateMultiROC_change(predictData, "Figure8KLM_S7F/ROC-multiLogic_PMD", c("ESCC_specificPMDs", "EAC_specificPMDs"))
DMRList=calculateMultiROC_change(predictData, "Figure8KLM_S7F/ROC-multiLogic_DMR", c("ESCC_DMRs", "EAC_DMRs"))
