library(dplyr)
library(readr)
library(valr)

args<-commandArgs(TRUE)
fragmentFile=args[1]
fragmentIndex=gsub("_.*", "", basename(fragmentFile))
print(fragmentIndex)
sampleList=read.table("Data/Figure7EFG_S5C/sampleInfo.txt", header=T, sep="\t")
targetTissueType=read.table("Data/Figure7EFG_S5C/GSE184462_tissueType.txt", header=T, sep="\t")
targetEpithelialCells=targetTissueType[targetTissueType$Cells%in%"Epithelial cells",]$Tissue_contribution
metaInfo=read.table("Data/Figure7EFG_S5C/GSE184462_metadata.tsv", sep="\t", stringsAsFactors = F, header=T)
targetMetaInfo=metaInfo[metaInfo$cell.type%in%targetEpithelialCells,]
targetMetaInfo$barcodeID=gsub(".*\\+", "", targetMetaInfo$cellID)
targetMetaInfo=targetMetaInfo[targetMetaInfo$sample%in%sampleList[sampleList$Gse_Link%in%fragmentIndex,]$Index,]

readDomian=function(inputFile){
  data=read.table(inputFile, sep="\t", stringsAsFactors =F)
  data=tibble(chrom=data[,1], start=as.numeric(data[,2]), end=as.numeric(data[,3]))
  return(data)
}
EAC_DMRs=readDomian("Data/MaskUnionPMDs_DMRs/hypoEAC_Tumor.bed")
ESCC_DMRs=readDomian("Data/MaskUnionPMDs_DMRs/hypoESCC_Tumor.bed")

genomeData=read_bed("Figure7EFG_S5C/hg38_chr_5000.bed")
genomeData$ID=paste0("ID", 1:nrow(genomeData))
fragments=read_tsv(fragmentFile, col_names = FALSE)


targetFragments=fragments[fragments$X4%in%targetMetaInfo$barcodeID,]
if(nrow(targetFragments)>0){
  targetFragments=tibble(chrom=targetFragments$X1, start=targetFragments$X2, end=targetFragments$X3,
                         barcordID=targetFragments$X4, number=targetFragments$X5)
  targetFragments=targetFragments[targetFragments$chrom%in%paste0("chr", c(1:22,"X", "Y", "M")),]
  targetFragments$length=targetFragments$end-targetFragments$start
  targetFragments=bed_sort(targetFragments)
  result=bed_intersect(genomeData, targetFragments)
  result=result[result$.overlap>=result$length.y*0.5, c(1:4,7)]
  colnames(result)=c("chrom", "start", "end", "ID", "barcode")
  backgroup_result=data.frame(table(result$barcode))
  colnames(backgroup_result)=c("barcode", "total")
  backgroup_result$total=backgroup_result$total/nrow(genomeData)

  getDomainCounts=function(domain, name){
    tmp=bed_intersect(genomeData, domain)
    tmp=unique(tmp[tmp$.overlap>5000*0.5,1:4])
    colnames(tmp)=c("chrom", "start", "end", "ID")
    tmpData=result[result$ID%in%tmp$ID,]
    if(nrow(tmpData)>0){
      tmpData_Result=data.frame(table(tmpData$barcode))
       colnames(tmpData_Result)=c("barcode", name)
       tmpData_Result[[name]]=tmpData_Result[[name]]/nrow(tmp)
    }else{
       tmpData_Result=data.frame(barcode= character(), name= numeric())
       colnames(tmpData_Result)=c("barcode", name)
    }
    return(tmpData_Result)
  }
  EAC_DMRs_result=getDomainCounts(EAC_DMRs, "EAC_DMRs")
  ESCC_DMRs_result=getDomainCounts(ESCC_DMRs, "ESCC_DMRs")


  merge_result=merge(backgroup_result, EAC_DMRs_result, by="barcode", all.x=T)
  merge_result=merge(merge_result, ESCC_DMRs_result, by="barcode", all.x=T)
  rownames(merge_result)=merge_result$barcode
  merge_result=merge_result[,-1]
  merge_result[is.na(merge_result)]=0
  merge_result$EAC_DMRs=merge_result$EAC_DMRs/merge_result$total
  merge_result$ESCC_DMRs=merge_result$ESCC_DMRs/merge_result$total
  merge_result=data.frame(barcodeID=rownames(merge_result), merge_result)
  merge_result=merge(targetMetaInfo[,c("tissue", "cell.type", "barcodeID")], merge_result, by="barcodeID")
  merge_result$sample=fragmentIndex
  write.table(merge_result, file=paste0("Data/Figure7EFG_S5C/DMRresult/", fragmentIndex, "_result.txt"), row.names = F, col.names = T, quote=F, sep="\t")
}
