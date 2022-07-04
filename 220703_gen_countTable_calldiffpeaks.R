#!/usr/bin/env Rscript
options(echo=T)
args = commandArgs(trailingOnly=TRUE)
library(edgeR)
library(stringr)
library(data.table)

output_dir<-args[1]
print(output_dir)

coverage_files<-list.files(output_dir,full.names = TRUE)
coverage_files_name<-str_split_fixed(
  list.files(output_dir,full.names = FALSE),
  "_ref_peaks_coverage.bed",2)[,1]
coverage_file<-fread(coverage_files[1])
names<-paste(coverage_file$V1,coverage_file$V2,coverage_file$V3,sep="_")
counts<-data.frame(names)
for(i in 1:length(coverage_files)) {
  coverage_file<-fread(coverage_files[i])
  counts[,i+1]<-coverage_file$V11
  names(counts)[i+1]<-coverage_files_name[i]
}
head(counts)
# 
# 
# counts<-counts[-(grep("JH", counts$name)), ] #Remove unscaffolded 
# counts<-counts[-(grep("chrUn", counts$name)), ]
# counts<-counts[-(grep("GL", counts$name)), ]
# group<-c("c","c","d","d")
# 
# #prior<-fread("~/Box/GLIS3_ATAC/summary/test_project_mm10_peaks_coverage.tsv")
# 
# 
# y<-DGEList(counts = counts[,c(2:5)],group=group,remove.zeros = TRUE,genes = counts[,c("names")])
# keep<-filterByExpr(y)
# y<-y[keep, , keep.lib.sizes=FALSE]
# y<-calcNormFactors((y))
# y<-estimateDisp(y)
# et<-exactTest(y)
# topTags(et)
# glis3DE<-(et$table)
# glis3DE$region<-et$genes$genes
# glis3DE$FDR<-p.adjust(glis3DE$PValue,method = "BH")
# glis3DE<-glis3DE[,c("region","logCPM","logFC","PValue","FDR")]
# 
# nrow(glis3DE[glis3DE$FDR<0.05 & glis3DE$logFC>0,])
# nrow(glis3DE[glis3DE$FDR<0.05 & glis3DE$logFC<0,])
# 
# volcanoData <- cbind(glis3DE$logFC, -log10(glis3DE$FDR))
# colnames(volcanoData) <- c("logFC", "-LogPval")
# DEGs <- glis3DE$FDR < 0.05 & glis3DE$logFC>0
# point.col <- ifelse(DEGs, "red", "black")
# plot(volcanoData, pch = 16, col = point.col, cex = 0.5, xlim=c(-5,5))
# 
# 
# glis3DE_up<-glis3DE[glis3DE$logFC>0 & glis3DE$FDR<0.05,]
# glis3DE_up_regions<-str_split_fixed(glis3DE_up$region,"_",3)
# glis3DE_up$chr<-glis3DE_up_regions[,1]
# glis3DE_up$start<-as.integer(glis3DE_up_regions[,2])
# glis3DE_up$end<-as.integer(glis3DE_up_regions[,3])
# glis3DE_up_bed<-glis3DE_up[,c("chr","start","end")]
# write.table(glis3DE_up_bed, "~/Box/GLIS3_ATAC/GLIS3_nolambda_qe-7_sh-30_peaks/homer_input/glis3DE_up_0.05.bed",quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")
# 
# glis3DE_up_01<-glis3DE[glis3DE$logFC>0 & glis3DE$FDR<0.01,]
# glis3DE_up_01_regions<-str_split_fixed(glis3DE_up_01$region,"_",3)
# glis3DE_up_01$chr<-glis3DE_up_01_regions[,1]
# glis3DE_up_01$start<-as.integer(glis3DE_up_01_regions[,2])
# glis3DE_up_01$end<-as.integer(glis3DE_up_01_regions[,3])
# glis3DE_up_01_bed<-glis3DE_up_01[,c("chr","start","end")]
# write.table(glis3DE_up_01_bed, "~/Box/GLIS3_ATAC/GLIS3_nolambda_qe-7_sh-30_peaks/homer_input/glis3DE_up_01.bed",quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")
# 
# glis3DE_up_001<-glis3DE[glis3DE$logFC>0 & glis3DE$FDR<0.001,]
# glis3DE_up_001_regions<-str_split_fixed(glis3DE_up_001$region,"_",3)
# glis3DE_up_001$chr<-glis3DE_up_001_regions[,1]
# glis3DE_up_001$start<-as.integer(glis3DE_up_001_regions[,2])
# glis3DE_up_001$end<-as.integer(glis3DE_up_001_regions[,3])
# glis3DE_up_001_bed<-glis3DE_up_001[,c("chr","start","end")]
# write.table(glis3DE_up_001_bed, "~/Box/GLIS3_ATAC/GLIS3_nolambda_qe-7_sh-30_peaks/homer_input/glis3DE_up_001.bed",quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")
# 
# 
# glis3DE_up_00001<-glis3DE[glis3DE$logFC>0 & glis3DE$FDR<0.00001,]
# glis3DE_up_00001_regions<-str_split_fixed(glis3DE_up_00001$region,"_",3)
# glis3DE_up_00001$chr<-glis3DE_up_00001_regions[,1]
# glis3DE_up_00001$start<-as.integer(glis3DE_up_00001_regions[,2])
# glis3DE_up_00001$end<-as.integer(glis3DE_up_00001_regions[,3])
# glis3DE_up_00001_bed<-glis3DE_up_00001[,c("chr","start","end")]
# write.table(glis3DE_up_00001_bed, "~/Box/GLIS3_ATAC/GLIS3_nolambda_qe-7_sh-30_peaks/homer_input/glis3DE_up_00001.bed",quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")
# 
# 
# glis3DE_down<-glis3DE[glis3DE$logFC<0 & glis3DE$FDR<0.05,]
# glis3DE_down_regions<-str_split_fixed(glis3DE_down$region,"_",3)
# glis3DE_down$chr<-glis3DE_down_regions[,1]
# glis3DE_down$start<-as.integer(glis3DE_down_regions[,2])
# glis3DE_down$end<-as.integer(glis3DE_down_regions[,3])
# glis3DE_down_bed<-glis3DE_down[,c("chr","start","end")]
# write.table(glis3DE_down_bed, "~/Box/GLIS3_ATAC/GLIS3_nolambda_qe-7_sh-30_peaks/glis3DE_down_0.05.bed",quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")
# 
# glis3DE_regions<-str_split_fixed(glis3DE$region,"_",3)
# glis3DE$chr<-glis3DE_regions[,1]
# glis3DE$start<-as.integer(glis3DE_regions[,2])
# glis3DE$end<-as.integer(glis3DE_regions[,3])
# glis3DE_bed<-glis3DE[,c("chr","start","end")]
# write.table(glis3DE_bed, "~/Box/GLIS3_ATAC/GLIS3_nolambda_qe-7_sh-30_peaks/homer_input/glis3DE_background.bed",quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")