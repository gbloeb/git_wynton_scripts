#!/usr/bin/env Rscript
options(echo=T)
args = commandArgs(trailingOnly=TRUE)
library(edgeR)
library(stringr)
library(data.table)

#counts<-read.csv("~/count_table.csv")
#args<-c("trial","ML3_rep2_Cont1_200k_S10","ML3_rep2_Cont2_200k_S11","ML3_rep2_Cont3_200k_S12")

output_dir<-args[1]

setwd(output_dir)
dir.create("homer_input")
dir.create("220703_gen_countTable_output")



if(length(args)>1){                 #Not yet fully implemented, would allow indicating the control samples in the call of Rscript
ctrl_samples<-args[2:length(args)]
}
print(output_dir)

coverage_files<-setdiff(list.files(), list.dirs(recursive = FALSE, full.names = FALSE))
coverage_files_name<-str_split_fixed(coverage_files, "_coverage",2)[,1]
coverage_file<-fread(coverage_files[1])
names<-paste(coverage_file$V1,coverage_file$V2,coverage_file$V3,sep="_")
counts<-data.frame(names)
for(i in 1:length(coverage_files)) {
  coverage_file<-fread(coverage_files[i])
  counts[,i+1]<-coverage_file$V11
  names(counts)[i+1]<-coverage_files_name[i]
}
head(counts)

write.csv(counts, paste0("220703_gen_countTable_output/","count_table.csv"), row.names = FALSE)

counts<-counts[-(grep("JH", counts$name)), ] #Remove unscaffolded 
counts<-counts[-(grep("chrUn", counts$name)), ]
counts<-counts[-(grep("GL", counts$name)), ]


if(length(args)==1) {  #If no control samples are given, will attempt to search for ones that have "ctrl" in name

  ctrl_index<-grep("ctrl", colnames(counts)[2:ncol(counts)], ignore.case=TRUE) #only works if control has ctrl in file names
  if(length(ctrl_index)<2){stop(paste("These columns of count matrix contain control please rerun with contol samples indicated:",ctrl_index)) }

} else {    #else will use control sample names to select controls
  ctrl_index<-match(ctrl_samples, colnames(counts)[2:ncol(counts)]) 
}

group<-vector(mode = "character",length = ncol(counts)-1)
group[ctrl_index]<-"c" #Mark controls 
group[-ctrl_index]<-"e"


print("Ctrl samples are:")
print(colnames(counts)[2:ncol(counts)][ctrl_index])
print("Experimental samples are:")
print(colnames(counts)[2:ncol(counts)][-ctrl_index])

design_matrix<-data.frame(colnames(counts)[2:ncol(counts)],group)
#prior<-fread("~/Box/GLIS3_ATAC/summary/test_project_mm10_peaks_coverage.tsv")
write.csv(design_matrix,file = paste0("220703_gen_countTable_output/","design_matrix.csv"), row.names = FALSE)

y<-DGEList(counts = counts[,c(2:ncol(counts))],group=group,remove.zeros = TRUE,genes = counts[,c("names")])
keep<-filterByExpr(y)

pdf(paste0("220703_gen_countTable_output/","MDSplot_group.pdf"))
plotMDS(y,labels = group, gene.selection = "common")
dev.off()

pdf(paste0("220703_gen_countTable_output/","MDSplot_sample.pdf"))
plotMDS(y, gene.selection = "common")
dev.off()

y<-y[keep, , keep.lib.sizes=FALSE]
y<-calcNormFactors((y))
y<-estimateDisp(y)
et<-exactTest(y)
topTags(et)
TestedPeaks<-(et$table)
TestedPeaks$region<-et$genes$genes
TestedPeaks$FDR<-p.adjust(TestedPeaks$PValue,method = "BH")
TestedPeaks<-TestedPeaks[,c("region","logCPM","logFC","PValue","FDR")]

volcanoData <- cbind(TestedPeaks$logFC, -log10(TestedPeaks$FDR))
colnames(volcanoData) <- c("logFC", "-LogPval")
DEGs <- TestedPeaks$FDR < 0.05 & TestedPeaks$logFC>0
point.col <- ifelse(DEGs, "red", "black")

pdf(paste0("220703_gen_countTable_output/","peak_volcanoData.pdf"))
plot(volcanoData, pch = 16, col = point.col, cex = 0.5, xlim=c(-5,5))
dev.off()


TestedPeaks_up<-TestedPeaks[TestedPeaks$logFC>0 & TestedPeaks$FDR<0.05,]
TestedPeaks_up_regions<-str_split_fixed(TestedPeaks_up$region,"_",3)
TestedPeaks_up$chr<-TestedPeaks_up_regions[,1]
TestedPeaks_up$start<-as.integer(TestedPeaks_up_regions[,2])
TestedPeaks_up$end<-as.integer(TestedPeaks_up_regions[,3])
TestedPeaks_up_bed<-TestedPeaks_up[,c("chr","start","end")]
write.table(TestedPeaks_up_bed, paste0(output_dir,"/homer_input/TestedPeaks_up_0.05.bed"),quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")

TestedPeaks_up_01<-TestedPeaks[TestedPeaks$logFC>0 & TestedPeaks$FDR<0.01,]
TestedPeaks_up_01_regions<-str_split_fixed(TestedPeaks_up_01$region,"_",3)
TestedPeaks_up_01$chr<-TestedPeaks_up_01_regions[,1]
TestedPeaks_up_01$start<-as.integer(TestedPeaks_up_01_regions[,2])
TestedPeaks_up_01$end<-as.integer(TestedPeaks_up_01_regions[,3])
TestedPeaks_up_01_bed<-TestedPeaks_up_01[,c("chr","start","end")]
write.table(TestedPeaks_up_01_bed, paste0(output_dir,"/homer_input/TestedPeaks_up_0.01.bed"),quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")

TestedPeaks_up_001<-TestedPeaks[TestedPeaks$logFC>0 & TestedPeaks$FDR<0.001,]
TestedPeaks_up_001_regions<-str_split_fixed(TestedPeaks_up_001$region,"_",3)
TestedPeaks_up_001$chr<-TestedPeaks_up_001_regions[,1]
TestedPeaks_up_001$start<-as.integer(TestedPeaks_up_001_regions[,2])
TestedPeaks_up_001$end<-as.integer(TestedPeaks_up_001_regions[,3])
TestedPeaks_up_001_bed<-TestedPeaks_up_001[,c("chr","start","end")]
write.table(TestedPeaks_up_001_bed, paste0(output_dir,"/homer_input/TestedPeaks_up_0.001.bed"),quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")


TestedPeaks_down<-TestedPeaks[TestedPeaks$logFC<0 & TestedPeaks$FDR<0.05,]
TestedPeaks_down_regions<-str_split_fixed(TestedPeaks_down$region,"_",3)
TestedPeaks_down$chr<-TestedPeaks_down_regions[,1]
TestedPeaks_down$start<-as.integer(TestedPeaks_down_regions[,2])
TestedPeaks_down$end<-as.integer(TestedPeaks_down_regions[,3])
TestedPeaks_down_bed<-TestedPeaks_down[,c("chr","start","end")]
write.table(TestedPeaks_down_bed, paste0(output_dir,"/homer_input/TestedPeaks_down_0.05.bed"),quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")

TestedPeaks_down<-TestedPeaks[TestedPeaks$logFC<0 & TestedPeaks$FDR<0.01,]
TestedPeaks_down_regions<-str_split_fixed(TestedPeaks_down$region,"_",3)
TestedPeaks_down$chr<-TestedPeaks_down_regions[,1]
TestedPeaks_down$start<-as.integer(TestedPeaks_down_regions[,2])
TestedPeaks_down$end<-as.integer(TestedPeaks_down_regions[,3])
TestedPeaks_down_bed<-TestedPeaks_down[,c("chr","start","end")]
write.table(TestedPeaks_down_bed, paste0(output_dir,"/homer_input/TestedPeaks_down_0.01.bed"),quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")

TestedPeaks_down<-TestedPeaks[TestedPeaks$logFC<0 & TestedPeaks$FDR<0.001,]
TestedPeaks_down_regions<-str_split_fixed(TestedPeaks_down$region,"_",3)
TestedPeaks_down$chr<-TestedPeaks_down_regions[,1]
TestedPeaks_down$start<-as.integer(TestedPeaks_down_regions[,2])
TestedPeaks_down$end<-as.integer(TestedPeaks_down_regions[,3])
TestedPeaks_down_bed<-TestedPeaks_down[,c("chr","start","end")]
write.table(TestedPeaks_down_bed, paste0(output_dir,"/homer_input/TestedPeaks_down_0.001.bed"),quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")

TestedPeaks_regions<-str_split_fixed(TestedPeaks$region,"_",3)
TestedPeaks$chr<-TestedPeaks_regions[,1]
TestedPeaks$start<-as.integer(TestedPeaks_regions[,2])
TestedPeaks$end<-as.integer(TestedPeaks_regions[,3])
TestedPeaks_bed<-TestedPeaks[,c("chr","start","end")]
write.table(TestedPeaks_bed, paste0(output_dir,"/homer_input/TestedPeaks_background.bed"),quote = FALSE, row.names = FALSE,col.names = FALSE, sep="\t")
