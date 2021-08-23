###ATAC-seq
options(stringsAsFactors = F)


library(DiffBind)
setwd("/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB")

#Reading in the peaksets
samples <- read.table("diffbind/spt5_ATAC_sampleinfo_DB.txt",header = T ,sep="\t")
samples

tamoxifen <- dba(sampleSheet=samples)
tamoxifen



#correlation heatmap can be generated which gives an initial clustering of 
#the samples using the cross-correlations of each row of the binding matrix
#!!!using occupancy (peak caller score) data
pdf(paste0("/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/",
           "Heatmap_sample_clustering_usingpeakscore.pdf"),width = 5,height = 5)
plot(tamoxifen)
dev.off()


#Counting reads
tamoxifen <- dba.count(tamoxifen)
tamoxifen


#For each sample, multiplying the value in the Reads column by the corresponding FRiP value
#will yield the number of reads that overlap a consensus peak. 
info <- dba.show(tamoxifen)
libsizes <- cbind(LibReads=info$Reads, FRiP=info$FRiP,
                    PeakReads=round(info$Reads * info$FRiP))

rownames(libsizes) <- info$ID
libsizes
write.table(libsizes,"/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/libsizes.txt",
            row.names = T,col.names = NA,sep ="\t",quote=F)


pdf(paste0("/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/",
           "Heatmap_sample_clustering_using_readcount.pdf"),width = 5,height = 5)
plot(tamoxifen)
dev.off()


tamoxifen <- dba.normalize(tamoxifen)
norm <- dba.normalize(tamoxifen, bRetrieve=TRUE)
norm

normlibs <- cbind(FullLibSize=norm$lib.sizes, NormFacs=norm$norm.factors,
                  NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs) <- info$ID
normlibs

#Establishing a model design and contrast
tamoxifen <- dba.contrast(tamoxifen,
                          reorderMeta=list(Condition="DMSO"),minMembers = 2)
#tamoxifen <- dba.contrast(tamoxifen,design="~Condition")
tamoxifen <- dba.analyze(tamoxifen)
dba.show(tamoxifen, bContrasts=TRUE)

#Retrieving the differentially bound sites
tamoxifen.DB <- dba.report(tamoxifen)
result <- as.data.frame(tamoxifen.DB)
write.table(result,"/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/spt5_db_grange.txt",
            row.names = F,quote=F,sep="\t")



####plot
#PCA plot
pdf(paste0("/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/",
           "PCA.pdf"),width = 7,height = 6)
dba.plotPCA(tamoxifen,DBA_TISSUE,label=DBA_CONDITION)
dev.off()

pdf(paste0("/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/",
           "PCA_DBpeaks.pdf"),width = 7,height = 6)
dba.plotPCA(tamoxifen, contrast=1, label=DBA_REPLICATE)
dev.off()

#MA plot
pdf(paste0("/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/",
           "MA_allpeaks.pdf"),width = 5,height = 4.5)
dba.plotMA(tamoxifen)
dev.off()


#Volcano plot
pdf(paste0("/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/",
           "Volcano_allpeaks.pdf"),width = 5,height = 4)
dba.plotVolcano(tamoxifen)
dev.off()

#Boxplots
pdf(paste0("/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/",
           "Boxplots.pdf"),width = 7,height = 6)
dba.plotBox(tamoxifen)
dev.off()

pvals <- dba.plotBox(tamoxifen)
write.table(pvals,"/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/spt5_diffbind_boxplot_pvalue.txt",
            row.names = T,col.names = NA,quote=F,sep="\t")


#Heatmap
hmap <- colorRampPalette(c("#005F8DFF","#06AFCBFF","#C6E4F0FF",
                           "#F9FAFC", "#FCE0CD", "#FFAF8A", "#FD705F", "#E12234","#840729"), alpha = TRUE)(50)
pdf(paste0("/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/",
           "Heatmap_diffbind_peak_readcount.pdf"),width = 10,height = 10)
dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE,
                               scale="row", colScheme = hmap)
dev.off()


hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
pdf(paste0("/share/home/Blueberry/Projects/Chenlab/Hushibin/SPT5/ATAC-seq/DLD1/06_analysis/06.1_SPT5_DB/diffbind/",
           "Heatmap_diffbind_peak_readcount.pdf"),width = 10,height = 10)
dba.plotHeatmap(tamoxifen, contrast=1, correlations=FALSE,
                               scale="row", colScheme = hmap)
dev.off()



