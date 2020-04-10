#' ---
#' title: "Bacillus subtilis 3610 Transcriptome analysis - biofilm"
#' date: "`r format(Sys.time(), '%d %B %Y')`"
#' output: html_notebook
#' ---
#' 
#' 
#+ setup-libs-folders
options(java.parameters = "-Xmx10000m")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(GenomicAlignments)
  library(Rsamtools)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(PoiClaClu)
  library(ggplot2)
  library(ggrepel)
  library(gage)
  library(pathview)
  library(biomaRt)
  library(xlsx)
  library(tidyverse)
})


folderPrefix <- "path_to_your_folder_prefix"


projectFolder <- "your_project_folder_name"
inFolder <- paste(folderPrefix, projectFolder, "bam", sep = "/")
outFolder <- paste(folderPrefix, projectFolder, "count", sep = "/")

setwd(file.path(folderPrefix, projectFolder))

# Java garbage collection
jgc <- function() {
  .jcall("java/lang/System", method = "gc")
}

#' 
#'  
#' 
#' 
#+ load-genome 
chrfile <- file.path(outFolder, "file_with_Bs3610_chrsizes")
chrominfo <- read.delim(chrfile, header = FALSE)
names(chrominfo) <- c("chrom", "length")
chrominfo$is_circular <- FALSE
chrominfo[1, "is_circular"] <- TRUE



gtffile <- file.path(outFolder, "Bs3610.gff")
gff_granges <- rtracklayer::import.gff3(gtffile, feature.type = c("gene", "pseudogene"))
names(gff_granges) <- gff_granges$locus_tag


gtfdescfile <- file.path(outFolder, "Bs3610.gff")
gff_desc <- rtracklayer::import.gff3(gtfdescfile, feature.type = c("gene", "pseudogene"))



#+ count-fragments-unix, eval=FALSE
p1 <- ScanBamParam()


# *._sorted.bam$

infiles <- file.path(inFolder, list.files(inFolder, patt="*._sorted.bam$"))
bamfiles <- BamFileList(infiles, yieldSize=2000000, asMates = TRUE)

sampleTable <- data.frame(filePath = infiles,
                          filename = sub(paste(inFolder, "/", sep = ""), "", infiles), stringsAsFactors = FALSE)
sampleTable$sample <- sub("_.*", "", sampleTable$filename)


save(sampleTable, file = file.path(outFolder, "Bs3610_biofilm_Transcriptome_sampleTable.Robj"))

se <- summarizeOverlaps(features=gff_granges, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        param=p1)

save(se, file = file.path(outFolder, "Bs3610_biofilm_Transcriptome_se_gff.Robj"))

#' 
#' 
#' 
#' 
#+ load-counts 
sample_desc <- read_tsv("Sample_Name_biofilm.txt") %>%
  separate(RealName, into = c("Time", "Replicate"))

time_points <- c("LC", "6H", "12H", "1D", "2D", "3D", "5D", "7D", "14D", "1M", "2M")


sampleTableFile <- file.path(outFolder, "Bs3610_biofilm_Transcriptome_sampleTable.Robj")
load(sampleTableFile)
sample_table_final <- sampleTable %>%
  mutate(sample = paste("B", sample, sep = "")) %>% 
  left_join(sample_desc, by=c("sample" = "RNAseq")) %>%
  mutate(Time = factor(Time, levels = time_points))


countFile <- file.path(outFolder, "Bs3610_biofilm_Transcriptome_se_gff.Robj")
load(countFile)
colData(se) <- DataFrame(sample_table_final)
colnames(se) <- sample_table_final$sample
rownames(se) <- rownames(se) %>% make.unique()

dds <- DESeqDataSet(se, design = ~ Time)

dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# Remove ribosomal sequences
dds <- dds[! grepl("RNA", rownames(dds)),]
nrow(dds)


#+ rlog-normalize
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)

#+ export-counts, eval=FALSE
# export raw counts
rawCounts <- assay(se)
write.table(rawCounts, file="Bs3610_biofilm_raw_counts.txt", sep="\t",
            row.names = TRUE,
            quote = FALSE)




# Plot PCA    
pca_data <- plotPCA(rld, intgroup = c("Time"), returnData = TRUE)
perc_var <- attr(pca_data, "percentVar")
pca_plot <- pca_data %>%
  left_join(sample_desc, by=c("name" = "RNAseq", "Time"))

ggplot(pca_plot, aes(x=PC1, y=PC2, color=Time)) +
  geom_point(size=6) +
  scale_color_manual(values = c('#9C3234', '#076D94', '#ED3237', '#61BF1A', '#AED8E8', '#4F8C0D', '#52C6F1', '#00AFEF', '#ECA7AB', '#0A96CC', '#606062'), guide="none") +
  scale_fill_manual(values = c('#9C3234', '#076D94', '#ED3237', '#61BF1A', '#AED8E8', '#4F8C0D', '#52C6F1', '#00AFEF', '#ECA7AB', '#0A96CC', '#606062'), guide="none") +
  labs(title = "",
       subtitle = "",
       x = paste0("PC1: ",round(perc_var[1] * 100),"% variance"),
       y = paste0("PC2: ",round(perc_var[2] * 100),"% variance"),
       fill = "Group",
       colour = "Time",
       shape = "Time") +
  theme(panel.grid = element_line(color = 'gray')) + 
  theme(panel.background = element_rect(fill = NA))


ggsave(filename = "Bs3610_biofilm_PCA.pdf", width = 10, height = 10, dpi = 300)



#Likelihood ratio test - LRT full dataset

ddsTS <- dds
ddsTS <- DESeq(ddsTS, test = "LRT", full = ~Time, reduced = ~1)
ddsTS <- results(ddsTS)

head(ddsTS[order(ddsTS$padj),], 10)

write.table(ddsTS, "biofilm_LRT_full_dataset.txt", sep="\t", quote=F, row.names=T)


#Likelihood ratio test - LRT reduced dataset (without LC, 1M, 2M)

se_reduced <- se[,c("B4", "B5", "B6", "B22", "B23", "B24", "B7", "B8", "B9", "B10", "B11", "B12", "B13", "B14", "B15", "B16", "B17", "B18", "B1", "B2", "B3", "B19", "B20", "B21")]

id_table <- mcols(gff_granges) %>%
  as.data.frame()

dds_reduced <- DESeqDataSet(se_reduced, design = ~ Time)
dds_reduced <- dds_reduced[ rowSums(counts(dds_reduced)) > 1, ]
dds_reduced <- DESeq(dds_reduced)

ddsTS_reduced <- dds_reduced
ddsTS_reduced <- DESeq(ddsTS_reduced, test = "LRT", full = ~Time, reduced = ~1)
ddsTS_reduced <- results(ddsTS_reduced)

head(ddsTS_reduced[order(ddsTS_reduced$padj),], 10)

write.table(ddsTS_reduced, "biofilm_LRT_reduced.txt", sep="\t", quote=F, row.names=T)


# Differential expression Excel generation

id_table <- mcols(gff_granges) %>%
  as.data.frame()

dds <- DESeqDataSet(se, design = ~ Time)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)


# Pairwise differential expression


dif_result <- function(V1, V2, ...) {
  gc()
  
  from <- V2
  to <- V1
  
  resSig <- results(dds, contrast=c("Time", from, to))
  resSig$locus_tag <- rownames(resSig)
  resSig <- resSig %>%
    as.data.frame() %>%
    left_join(id_table[, c("locus_tag", "gene_biotype", "gene", "Name")], by = c("locus_tag")) %>%
    column_to_rownames("locus_tag")
  
  sheet <- createSheet(wb, sheetName = paste0(from, " vs. ", to))
  addDataFrame(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ], 
               sheet,
               startRow=1,
               startColumn=1,
               colnamesStyle=cs3,
               rownamesStyle=cs1)
  autoSizeColumn(sheet, 1:(length(resSig)+1))
  
}

wb <- createWorkbook()
cs1 <- CellStyle(wb) + Font(wb, isItalic=TRUE)           # rowcolumns
cs3 <- CellStyle(wb) + Font(wb, isBold=TRUE) + Border()  # header

combn(time_points,2) %>% t() %>% as_tibble() %>% pwalk(dif_result)



saveWorkbook(wb, file = "Bs3610_biofilm_all_vs_all_diffex.xlsx")
