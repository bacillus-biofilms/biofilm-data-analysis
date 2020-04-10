# set folder prefix and project folder and set them as working directory

folderPrefix <- "path_to_your_folder_prefix"
projectFolder <- "your_project_folder_name"

setwd(file.path(folderPrefix, projectFolder))

# read necessary libraries. Please make sure they are all previously installed.

library(biomartr)
library(coRdon)
library(tidyverse)
library(plyr)
library(stringr)
library(myTAI)
library(orthologr)



# make and save input for TAI

# read file with locus_tag and raw counts divided by gene length and resolved replicates

preprocessed <- read.table(file="name_of_file_with_locus_tag_and_resolved_replicates", header=T, sep = '\t')


# read table containing protein_id, locus_tag and gene_length obtained from feature table file from GenBank

pid_lt_gl <- read.table(file="name_of_file_with_protein_id_locus_tag_gene_length", header=T, sep = '\t')


# read table containing protein_id and ps obtained from phylostratigraphy

map <- read.table(file="name_of_file_with_protein_id_and_ps", header=T, sep = '\t')


# join map with locus_tag and raw counts divided by gene length and replicates resolved

map_pid_lt_gt <- join (map, pid_lt_gl, by="protein_id", type="inner")
input_for_TAI <- join (map_pid_lt_gt, preprocessed, by="locus_tag", type="inner")


# remove unnecessary columns and rename columns

input_for_TAI <- input_for_TAI[,c(2,3,5:15)]
colnames(input_for_TAI) <- c("ps", "locus_tag", "LC", "6H", "12H", "1D", "2D", "3D", "5D", "7D", "14D", "1M", "2M")


# save input for TAI

write.table(input_for_TAI, "name_your_input_for_TAI", sep="\t", quote=F, row.names=F)


# plot TAI and save

TAI_plot <- PlotSignature(input_for_TAI, measure = "TAI",
                          TestStatistic = "FlatLineTest", modules = NULL, permutations = 1000,
                          lillie.test = FALSE, p.value = TRUE, shaded.area = FALSE,
                          custom.perm.matrix = NULL, xlab = "Biofilm stages",
                          ylab = "TAI", main = "", lwd = 1, alpha = 0.1,
                          y.ticks = 10)

ggsave(filename = "name_your_TAI_plot", width = 10, height = 10, dpi = 300)










# dNdS calculation

# calculate dN and dS between Bacillus subtilis subsp. subtilis str. NCIB3610 and Bacillus licheniformis DSM 13 = ATCC 14580 and save

dNdS_calc <- dNdS(query_file      = "name_of_Bs3610_protein_file",
                 subject_file    = "name_of_BlDSM13_protein_file",
                 seq_type        = "protein",
                 delete_corrupt_cds = FALSE, 
                 ortho_detection = "RBH",
                 aa_aln_type     = "pairwise",
                 aa_aln_tool     = "NW",
                 codon_aln_tool  = "pal2nal",
                 dnds_est.method = "YN",
                 store_locally = F,
                 comp_cores      = 1 )

write.table(x         = dNdS_calc, 
            file      = "name_your_dNdS_Bs_Bl_output_file", 
            sep       = "\t", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote     = FALSE)






# make and save input for TdNI

dNdS_calc <- read.table(file="path_to_file_containing_dN_and_dS_values_for_Bs_and_Bl", header=T, sep = '\t')
protein_id <- str_extract(dNdS_calc$query_id, 'AQZ[^_]*')
dN <- cbind (protein_id, dNdS_calc)
dN <- dN[,c("dN", "protein_id")]

dN_locus_tag <- join (dN, pid_lt_gl, by="protein_id", type="inner")

input_for_TdNI <- join (dN_locus_tag, input_for_TAI, by="locus_tag", type="inner")
input_for_TdNI <- input_for_TdNI[,c("dN", "locus_tag", "LC", "6H", "12H", "1D", "2D", "3D", "5D", "7D", "14D", "1M", "2M")]

write.table(input_for_TdNI, "name_your_input_file_for_TdNI", sep="\t", quote=F, row.names=F)


# plot TdNI

TdNI_plot <- PlotSignature(input_for_TdNI, measure = "TDI",
                          TestStatistic = "FlatLineTest", modules = NULL, permutations = 1000,
                          lillie.test = FALSE, p.value = TRUE, shaded.area = FALSE,
                          custom.perm.matrix = NULL, xlab = "Biofilm stages",
                          ylab = "TdNI", main = "", lwd = 1, alpha = 0.1,
                          y.ticks = 10)

ggsave(filename = "name_your_TdNI_plot", width = 10, height = 10, dpi = 300)






# make and save input for TdSI

dS <- cbind (protein_id, dNdS_calc)
dS <- dS[,c("dS", "protein_id")]

dS_locus_tag <- join (dS, pid_lt_gl, by="protein_id", type="inner")

input_for_TdSI <- join (dS_locus_tag, input_for_TAI, by="locus_tag", type="inner")
input_for_TdSI <- input_for_TdSI[,c("dS", "locus_tag", "LC", "6H", "12H", "1D", "2D", "3D", "5D", "7D", "14D", "1M", "2M")]
input_for_TdSI <- input_for_TdSI[complete.cases(input_for_TdSI), ]

write.table(input_for_TdSI, "name_your_input_file_for_TdSI", sep="\t", quote=F, row.names=F)


# plot TdSI and save

TdSI_plot <- PlotSignature(input_for_TdSI, measure = "TDI",
                           TestStatistic = "FlatLineTest", modules = NULL, permutations = 1000,
                           lillie.test = FALSE, p.value = TRUE, shaded.area = FALSE,
                           custom.perm.matrix = NULL, xlab = "Biofilm stages",
                           ylab = "TdSI", main = "", lwd = 1, alpha = 0.1,
                           y.ticks = 10)

ggsave(filename = "name_your_TdSI_plot", width = 10, height = 10, dpi = 300)









# TCBI calculation

# calculate codon usage bias for Bs3610

bacillus_genes <- read_cds("GCA_002055965.1_ASM205596v1_cds_from_genomic.fna")
bacillus_ribosomal_genes <- read_cds("ribosomal_Bs_3610.fasta")

codon_table = Biostrings::DNAStringSet(bacillus_genes) %>% codonTable()
codon_table_ribosomal = Biostrings::DNAStringSet(bacillus_ribosomal_genes) %>% codonTable()

ref_list = list()
ref_list[["ribosomal_61"]] = codon_table_ribosomal

milc <- MILC(codon_table, subsets = ref_list)

milc <- milc %>% as.data.frame() %>%
  add_column(.before = 1, gene = names(bacillus_genes))
locus_tag <- str_extract(milc$gene, 'B4U62_[^]]*')
milc <- cbind(locus_tag, milc)
milc <- milc[,c("ribosomal_61", "locus_tag")]


# make input for TCBI and save

input_for_TCBI <- join (milc, input_for_TAI, by="locus_tag", type="inner")
input_for_TCBI <- input_for_TCBI[,c("ribosomal_61", "locus_tag", "LC", "6H", "12H", "1D", "2D", "3D", "5D", "7D", "14D", "1M", "2M")]

write.table(input_for_TCBI, "name_your_input_file_for_TCBI", sep="\t", quote=F, row.names=F)


# plot TCBI and save

TCBI_plot <- PlotSignature(input_for_TCBI, measure = "TAI",
                           TestStatistic = "FlatLineTest", modules = NULL, permutations = 1000,
                           lillie.test = FALSE, p.value = TRUE, shaded.area = FALSE,
                           custom.perm.matrix = NULL, xlab = "Biofilm stages",
                           ylab = "TCBI", main = "", lwd = 1, alpha = 0.1,
                           y.ticks = 10)

ggsave(filename = "name_your_TCBI_plot", width = 10, height = 10, dpi = 300)