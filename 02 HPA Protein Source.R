library(dplyr)
library(reshape2)

Proteins_KPMP <- readxl::read_excel("~/01 Correlation/Output/Cor_Sigificant_Positive.xlsx")
cor_col <- colnames(Proteins_KPMP)[which(grepl("Spearman$", colnames(Proteins_KPMP)))]
p_value_col <- colnames(Proteins_KPMP)[which(grepl("value$", colnames(Proteins_KPMP)))]

Proteins_KPMP <- select(Proteins_KPMP, c("Protein", "Gene", all_of(cor_col)))
# Keep Unique Gene Name
Proteins_KPMP <- filter(Proteins_KPMP, !Protein %in% c("SIGLEC14_8248_222", "SIGLEC14_5125_6",
                                                       "C1QTNF4_25964_12", 
                                                       "TIMP1_23173_3", "TIMP1_2211_9", 
                                                       "PNP_10039_32"))
Proteins_KPMP <- melt(Proteins_KPMP, id = c("Protein", "Gene"))
Proteins_KPMP <- na.omit(Proteins_KPMP)

Proteins_KPMP$variable <- gsub("_Spearman", "", Proteins_KPMP$variable)

library(data.table)
# https://www.proteinatlas.org/about/download
protein <- read.table("Data/proteinatlas.tsv", sep = '\t', header = TRUE)
protein <- filter(protein, Gene %in% Proteins_KPMP$Gene)
test <- colnames(protein) %>% as.data.frame()
protein_enrichment <- select(protein, c("Gene", "Ensembl", "Uniprot", 
                                        "RNA.tissue.cell.type.enrichment", "RNA.tissue.specific.nTPM"))

protein_enrichment <- merge(Proteins_KPMP, protein_enrichment, by.x = "Gene", by.y = "Gene", all.x = T)
writexl::write_xlsx(protein_enrichment, "Output/proteinatlas_significant_protein_enrichment_positive_only.05.22.2025.xlsx")
