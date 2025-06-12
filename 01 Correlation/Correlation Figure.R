library(readxl)
library(dplyr)
library(tidyr)

library(ggplot2)
library(corrplot)
library(RColorBrewer)

Data <- read_excel("Combined_Bulk_scRNA_ANMLNormalized.xlsx", sheet = "Spearman")
Names <- data.frame(ColName = colnames(Data))
Names <- Names %>% filter(!ColName %in% c("Protein", "Gene")) %>% 
  separate(ColName, into = c("Bulk_Names", "Other"), sep = "_", remove = F)
Bulk_Names <- unique(Names$Bulk_Names)

Protein_Hits <- data.frame(Protein = NA, Gene = NA)

for(Name in Bulk_Names){
  Data_Name <- select(Data, c("Protein", "Gene", Names$ColName[which(Names$Bulk_Names == Name)]))
  p_index <- which(grepl("value$", colnames(Data_Name)))
  # hits defined as Spearman p<0.001 with positive correlation
  Sig_Inx <- which(Data_Name[, p_index] < 0.001)
  Data_Name[Data_Name < 0] <- NA
  
  Name_Hits <-  Data_Name[Sig_Inx, ] %>% na.omit()
  Protein_Hits <- merge(Protein_Hits, Name_Hits, by.x = c("Protein", "Gene"), by.y = c("Protein", "Gene"), all.x = T, all.y = T)
}

Protein_Hits <- filter(Protein_Hits, !is.na(Protein)) %>% unique()

Protein_Hits %>% writexl::write_xlsx("Output/Cor_Sigificant_Positive.xlsx")

# For same Protein Name, keep the most significant one
# Selection <- Protein_Hits %>% group_by(Gene) %>% mutate(count = n()) %>% filter(count > 1)
Protein_Hits <- filter(Protein_Hits, !Protein %in% c("SIGLEC14_8248_222", "SIGLEC14_5125_6",
                                                     "TIMP1_23173_3", "TIMP1_2211_9"))
row.names(Protein_Hits) <- Protein_Hits$Gene

CorFigure_DT <- select(Protein_Hits, -c("Protein", "Gene"))

cor_col <- colnames(CorFigure_DT)[which(grepl("Spearman$", colnames(CorFigure_DT)))]
p_value_col <- colnames(CorFigure_DT)[which(grepl("value$", colnames(CorFigure_DT)))]

CorFigure_DT_cor <- select(CorFigure_DT, c(all_of(cor_col)))
CorFigure_DT_cor <- t(CorFigure_DT_cor)
rname <- rownames(CorFigure_DT_cor) %>% as.data.frame()
rname <- separate(rname,  . , c("Part", "Spearman"), sep = "_")
rownames(CorFigure_DT_cor) <- rname$Part

CorFigure_DT_p <- select(CorFigure_DT, c(all_of(p_value_col)))
CorFigure_DT_p <- t(CorFigure_DT_p)
rname <- rownames(CorFigure_DT_p) %>% as.data.frame()
rname <- separate(rname,  . , c("Part", "Spearman", "p-value"), sep = "_")
rownames(CorFigure_DT_p) <- rname$Part

CorFigure_DT <-  list(r = CorFigure_DT_cor, p = CorFigure_DT_p)

png(filename = "Output/Cor_heatmap_sigificant positive.05.22.2025.png", res=1200, width = 30000, height = 6800)
corrplot(CorFigure_DT$r, p.mat = CorFigure_DT$p, 
         is.corr = FALSE,
         type = "full", order = "original",
         col = colorRampPalette(c("#fff5f0", "#fed976", "#feb24c", "#fc4e2a", "#800026"))(5),
         na.label = " ", 
         insig = "label_sig", sig.level = c(0.00001, 0.0001, 0.001), pch.cex = 1.5,pch.col = "grey60",
         tl.srt = 45, tl.col = "black", cl.pos = "b",
         col.lim = c(0.5, 1)
         )
dev.off()

# Present Other Coefficients
Protein_Hits_Other <- filter(Data, Protein %in% Protein_Hits$Protein) %>% as.data.frame()
Protein_Hits_Other %>% writexl::write_xlsx("Output/Cor_Show_all.05.28.2025.xlsx")

row.names(Protein_Hits_Other) <- Protein_Hits_Other$Gene
Protein_Hits_Other <- select(Protein_Hits_Other, -c("Protein", "Gene"))

Protein_Hits_Other_cor <- select(Protein_Hits_Other, c(all_of(cor_col)))
Protein_Hits_Other_cor <- t(Protein_Hits_Other_cor)
rname <- rownames(Protein_Hits_Other_cor) %>% as.data.frame()
rname <- tidyr::separate(rname,  . , c("Part", "Spearman"), sep = "_")
rownames(Protein_Hits_Other_cor) <- rname$Part
# Show Positive only
Protein_Hits_Other_cor[Protein_Hits_Other_cor < 0] <- NA
which(is.na(Protein_Hits_Other_cor), arr.ind = T)

Protein_Hits_Other_p <- select(Protein_Hits_Other, c(all_of(p_value_col)))
Protein_Hits_Other_p <- t(Protein_Hits_Other_p)
row.names(Protein_Hits_Other_p) <- row.names(Protein_Hits_Other_cor)
Protein_Hits_Other_p[which(is.na(Protein_Hits_Other_cor), arr.ind = T)] <- NA

Protein_Hits_Other <-  list(r = Protein_Hits_Other_cor, p = Protein_Hits_Other_p)

png(filename = "Output/Cor_heatmap_show_all.05.22.2025.Positive only.png", res=1200, width = 30000, height = 6800)
corrplot(Protein_Hits_Other$r, 
         p.mat = Protein_Hits_Other$p, 
         type = "full", order = "original",
         col = colorRampPalette(c("#08306b",  "#2171b5", "#6baed6", "#c6dbef", "white", 
                                  "#fff5f0", "#fed976", "#feb24c", "#fc4e2a", "#800026"))(10),
         na.label = " ", 
         insig = "label_sig", sig.level = c(0.00001, 0.0001, 0.001), pch.cex = 1.5, pch.col = "grey90",
         col.lim = c(0,1),
         tl.srt = 45, tl.col = "black", cl.pos = "b")
dev.off()

# Summary table ----
library(gtsummary)
Tbl_DT <- Protein_Hits %>% 
  select(-c("Protein", "Gene")) %>% 
  apply(2, as.numeric) %>% 
  as.data.frame() 

count <- Protein_Hits %>% 
  select(-c("Protein", "Gene")) %>% 
  sapply(is.na) %>% 
  as.data.frame() %>% 
  sapply(sum, na.rm = T) * (-1) + 94

Min <- Protein_Hits %>% 
  select(-c("Protein", "Gene")) %>% 
  sapply(min, na.rm = T)

Max <- Protein_Hits %>% 
  select(-c("Protein", "Gene")) %>% 
  sapply(max, na.rm = T)

SummaryTable <- cbind(count, Min, Max) %>% as.data.frame()
SummaryTable$Min <- sprintf("%.3f", SummaryTable$Min)
SummaryTable$Max <- sprintf("%.3f", SummaryTable$Max)
SummaryTable$Colnames <- rownames(SummaryTable)

SummaryTable <- separate(SummaryTable, Colnames, c("Part", "Cal", "Cal2"), sep = "_")
SummaryTable <- filter(SummaryTable, is.na(Cal2)) %>% select(c("Part", "count", "Min", "Max"))

SummaryTable$count <- paste0("N = ", SummaryTable$count)
SummaryTable$MM <- paste0("(", SummaryTable$Min, ", ", SummaryTable$Max, ")")

writexl::write_xlsx(SummaryTable, "Output/SummaryTable.xlsx")
