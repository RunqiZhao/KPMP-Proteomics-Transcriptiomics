# Tissue EnrichC
# BiocManager::install("TissueEnrich")
library(TissueEnrich)
library(tidyr)
library(dplyr)
library(reshape2)

Proteins_KPMP <- readxl::read_excel("~/02 Enrichment/Output/proteinatlas_significant_protein_enrichment_positive_only.05.22.2025.xlsx")
Proteins_KPMP$Gene[which(!is.na(Proteins_KPMP$Ensembl))]
geneList <- unique(Proteins_KPMP$Gene)
geneList <- na.omit(geneList)
gs <- GeneSet(geneIds = geneList,
              organism = 'Homo Sapiens',
              geneIdType = SymbolIdentifier())

output <- teEnrichment(gs, tissueSpecificGeneType = 1)

enrichmentOutput <- setNames(data.frame(assay(output[[1]]),
                                        row.names = rowData(output[[1]])[,1]),
                             colData(output[[1]])[,1])
enrichmentOutput$Tissue <- row.names(enrichmentOutput)

#plot the P-Values of enrichment
ggplot(enrichmentOutput,aes(x = reorder(Tissue, -Log10PValue),
                            y = Log10PValue,
                            label = Tissue.Specific.Genes,fill = Tissue))+
  geom_bar(stat = 'identity')+
  labs(x = '', y = '-LOG10(P-Value)')+
  theme_bw()+
  theme(legend.position='none')+
  theme(plot.title = element_text(hjust = 0.5,size = 20),
        axis.title = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# search the expression of each protein/Gene
exp_all <- NA
for (i in 1:dim(enrichmentOutput)[1]){
  seExp <- output[[2]][[i]]
  exp <- setNames(data.frame(assay(seExp), 
                             row.names = rowData(seExp)[,1]), colData(seExp)[,1])
  exp$Gene <- row.names(exp)
  
  exp <- melt(exp, id = "Gene")
  colnames(exp) <- c("Gene", "key", "value")
  
  exp_all <- rbind(exp_all, exp) %>% as.data.frame()
}
exp_all <- unique(exp_all) %>% na.omit()

png("Output/TISSUE-enrichment of HPA.tiff", res = 1200, width = 12000, height = 10800)
ggplot(exp_all, aes(key, Gene)) + 
  geom_tile(aes(fill = value),
            colour = "white") + 
  scale_fill_gradient(low = "white",
                      high = "red")+
  labs(x='', y = '')+
  theme_bw()+
  guides(fill = guide_legend(title = "Log2(TPM)"))+
  theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),panel.grid.major= element_blank(),panel.grid.minor = element_blank())
dev.off()
