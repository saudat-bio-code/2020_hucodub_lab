##### Saudat
setwd("/Users/saudaalishaeva/Desktop/Bekpen_Lab_Project/data/")
#setwd("/Volumes/TOSHIBA/archieved_Work/HuCoDup/PhyRepID-master/pipeline/logs/")
library(clusterProfiler)
library("biomaRt")
library(readr)
library(org.Hs.eg.db)
library(data.table)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene) #always check annotation
OrgDb <- org.Hs.eg.db


#library(DOSE)
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevelsStyle(background_info) <- "UCSC"

background_info = read_tsv("phyrepid_results_human_full_lineage.tsv")

background_info = background_info %>%
  separate(identifier, c("tree", "geneId", "other"), "_")

#setwd("/Volumes/TOSHIBA/archieved_Work/HuCoDup/PhyRepID-master/pipeline/logs/")
#library(jsonlite)
#orthologs <- fromJSON("orthologs_filtered.json", flatten=TRUE)
df_orth = as.data.frame(names(orthologs))

df_orth = df_orth %>%
  separate(ens, c("link", "id"), "/ensembl/")
bgr_orthologs = as.data.frame(df_orth$id)
write_csv(bgr_orthologs,"bgr_orthologs.csv")

### always check gene names

ensembl=useMart("ensembl")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
biomart_annotated_2<-getBM(attributes=c('ensembl_gene_id', 'entrezgene_id'), 
                         filters = 'ensembl_gene_id', 
                         values = background_info$geneId, uniqueRows = F,
                         mart = ensembl)

background_info$geneId = 
#somehow you need to add results with p values to the df_ann_results

not_redundant_files<-background_info[(background_info$geneId %in% biomart_annotated_2$ensembl_gene_id), ]
#not_redundant_files$entrez<-biomart_annotated_2[( biomart_annotated_2$hgnc_symbol %in% background_info$gene_symbol ),]
not_redundant_files = not_redundant_files[!duplicated(not_redundant_files$geneId), ]
not_redundant_files = as.data.frame(not_redundant_files)
rownames(not_redundant_files) = not_redundant_files$geneId
biomart_annotated_3 = biomart_annotated_2[!duplicated(biomart_annotated_2$ensembl_gene_id), ]
rownames(biomart_annotated_3) =  biomart_annotated_3$ensembl_gene_id
merged_annotat_results <- merge(not_redundant_files,biomart_annotated_3, by= 0, all=TRUE) #by=0 means it will merge by raw names
merged_annotat_results_rmv_na = merged_annotat_results[complete.cases(merged_annotat_results[ ,16]),]
#Prepare an input file with ranked genes
geneList <- as.vector(merged_annotat_results_rmv_na$duplication_score)
names(geneList) <- as.vector(merged_annotat_results_rmv_na$entrezgene_id) #geneId stands for entrezgene ID 
#gene <- df_ann_DE_results$geneId

positive <- subset(merged_annotat_results_rmv_na, duplication_score >= 0)
geneList <- as.vector(positive$duplication_score)
names(geneList) <- as.vector(positive$entrezgene_id) #geneId stands for entrezgene ID 
positive_short = positive$entrezgene_id
positive_short = as.data.frame(positive_short)
write_csv(positive_short, "positive.csv")

negative <- subset(merged_annotat_results_rmv_na, duplication_score < 0)
geneList <- as.vector(negative$duplication_score)
names(geneList) <- as.vector(negative$entrezgene_id) #geneId stands for entrezgene ID 

###############################################    GO   ###############################################
geneList <- sort(geneList, decreasing = T)
#GO over-represintation test
ego <- enrichGO(gene          = unique(gene),
                universe      = names(unique(geneList)),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH")
head(ego)
emapplot(ego)
#Gene Set Enrichment Analysis GSEA
ego3 <- gseGO(geneList     = geneList[unique(names(geneList))],
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000)
head(ego3)
emapplot(ego3)

###############################################  KEGG  ###############################################
#kegg over-repres
kk <- enrichKEGG(gene         = unique(gene),
                 universe      = names(unique(geneList)),
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)
head(kk)
#kegg gsea
kk2 <- gseKEGG(geneList     = unique(geneList),
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
emapplot(kk2)

#disease enrichment
disEnr <- enrichDO(gene    = unique(gene),
                   ont           = "DO",
                   pAdjustMethod = "BH",
                   universe      = names(geneList))
head(disEnr)