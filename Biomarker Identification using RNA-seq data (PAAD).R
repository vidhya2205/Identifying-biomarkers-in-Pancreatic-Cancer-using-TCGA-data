setwd("D:/HackBio/Advanced genomicss course- bioinformatics for cancer biology/Proj1")
#Problem Statement - Identify potential biomarkers for early detection of pancreatic cancer by analyzing TCGA differential gene expression data. This project focuses on finding genes that are significantly expressed in early-stage pancreatic cancers compared to normal tissues.

#Load the packages
library("TCGAbiolinks")
library("edgeR")
library("EDASeq")
library("gplots")
#library("sesameData")
library("SummarizedExperiment")
library("biomaRt")

#Step 1 - Data Acquisition: Downlading data from TCGA (Though there are 179 tumor tissues only 4 normal)
getProjectSummary("TCGA-PAAD")

paadQ<- GDCquery(project = "TCGA-PAAD", 
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification")
GDCdownload(paadQ)
paad.data<-GDCprepare(paadQ)
View(paad.data)

table(paad.data$tissue_type) # type of tissue 179 tumor and only 4 normal
table(paad.data$ajcc_pathologic_stage) # early stage cancer can be considered as stage I, IA and IB

#select the unstranded dataset
paad.raw.data<-assays(paad.data) 
dim(paad.raw.data$unstranded) 
View(paad.raw.data$unstranded)

#metadata for further selecting process
simpleMeta<-data.frame("barcode"= paad.data$barcode, "Tissue"= paad.data$tissue_type, "Stage" =paad.data$ajcc_pathologic_stage)
View(simpleMeta)

#Selecting Tumrors with stage 1 (early stage) & normal tissue for further analysis

selectedBarcodes<- c(subset(simpleMeta, Tissue == "Tumor" & (Stage == "Stage I" | Stage == "Stage IA" | Stage == "Stage IB"))$barcode[c(1:10)],subset(simpleMeta, Tissue == "Normal")$barcode) # To get a sample ID
selectedBarcodes

selectedData<-paad.raw.data$unstranded[,c(selectedBarcodes)]
dim(selectedData)
View(selectedData) 

#Step 2: Normalization and Filteration
normData<- TCGAanalyze_Normalization(tabDF = selectedData, geneInfo = geneInfoHT, method= "geneLength")
filtData<- TCGAanalyze_Filtering(tabDF = normData,
                                 method = "quantile",
                                 qnt.cut = 0.25)
View(filtData)
dim(filtData)

#Step 3: Differential gene expression analysis
selectResults<-TCGAanalyze_DEA(mat1 = filtData[, c(selectedBarcodes)[1:10]],
                               mat2 = filtData[, c(selectedBarcodes)[11:14]],
                               Cond1type = "Stage 1 Pancreatic adenocarcinoma",
                               Cond2type = "Normal",
                               pipeline = "edgeR", 
                               fdr.cut = 0.01, #false discovery rate the % of something being false <0.01
                               logFC.cut = 2) #fold change is 2
View(selectResults)

#Differnetial expression levels for tumour vs normal
selectResults.levels<-
  TCGAanalyze_LevelTab(selectResults, "Stage 1 Pancreatic adenocarcinoma", "Normal", 
                       filtData[,c(selectedBarcodes)[1:10]],
                       filtData[,c(selectedBarcodes)[11:14]])
View(selectResults.levels) 

#Obtaing biomarkers by identifying top dysregulated(up/down regulated genes)
upreg.genes<- rownames(subset(selectResults.levels, logFC > 2))
dnreg.genes<- rownames(subset(selectResults.levels, logFC < -2))
mart<-useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

upreg.genes<- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = upreg.genes, mart = mart)$hgnc_symbol
dnreg.genes<- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), filters = 'ensembl_gene_id', values = dnreg.genes, mart = mart)$hgnc_symbol

#Performing Enrichment anlysis
up.EA<- TCGAanalyze_EAcomplete(TFname = "Upregulated", upreg.genes) 
dn.EA<- TCGAanalyze_EAcomplete(TFname = "Downregulated", dnreg.genes)

TCGAvisualize_EAbarplot(tf = rownames(up.EA$ResBP),#Rownames
                        GOBPTab = up.EA$ResBP, #results for BP
                        GOMFTab = up.EA$ResMF, #results for MF
                        GOCCTab = up.EA$ResCC, #results for CC
                        PathTab = up.EA$ResPat, #results for PAthway
                        nRGTab = upreg.genes, #number of genes in the list
                        nBar = 5, #max number of bars is 5 but can be increased to 10
                        text.size = 2, # 2 or 1.5
                        fig.width = 30, # size of figure
                        fig.height = 15) #generates a pdf in the working directory

TCGAvisualize_EAbarplot(tf = rownames(dn.EA$ResBP),#Rownames
                        GOBPTab = dn.EA$ResBP, #results for BP
                        GOMFTab = dn.EA$ResMF, #results for MF
                        GOCCTab = dn.EA$ResCC, #results for CC
                        PathTab = dn.EA$ResPat, #results for PAthway
                        nRGTab = dnreg.genes, #number of genes in the list
                        nBar = 5, #max number of bars is 5 but can be increased to 10
                        text.size = 2, # 2 or 1.5
                        fig.width = 30, # size of figure
                        fig.height = 15) #generates a pdf in the working directory