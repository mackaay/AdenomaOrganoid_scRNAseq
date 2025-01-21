library(dplyr)
library(Seurat)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(harmony)
library(RColorBrewer)
library(stringr)
library(DoubletFinder)
library(pheatmap)

fine.color <- brewer.pal(8, "Set1")
fine.color <- colorRampPalette(brewer.pal(8, "Set1"))(6)

pal.rename <- colorRampPalette(brewer.pal(8, "Accent"))(17)
names(pal.rename) <-c("dying","HMGB1+Enterocyte", "LYZ+Enterocyte", "BEST4+Enterocyte", "Enterocyte1", "Enterocyte2", "immature Enterocyte", 
                      "Goblet1", "Goblet2", "immature Goblet","TA1", "TA2", "TA3", "Progenitor","ISC",  
                      "EEC", "Tuft" ) 

###0. reading files ####
pbmc.data <- Read10X(data.dir = "/datasets/work/hb-exvivotcell/work/Data/level2/IRfsp_boost/N_Vehicle//outs/filtered_feature_bc_matrix//")
ctrl <- CreateSeuratObject(counts = pbmc.data, project = "N_Vehicle", min.cells = 3, min.features = 100)
ctrl
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-")
ctrl[["percent.rb"]] <- PercentageFeatureSet(ctrl, pattern = "^RP[SL]")
VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2, pt.size = -1, cols = fine.color)
plot1 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = fine.color, pt.size = 0.5)#+ylim(c(0,25))
plot2 <- FeatureScatter(ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = fine.color, pt.size = 0.5)
plot1 + plot2
ncol(subset(ctrl, subset = percent.mt <25)) /ncol(ctrl)
ctrl <- subset(ctrl, subset = nFeature_RNA > 100 &  percent.mt < 25)
#ctrl <- ctrl[grep("^MT-|^RP[SL]", rownames(ctrl), invert = T),]
VlnPlot(ctrl, features = "nFeature_RNA", pt.size = -1)

pbmc.data <- Read10X(data.dir = "/datasets/work/hb-exvivotcell/work/Data/level2/IRfsp_boost/N_LPS/outs/filtered_feature_bc_matrix/")
lps <- CreateSeuratObject(counts = pbmc.data, project = "N_LPS", min.cells = 3, min.features = 100)
lps[["percent.mt"]] <- PercentageFeatureSet(lps, pattern = "^MT-")
lps[["percent.rb"]] <- PercentageFeatureSet(lps, pattern = "^RP[SL]")
VlnPlot(lps, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2, pt.size = -1, cols = fine.color)
plot1 <- FeatureScatter(lps, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = fine.color, pt.size = 0.5)#+ylim(c(0,25))
plot2 <- FeatureScatter(lps, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = fine.color, pt.size = 0.5)
plot1 + plot2
ncol(subset(lps, subset = percent.mt <25)) /ncol(lps)
lps <- subset(lps, subset = nFeature_RNA > 100 &  percent.mt < 25)
VlnPlot(lps, features = "nFeature_RNA", pt.size = -1)#+ylim(c(0,5000))
#lps <- lps[grep("^MT-|^RP[SL]", rownames(lps), invert = T),]

pbmc.data <- Read10X(data.dir = "/datasets/work/hb-exvivotcell/work/Data/level2/IRfsp_boost/N_Combo/outs/filtered_feature_bc_matrix/")
combo <- CreateSeuratObject(counts = pbmc.data, project = "N_Combo", min.cells = 3, min.features = 100)
combo[["percent.mt"]] <- PercentageFeatureSet(combo, pattern = "^MT-")
combo[["percent.rb"]] <- PercentageFeatureSet(combo, pattern = "^RP[SL]")
VlnPlot(combo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2, pt.size = -1, cols = fine.color)
plot1 <- FeatureScatter(combo, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = fine.color, pt.size = 0.5)#+ylim(c(0,25))
plot2 <- FeatureScatter(combo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = fine.color, pt.size = 0.5)
plot1 + plot2
ncol(subset(combo, subset = percent.mt <25)) /ncol(combo)
VlnPlot(combo, features = "nFeature_RNA", pt.size = -1)
combo <- subset(combo, subset = nFeature_RNA > 500 &  percent.mt < 25)
#combo <- combo[grep("^MT-|^RP[SL]", rownames(combo), invert = T),]
plot1 <- VlnPlot(ctrl, features = "nFeature_RNA", pt.size = -1) + ylim(c(0,7000))
plot2 <- VlnPlot(lps, features = "nFeature_RNA", pt.size = -1)+ ylim(c(0,7000))
plot3 <- VlnPlot(combo, features = "nFeature_RNA", pt.size = -1)+ ylim(c(0,7000))
plot1 + plot2+ plot3

####Adenoma#####
pbmc.data <- Read10X(data.dir = "/datasets/work/hb-exvivotcell/work/Data/level2/IRfsp_boost/A_Vehicle/outs/filtered_feature_bc_matrix/")
Actrl <- CreateSeuratObject(counts = pbmc.data, project = "A_Vehicle", min.cells = 3, min.features = 100)
Actrl
Actrl[["percent.mt"]] <- PercentageFeatureSet(Actrl, pattern = "^MT-")
Actrl[["percent.rb"]] <- PercentageFeatureSet(Actrl, pattern = "^RP[SL]")
VlnPlot(Actrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2, pt.size = -1, cols = fine.color)
plot1 <- FeatureScatter(Actrl, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = fine.color, pt.size = 0.5)#+ylim(c(0,25))
plot2 <- FeatureScatter(Actrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = fine.color, pt.size = 0.5)
plot1 + plot2
ncol(subset(Actrl, subset = percent.mt <25)) /ncol(Actrl)
Actrl <- subset(Actrl, subset = nFeature_RNA > 100 &  percent.mt < 25)
#Actrl <- Actrl[grep("^MT-|^RP[SL]", rownames(Actrl), invert = T),]
VlnPlot(Actrl, features = "nFeature_RNA", pt.size = -1)

pbmc.data <- Read10X(data.dir = "/datasets/work/hb-exvivotcell/work/Data/level2/IRfsp_boost/A_LPS/outs/filtered_feature_bc_matrix/")
Alps <- CreateSeuratObject(counts = pbmc.data, project = "A_LPS", min.cells = 3, min.features = 100)
Alps[["percent.mt"]] <- PercentageFeatureSet(Alps, pattern = "^MT-")
Alps[["percent.rb"]] <- PercentageFeatureSet(Alps, pattern = "^RP[SL]")
VlnPlot(Alps, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2, pt.size = -1, cols = fine.color)
plot1 <- FeatureScatter(Alps, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = fine.color, pt.size = 0.5)#+ylim(c(0,25))
plot2 <- FeatureScatter(Alps, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = fine.color, pt.size = 0.5)
plot1 + plot2
ncol(subset(Alps, subset = percent.mt <25)) /ncol(Alps)
Alps <- subset(Alps, subset = nFeature_RNA > 100 &  percent.mt < 25)
VlnPlot(Alps, features = "nFeature_RNA", pt.size = -1)#+ylim(c(0,5000))
#Alps <- Alps[grep("^MT-|^RP[SL]", rownames(Alps), invert = T),]

pbmc.data <- Read10X(data.dir = "/datasets/work/hb-exvivotcell/work/Data/level2/IRfsp_boost/A_Combo/outs/filtered_feature_bc_matrix/")
Acombo <- CreateSeuratObject(counts = pbmc.data, project = "A_Combo", min.cells = 3, min.features = 100)
Acombo[["percent.mt"]] <- PercentageFeatureSet(Acombo, pattern = "^MT-")
Acombo[["percent.rb"]] <- PercentageFeatureSet(Acombo, pattern = "^RP[SL]")
VlnPlot(Acombo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 2, pt.size = -1, cols = fine.color)
plot1 <- FeatureScatter(Acombo, feature1 = "nCount_RNA", feature2 = "percent.mt", cols = fine.color, pt.size = 0.5)#+ylim(c(0,25))
plot2 <- FeatureScatter(Acombo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = fine.color, pt.size = 0.5)
plot1 + plot2
ncol(subset(Acombo, subset = percent.mt <25)) /ncol(Acombo)
VlnPlot(Acombo, features = "nFeature_RNA", pt.size = -1)
Acombo <- subset(Acombo, subset = nFeature_RNA > 100 &  percent.mt < 25)
#Acombo <- Acombo[grep("^MT-|^RP[SL]", rownames(Acombo), invert = T),]
plot1 <- VlnPlot(ctrl, features = "nFeature_RNA", pt.size = -1) + ylim(c(0,7000))
plot2 <- VlnPlot(lps, features = "nFeature_RNA", pt.size = -1)+ ylim(c(0,7000))
plot3 <- VlnPlot(Acombo, features = "nFeature_RNA", pt.size = -1)+ ylim(c(0,7000))
plot1 + plot2+ plot3



###1. batch correction #####
seu.list <- list(ctrl, lps, combo, Actrl, Alps, Acombo)
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x) # normalize data
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000) # find 2k feature genes
})
features <- SelectIntegrationFeatures(object.list = seu.list)
data.anchors <- FindIntegrationAnchors(object.list = seu.list,
                                       anchor.features = features,
                                       reduction = "cca",
                                       verbose=FALSE)
data.combined <- IntegrateData(anchorset = data.anchors, verbose=FALSE)

DefaultAssay(data.combined) <- "integrated"
data.combined <- ScaleData(data.combined, verbose = TRUE)
data.combined <- RunPCA(data.combined, npcs = 20, verbose = TRUE)
data.combined <- RunUMAP(data.combined, reduction = "pca",
                         dims = 1:20, verbose = TRUE)
data.combined <- FindNeighbors(data.combined, reduction = "pca", dims = 1:20,verbose = TRUE)
data.combined <- FindClusters(data.combined, resolution = c(0.4, 0.6, 0.8,1), verbose = TRUE)
Idents(data.combined) <- data.combined[["integrated_snn_res.1"]] 
DimPlot(data.combined, reduction = "umap", label = T,  pt.size = 0.5)
FeaturePlot(data.combined, features = "percent.mt")
FeaturePlot(data.combined, features = "percent.rb")
FeaturePlot(data.combined, features = "nFeature_RNA")
Idents(data.combined) <- data.combined[["orig.ident"]] 
DimPlot(data.combined, reduction = "umap", label = T,  pt.size = 0.5)


###2. Annotation####
#https://www.nature.com/articles/s41588-022-01088-x
DefaultAssay(data.combined) <- "RNA"
FeaturePlot(data.combined, features = c("EPCAM", "VIM", "PTPRC", "CD80"), pt.size = 1, ncol = 2)
tmp <- subset(data.combined, subset = EPCAM >0)
ncol(tmp)/ncol(data.combined)
ncol(subset(data.combined, subset = VIM >0))/ncol(data.combined)
FeatureScatter(data.combined, feature1 = "VIM", feature2 = "EPCAM")
tmp <- subset(data.combined, subset = VIM >0)
tmp$orig.ident
pie(table(tmp$orig.ident), col = fine.color)
pie(table(data.combined$orig.ident), col = fine.color)
FeatureScatter(data.combined, feature1 = "nCount_RNA", feature2 = "percent.mt",  pt.size = 0.5)

Idents(data.combined) <- "orig.ident"
DimPlot(data.combined, reduction = "umap", label = F, pt.size = 1, cols = fine.color) +ggtitle("Renaming All")+ #NoLegend() 

save(data.combined, file = "./34all_organoid_seurat.rds")

DotPlot(data.combined, features = c("DCLK1","HTR3C", "HTR3E" , "B4GALNT4", "IL13")) + RotatedAxis()+ggtitle("Tuft")
DotPlot(data.combined, features = c("KLK1", "ITLN1", "WFDC2" , "CLCA1", "RETNLB")) + RotatedAxis()+ggtitle("immature Goblet")
DotPlot(data.combined, features = c("MUC2", "TFF1", "FCGBP" , "TBX10")) + RotatedAxis()+ggtitle("mature Goblet")
DotPlot(data.combined, features = c("MUC2", "TFF1", "FCGBP" , "TBX10","KLK1", "ITLN1", "WFDC2" , "CLCA1", "RETNLB"))+ RotatedAxis()+ggtitle("Goblet")
DotPlot(data.combined, features = c("CA1","RAB6B")) + RotatedAxis()+ggtitle("Enterocytes")
DotPlot(data.combined, features = c("CRYBA2","SCGN", "FEV", "RGS2", "CDKN1A")) + RotatedAxis()+ggtitle("EEC")
DotPlot(data.combined, features = c("BEST4", "CA7", "CA4", "OTOP2" ,  "CFTR", "NOTCH2", "NOTCH2NL","TRAF4","GUCA2A","GUCA2B")) + RotatedAxis()+ggtitle("BEST4+ Enterocyte")
DotPlot(data.combined, features = c("BEST4")) + RotatedAxis()+ggtitle("BEST4+ Enterocyte")
DotPlot(data.combined, features = c("SMOC2", "RGMB", "LGR5" , "ASCL2", "SOX9")) + RotatedAxis()+ggtitle("stem cell")
DotPlot(data.combined, features = c("VIL1", "NLRP6", "IL18", "REG3A", "CD25", "REG1", "SLC7A7", 
                                    "SLC7A8", "SLC2A2", "MGAM", "SLC5A1", "DGAT1", "SLC27A4", "NPC1L1", 
                                    "APOA4", "NEAT1", "MALAT1", "ENPP3", "ADA", "SLC28A2", "NT5E")) + RotatedAxis()+ggtitle("Villi cell")
cellsignature <- read.csv("./datasets/celltype_signature.csv")
DotPlot(data.combined, features = cellsignature$gene[grep("EEC", cellsignature$group)]) + RotatedAxis()+ggtitle("EEC")
DotPlot(data.combined, features = cellsignature$gene[grep("Enterocytes", cellsignature$group)]) + RotatedAxis()+ggtitle("Enterocytes")
DotPlot(data.combined, features = cellsignature$gene[grep("Goblet", cellsignature$group)]) + RotatedAxis()+ggtitle("Goblet")
DotPlot(data.combined, features = cellsignature$gene[grep("Paneth", cellsignature$group)]) + RotatedAxis()+ggtitle("Paneth")
DotPlot(data.combined, features = cellsignature$gene[grep("qISCs", cellsignature$group)]) + RotatedAxis()+ggtitle("qISCs")
DotPlot(data.combined, features = cellsignature$gene[grep("TA", cellsignature$group)]) + RotatedAxis()+ggtitle("TA")
DotPlot(data.combined, features = cellsignature$gene[grep("RSCs", cellsignature$group)]) + RotatedAxis()+ggtitle("RSCs")
DotPlot(data.combined, features = cellsignature$gene[grep("Tuft", cellsignature$group)]) + RotatedAxis()+ggtitle("Tuft")
DotPlot(data.combined, features = cellsignature$gene[grep("ISCs", cellsignature$group)]) + RotatedAxis()+ggtitle("ISCs")
DotPlot(data.combined, features = c("LGR5", "SMOC2", "RGMB", "PTPRO", "EPHB2", "LRIG1")) + RotatedAxis()+ggtitle("ISCs")
DotPlot(data.combined, features = c("KRT20", "FABP2", "MKI67", "OLFM4", "CLU", "DCLK1", "MUC2", "LYZ")) + RotatedAxis()+ggtitle("")
DotPlot(data.combined, features = c("KRT20", "MUC5B", "CHGA","MUC2", "IL4R" , "IL13RA1" , "IL17RB")) + RotatedAxis()+ggtitle("")
DotPlot(data.combined, features = c("PROM1", "MSI1", "OLFM4", "EPHB2")) + RotatedAxis()+ggtitle("TA")
DotPlot(data.combined, features = c("TNFRSF11A","TRAF6", "RANKL", "SPIB", "RANK", "GP2","ANXA5", "CCL20", "TNFAIP2", "CCL9","PGLRP1")) + RotatedAxis()+ggtitle("M like")

Idents(data.combined) <- data.combined[["seurat_clusters"]] 
DotPlot(data.combined, features = c("VIL1", "DGAT1", "NT5E","SLC5A1",#Villi
                                    "CA1","ALDOB", "PRAP1", "PLAC8", "FABP1","FABP2","KRT20", "BEST4","CA7", "OTOP2", "OTOP3", "CCL20", "CCL23",
                                    "CA4","TRAF4","GUCA2A","GUCA2B", #Enterocyte
                                    "MUC2","TFF1","KLK1","SPINK4", #Goblet
                                    "MKI67", "UBE2C", "HMMR","HIST1H1B", "PROM1", "MSI1", "OLFM4", "EPHB2",#TA
                                    "LYZ",#Paneth
                                    "PLCG2","DCLK1", "POU2F3", #Tuft
                                    "SMOC2", "RGMB", "LRIG1","BMI1", "CDK6", "RNF43", "MBNL2","LGR5" ,  "SOX9", #ISC
                                    "RGS2", "CDKN1A"#EE
)) + RotatedAxis()+ggtitle("")
DotPlot(data.combined, features = c("VIL1", "DGAT1", "NT5E","SLC5A1",#Villi
                                    "CA1","ALDOB", "PRAP1", "PLAC8", "FABP1","FABP2","KRT20",#Enterocyte
                                    "MUC2", "TFF1","KLK1","SPINK4", #Goblet
                                    "MKI67", "UBE2C", "HMMR","HIST1H1B",#TA
                                    "LYZ",#Paneth
                                    "PLCG2","DCLK1", "POU2F3", #Tuft
                                    "SMOC2", "RGMB", "LRIG1","BMI1", "CDK6", "RNF43", "MBNL2","LGR5" ,  "SOX9", #ISC
                                    "RGS2", "CDKN1A"#EE
                                    ), scale = T, 
                                    cols = c("white", "red")) + RotatedAxis()+ggtitle("")
FeaturePlot(data.combined, features = "LGR5", pt.size = 2)
FeaturePlot(data.combined, features = c("GP2", "ANXA5"), pt.size = 2)
FeaturePlot(data.combined, features = "LYZ", pt.size = 2)
VlnPlot(data.combined, features = c("KRT8","FABP1", "ANXA2","ANXA5","FAS", "FASL", "BCL2", "BAX", "HMGB1", "GSDMD",
                                    "CASP1", "CASP11", "CASP8","CASP7","CASP3",  "NLRP3", "NLRP6", "NLRC4", "NLRC9B", "LYZ"),  pt.size = 0, ncol = 6)

VlnPlot(data.combined, features = c("KRT8","FABP1", "ANXA2","ANXA5","NOX4"), 
        cols = pal.rename, pt.size = -1, ncol = 2)

FeatureScatter(data.combined, feature1 = "nFeature_RNA", feature2 = "percent.mt",  pt.size = 1)+xlim(c(0,10000))
FeatureScatter(data.combined[,grep("0|3", Idents(data.combined))], feature1 = "nFeature_RNA", feature2 = "percent.mt",  pt.size = 1)+xlim(c(0,10000))
VlnPlot(data.combined, features = c("nFeature_RNA", "percent.mt"), cols = pal.rename, pt.size = 0)

Idents(data.combined) <- "orig.ident"
DimPlot(data.combined, reduction = "umap", label = T,  pt.size = 0.5)
DimPlot(data.combined[,grep("N_Vehicle", Idents(data.combined))], reduction = "umap", label = T,  pt.size = 0.5) + ggtitle("N_Vehicle")
DimPlot(data.combined[,grep("N_LPS", Idents(data.combined))], reduction = "umap", label = T,  pt.size = 0.5)+ ggtitle("N_LPS")
DimPlot(data.combined[,grep("N_Combo", Idents(data.combined))], reduction = "umap", label = T,  pt.size = 0.5)+ ggtitle("N_Combo")
DimPlot(data.combined[,grep("A_Ctrl", Idents(data.combined))], reduction = "umap", label = T,  pt.size = 0.5)+ ggtitle("A_Ctrl")
DimPlot(data.combined[,grep("A_LPS", Idents(data.combined))], reduction = "umap", label = T,  pt.size = 0.5)+ ggtitle("A_LPS")
DimPlot(data.combined[,grep("A_Combo", Idents(data.combined))], reduction = "umap", label = T,  pt.size = 0.5)+ ggtitle("A_Combo")

nctrl <- data.combined[,grep("N_Vehicle", Idents(data.combined))]
Idents(nctrl) <- "integrated_snn_res.0.6"
FeatureScatter(nctrl[,grep("1|2|3", Idents(nctrl))], feature1 = "nFeature_RNA", feature2 = "percent.mt", cols = fine.color,  pt.size = 1)+xlim(c(0,10000))




#####Find all markers#####
Idents(data.combined)
pbmc.markers <- FindAllMarkers(data.combined[grep("^MT-|^RP[SL]|AC0", rownames(data.combined), invert = T),], only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100
as.data.frame(top100[grep("13|3", top100$cluster),])







###3.inferCNV####
#BiocManager::install("infercnv")
library(infercnv)

counts_matrix = as.matrix(data.combined@assays$RNA@counts[,colnames(data.combined)]) 
write.table(round(counts_matrix, digits=3), file='./plots/infercnv/all_infercnv.counts.matrix', quote=F, sep="\t")    

annotations_file <- Idents(data.combined)



###4. velocity####
# 获取每个细胞的barcode
tmp <- Cells(data.combined)
tmp <- paste("N_Vehicle", tmp, sep = "-")
tmp[grep("-1_2", tmp)] <- gsub("N_Vehicle", "N_LPS", tmp[grep("-1_2", tmp)] )
tmp[grep("-1_3", tmp)] <- gsub("N_Vehicle", "N_Combo", tmp[grep("-1_3", tmp)] )
tmp[grep("-1_4", tmp)] <- gsub("N_Vehicle", "A_Vehicle", tmp[grep("-1_4", tmp)] )
tmp[grep("-1_5", tmp)] <- gsub("N_Vehicle", "A_LPS", tmp[grep("-1_5", tmp)] )
tmp[grep("-1_6", tmp)] <- gsub("N_Vehicle", "A_Combo", tmp[grep("-1_6", tmp)] )
tmp <- gsub("-1_1", "", tmp)
tmp <- gsub("-1_2", "", tmp)
tmp <- gsub("-1_3", "", tmp)
tmp <- gsub("-1_4", "", tmp)
tmp <- gsub("-1_5", "", tmp)
tmp <- gsub("-1_6", "", tmp)
write.csv(tmp, file = "./velocity/velocity_cellID_obs_all.csv", row.names = FALSE)
# 获得每个细胞的UMAP或TSNE坐标，使用 Embeddings函数
tmp2 <- Embeddings(data.combined, reduction = "umap")
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_embeddings_all.csv")
# 提取每个细胞的cluster信息
tmp2<- data.combined@meta.data[, 'integrated_snn_res.1', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_clusters_all.csv")
# 提取每个细胞的celltype信息
tmp2<- data.combined@meta.data[, 'orig.ident', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_celltype_all.csv")


Idents(data.combined) <- "integrated_snn_res.1"
nctrl <- data.combined[,grep("N_Vehicle", data.combined$orig.ident)]
tmp <- Cells(nctrl)
tmp <- paste("N_Vehicle", tmp, sep = "-")
tmp <- gsub("-1_1", "", tmp)
write.csv(tmp, file = "./velocity/velocity_cellID_obs_nctrl.csv", row.names = FALSE)
# 获得每个细胞的UMAP或TSNE坐标，使用 Embeddings函数
tmp2 <- Embeddings(nctrl, reduction = "umap")
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_embeddings_nctrl.csv")
# 提取每个细胞的cluster信息
tmp2<- nctrl@meta.data[, 'seurat_clusters', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_clusters_nctrl.csv")
# 提取每个细胞的celltype信息
tmp2<- nctrl@meta.data[, 'orig.ident', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_celltype_nctrl.csv")





###5. Renaming####
Idents(data.combined) <- data.combined[['integrated_snn_res.1']]
DimPlot(data.combined,  label = T)#+NoLegend()

new.cluster.ids <- c("dying", "TA1", "HMGB1+Enterocyte", "immature Goblet", "Progenitor", "LYZ+Enterocyte", "EEC", "Enterocyte1", "TA2", 
                     "ISC", "Enterocyte2", "immature Goblet", "immature Enterocyte", "Goblet1", "Tuft", "TA3", "Goblet2", "BEST4+Enterocyte")
names(new.cluster.ids) <- levels(data.combined)
data.combined <- RenameIdents(data.combined, new.cluster.ids)
data.combined[['new.cluster.ids']] <- Idents(data.combined)
Idents(data.combined) <- data.combined[['new.cluster.ids']]
Idents(data.combined) <- factor(Idents(data.combined), levels = c("dying","HMGB1+Enterocyte", "LYZ+Enterocyte", "BEST4+Enterocyte", "Enterocyte1", "Enterocyte2", "immature Enterocyte", 
                                                                   "Goblet1", "Goblet2", "immature Goblet","TA1", "TA2", "TA3", "Progenitor","ISC",  
                                                                  "EEC", "Tuft" ) )
data.combined[['new.cluster.ids']] <- Idents(data.combined)
pal.rename <- colorRampPalette(brewer.pal(8, "Accent"))(17)
names(pal.rename) <-c("dying","HMGB1+Enterocyte", "LYZ+Enterocyte", "BEST4+Enterocyte", "Enterocyte1", "Enterocyte2", "immature Enterocyte", 
                      "Goblet1", "Goblet2", "immature Goblet","TA1", "TA2", "TA3", "Progenitor","ISC",  
                      "EEC", "Tuft" ) 
DimPlot(data.combined, reduction = "umap", label = F, pt.size = 1, cols = pal.rename) +ggtitle("Renaming All")#+ NoLegend()

DotPlot(data.combined, features = c("VIL1", "DGAT1", "NT5E","SLC5A1","TMPRSS2",#Villi
                                    "CA1","ALDOB", "PRAP1", "PLAC8", "FABP1","FABP2","KRT20",#Enterocyte
                                    "MUC2", "TFF1","KLK1","SPINK4", #Goblet
                                    "LYZ",#Paneth
                                    "PLCG2","DCLK1", "POU2F3", "IL4RA" , "IL13RA1" , "IL17RB", #Tuft
                                    "MKI67", "UBE2C", "HMMR","HIST1H1B", "EPHB2",#TA
                                    "SMOC2", "RGMB", "LRIG1","BMI1", "CDK6", "RNF43", "MBNL2","CFTR","LGR5" ,  "SOX9", #ISC
                                    "RGS2", "CDKN1A"#EE 
                                    )) + RotatedAxis()+ggtitle("")
DotPlot(data.combined, features = c("TRAF6", "TNFRSF11A","SPIB" , "GP2", "TLR4", 
                                    "CCL20", "NTRK2", "CCL23", "SOX8", "TNFRSF11B", "FOLR3")) + RotatedAxis()+ggtitle("M like cell")

save(data.combined, file = "./34all_organoid_seurat_renaming.rds")

VlnPlot(data.combined, features = "TLR4", cols = pal.rename, pt.size = 0)


gene.sig <- c("KLK12", "NEUROG3", "CLPS", "RCN2", "TUBB", "PRDX1", "SOX4", "MDK", "NPM1", "YBX1", "STMN1", "AFP", "TUBA1A", "PAX4", "SMIM6", 
              "CO7A2L", "SCG2", "GCH1", "GC", "NEUROD1", "RGS2", "KCTD12", "FEV", "CHGA", "CRYBA2", "SCT", "DDC", "LGALS4", "TPH1", 
              "SLC18A1", "SLC38A11", "CES1", "REG4", "CHGB", "TFF3", "TAC1", "CDHR3", "IGFBP3", "LMX1A", "TFF1", "GCLC", "PRAC1",
              "CDKN2A", "COL5A2", "RXFP4", "FXYD2", "NPW", "HOXB13")
secretory <- data.combined[,grep("ISC|Proge|Goblet", data.combined$new.cluster.ids)]
#secretory.df <- as.data.frame(secretory@assays$RNA@data)
Idents(secretory) <- 'new.cluster.ids'
secretory.df <- as.data.frame( AverageExpression(secretory)$RNA)

pdata <- data.frame(id = colnames(secretory.df), cell = colnames(secretory.df))
rownames(pdata) <- pdata$id
mat_col <- data.frame(group = pdata$cell
                      )
rownames(mat_col) <- colnames(secretory.df)
mat_colors <- list(group = pal.rename[c(8:10,14,15)])
breaksList = seq(0, 1.5, by = 0.01)
pheatmap(secretory.df[c("SMOC2", "RGMB", "LRIG1","BMI1", "CDK6", "RNF43", "MBNL2",  "SOX9", "LGR5",
                        "MUC2", "TFF1","SPINK4", "LGALS4", "REG4", "TFF3"),order(colnames(secretory.df))], 
         cluster_cols = F,cluster_rows = F, display_numbers = F, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList ,
         border_color      = NA, scale = "row",
         show_colnames     = FALSE,show_rownames     = T,
         annotation_col    = mat_col,annotation_colors = mat_colors,
         main = "Secretory cells" )

grep("REG3", rownames(data.combined), value = T)




###6. Percentage #####
load(file = "./34all_organoid_seurat_renaming.rds")
idx.cluster <- 17
group.df <- data.frame(Percentage =rep(0, 6*idx.cluster), 
                       group=c(rep("N_Vehicle", idx.cluster), rep("N_LPS", idx.cluster), rep("N_Combo", idx.cluster), 
                               rep("A_Vehicle", idx.cluster), rep("A_LPS", idx.cluster), rep("A_Combo", idx.cluster)))
Idents(data.combined) <- data.combined[["orig.ident"]] 
tmp <- data.combined[,grep("N_Vehicle", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[1:idx.cluster] <- table(Idents(tmp))/ncol(tmp)
tmp <- data.combined[,grep("N_LPS", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[(idx.cluster*1+1):(idx.cluster*2)] <- table(Idents(tmp))/ncol(tmp)
tmp <- data.combined[,grep("N_Combo", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[(idx.cluster*2+1):(idx.cluster*3)] <- table(Idents(tmp))/ncol(tmp)
tmp <- data.combined[,grep("A_Vehicle", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[c((idx.cluster*3+1):59, 61:(idx.cluster*4))] <- table(Idents(tmp))/ncol(tmp)
tmp <- data.combined[,grep("A_LPS", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[(idx.cluster*4+1):(idx.cluster*5)] <- table(Idents(tmp))/ncol(tmp)
tmp <- data.combined[,grep("A_Combo", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[c((idx.cluster*5+1):93,95:(idx.cluster*6))] <- table(Idents(tmp))/ncol(tmp)

tmp <- data.combined[,grep("N_Vehicle", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
group.df$Cluster <- rep(c(levels(Idents(tmp))), 6)
group.df$Cluster <- factor(group.df$Cluster, levels = c("dying","HMGB1+Enterocyte", "LYZ+Enterocyte", "BEST4+Enterocyte", "Enterocyte1", "Enterocyte2", "immature Enterocyte", 
                                                        "Goblet1", "Goblet2", "immature Goblet","TA1", "TA2", "TA3", "Progenitor","ISC",  
                                                        "EEC", "Tuft" ))
group.df$Percentage <- round(group.df$Percentage, 4)
group.df$group <- factor(group.df$group, levels = c("N_Vehicle", "N_LPS", "N_Combo","A_Vehicle", "A_LPS", "A_Combo"))

#https://r-charts.com/colors/
group.df[,] %>%
  ggplot(data = ., mapping = aes(x = group, y = Percentage, fill = Cluster)) +
  geom_col() +
  geom_text(mapping = aes(label = paste0(Percentage*100, "%") ),              # converting the values to percent
            size = 5,                                             # size of the font
            position = position_stack(vjust = 0.5)) +             # positioning in the middle
  #scale_fill_brewer(palette = "Accent") +                           # coloring the plot
  scale_fill_manual(values = pal.rename)+
  #facet_grid(.~sex) +
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "Percentage (All)",        # title
       fill = "Cluster") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  ) #+ 
png("./plots/32percentage.png", width = 15, height = 15, units = "cm", res = 600)
group.df[,] %>%
  ggplot(data = ., mapping = aes(x = group, y = Percentage, fill = Cluster)) +
  geom_col() +
  #scale_fill_brewer(palette = "Accent") +                           # coloring the plot
  scale_fill_manual(values = pal.rename)+
  #facet_grid(.~sex) +
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "Percentage (All)",        # title
       fill = "Cluster") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  ) #+ 
dev.off()

write.csv(group.df, file = "./Percentage_all.csv")



#####scVelo####
tmp <- Cells(data.combined)
tmp <- paste("N_Vehicle", tmp, sep = "-")
tmp[grep("-1_2", tmp)] <- gsub("N_Vehicle", "N_LPS", tmp[grep("-1_2", tmp)] )
tmp[grep("-1_3", tmp)] <- gsub("N_Vehicle", "N_Combo", tmp[grep("-1_3", tmp)] )
tmp[grep("-1_4", tmp)] <- gsub("N_Vehicle", "A_Vehicle", tmp[grep("-1_4", tmp)] )
tmp[grep("-1_5", tmp)] <- gsub("N_Vehicle", "A_LPS", tmp[grep("-1_5", tmp)] )
tmp[grep("-1_6", tmp)] <- gsub("N_Vehicle", "A_Combo", tmp[grep("-1_6", tmp)] )
tmp <- gsub("-1_1", "", tmp)
tmp <- gsub("-1_2", "", tmp)
tmp <- gsub("-1_3", "", tmp)
tmp <- gsub("-1_4", "", tmp)
tmp <- gsub("-1_5", "", tmp)
tmp <- gsub("-1_6", "", tmp)
write.csv(tmp, file = "./velocity/velocity_cellID_obs_all.csv", row.names = FALSE)
# 获得每个细胞的UMAP或TSNE坐标，使用 Embeddings函数
tmp2 <- Embeddings(data.combined, reduction = "umap")
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_embeddings_all.csv")
# 提取每个细胞的cluster信息
tmp2<- data.combined@meta.data[, 'new.cluster.ids', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_clusters_all_rename.csv")
# 提取每个细胞的celltype信息
tmp2<- data.combined@meta.data[, 'orig.ident', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_celltype_all.csv")








#####Find all markers#####
Idents(data.combined)
pbmc.markers <- FindAllMarkers(data.combined[grep("^MT-|^RP[SL]|AC0", rownames(data.combined), invert = T),], only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top100

DoHeatmap(data.combined, features = top10$gene)




Idents(data.combined)
pbmc.markers <- FindAllMarkers(data.combined[grep("^MT-|^RP[SL]|AC0", rownames(data.combined), invert = T),grep("ISC", Idents(data.combined))], only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

isc.avg <- as.data.frame( AverageExpression(data.combined[,grep("ISC", Idents(data.combined))])$RNA)
#rownames(isc.avg) <- gsub("RNA.", "", rownames(isc.avg))
library(pheatmap)
#pal.group
#annotation_col = data.frame(Cell = pData$group)
#rownames(annotation_col) = colnames(lcpm)
#names(pal.group) <- levels(as.factor( pData$group))
breaksList = seq(-1, 1, by = 0.01)
pheatmap(isc.avg[top10$gene, ], 
         scale = "row",  clustering_distance_cols = "euclidean", border_color = NA,
         clustering_method = "complete",show_rownames = T, show_colnames = T,
         fontsize_row = 10, fontsize_col =10, cluster_cols = F,  cluster_rows = F,   annotation_legend = T,
         color = colorRampPalette(c("blue", "white","red"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList, 
         #annotation_colors = list(Group = pal.group ),annotation_col = annotation_col,
         main = "ISG featured genes"
)
dev.off()
pheatmap(isc.avg[c("TUBB2B", "TUBA1A", "TUBB3", "TUBA3D", "SOX4" ,"DPYSL2", "S100A4",top10$gene[11:30]), ], 
         scale = "row",  clustering_distance_cols = "euclidean", border_color = NA,
         clustering_method = "complete",show_rownames = T, show_colnames = T,
         fontsize_row = 10, fontsize_col =12, cluster_cols = F,  cluster_rows = F,   annotation_legend = T,
         color = colorRampPalette(c("blue", "white","red"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList, 
         #annotation_colors = list(Group = pal.group ),annotation_col = annotation_col,
         main = "ISG featured genes"
)

FeaturePlot(data.combined, features = "LGR5")


###7. Removinig Dying cells####
Idents(data.combined) <- 'new.cluster.ids'
data.rm <- data.combined[,grep("dying", Idents(data.combined), invert = T)]
DimPlot(data.rm, reduction = "umap", label = F, pt.size = 0.5, cols = pal.rename) +
  ggtitle("Removed dying cells")#+ NoLegend()


nctrl <- data.rm[,grep("N_Vehicle", Idents(data.rm))]
Idents(nctrl) <- "new.cluster.ids"
nctrl
Idents(nctrl) <- "orig.ident"
tmp <- Cells(nctrl)
tmp <- paste("N_Vehicle", tmp, sep = "-")
tmp <- gsub("-1_1", "", tmp)
write.csv(tmp, file = "./velocity/velocity_cellID_obs_nctrl_rm.csv", row.names = FALSE)
# 获得每个细胞的UMAP或TSNE坐标，使用 Embeddings函数
tmp2 <- Embeddings(nctrl, reduction = "umap")
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_embeddings_nctrl_rm.csv")
# 提取每个细胞的cluster信息
tmp2<- nctrl@meta.data[, 'new.cluster.ids', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_clusters_nctrl_rm.csv")
# 提取每个细胞的celltype信息
tmp2<- nctrl@meta.data[, 'orig.ident', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_celltype_nctrl_rm.csv")


actrl <- data.rm[,grep("A_Ctrl", Idents(data.rm))]
Idents(actrl) <- "new.cluster.ids"
actrl
Idents(actrl) <- "orig.ident"
tmp <- Cells(actrl)
tmp <- paste("A_Ctrl", tmp, sep = "-")
tmp <- gsub("-1_4", "", tmp)
write.csv(tmp, file = "./velocity/velocity_cellID_obs_actrl_rm.csv", row.names = FALSE)
# 获得每个细胞的UMAP或TSNE坐标，使用 Embeddings函数
tmp2 <- Embeddings(actrl, reduction = "umap")
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_embeddings_actrl_rm.csv")
# 提取每个细胞的cluster信息
tmp2<- actrl@meta.data[, 'seurat_clusters', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_clusters_actrl_rm.csv")
# 提取每个细胞的celltype信息
tmp2<- actrl@meta.data[, 'orig.ident', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_celltype_actrl_rm.csv")





##
load(file = "./34all_organoid_seurat_renaming.rds")
Idents(data.combined) <- "orig.ident"
nctrl <- data.combined[,grep("N_Vehicle", Idents(data.combined))]
nctrl
Idents(nctrl) <- "orig.ident"
tmp <- Cells(nctrl)
tmp <- paste("N_Vehicle", tmp, sep = "-")
tmp <- gsub("-1_1", "", tmp)
write.csv(tmp, file = "./velocity/velocity_cellID_obs_nctrl.csv", row.names = FALSE)
# 获得每个细胞的UMAP或TSNE坐标，使用 Embeddings函数
tmp2 <- Embeddings(nctrl, reduction = "umap")
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_embeddings_nctrl.csv")
# 提取每个细胞的cluster信息
tmp2<- nctrl@meta.data[, 'new.cluster.ids', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_clusters_nctrl.csv")
# 提取每个细胞的celltype信息
tmp2<- nctrl@meta.data[, 'orig.ident', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_celltype_nctrl.csv")


actrl <- data.combined[,grep("A_Ctrl", Idents(data.combined))]
actrl
Idents(actrl) <- "orig.ident"
tmp <- Cells(actrl)
tmp <- paste("A_Ctrl", tmp, sep = "-")
tmp <- gsub("-1_4", "", tmp)
write.csv(tmp, file = "./velocity/velocity_cellID_obs_actrl.csv", row.names = FALSE)
# 获得每个细胞的UMAP或TSNE坐标，使用 Embeddings函数
tmp2 <- Embeddings(actrl, reduction = "umap")
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_embeddings_actrl.csv")
# 提取每个细胞的cluster信息
tmp2<- actrl@meta.data[, 'new.cluster.ids', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_clusters_actrl.csv")
# 提取每个细胞的celltype信息
tmp2<- actrl@meta.data[, 'orig.ident', drop = FALSE]
rownames(tmp2) <- tmp
write.csv(tmp2, file = "./velocity/velocity_cell_celltype_actrl.csv")


####Mesenchymal Cell####
tmp <- subset(data.rm, subset = VIM >0)
tmp$orig.ident
pie(table(tmp$orig.ident), col = fine.color)
pie(table(data.rm$orig.ident), col = fine.color)
ncol(tmp)/ncol(data.rm)
ncol(subset(data.rm, subset = VIM >0))/ncol(data.rm)



Idents(data.rm) <- 'orig.ident'
data.rm.avg <- as.data.frame( AverageExpression(data.rm[,])$RNA)
#rownames(isc.avg) <- gsub("RNA.", "", rownames(isc.avg))
library(pheatmap)
#pal.group
#annotation_col = data.frame(Cell = pData$group)
#rownames(annotation_col) = colnames(lcpm)
#names(pal.group) <- levels(as.factor( pData$group))
c
dev.off()


####tight junction, Villi, Cytokine, ####
load(file = "./34all_organoid_seurat_renaming.rds")
Idents(data.combined) <- 'new.cluster.ids'
VlnPlot(data.combined, features = c("TJP1","TJP2","TJP3","HIF1A","CLDN3", "OCLN","F11R","MYH1","ACTB"), 
        cols = pal.rename, pt.size = -1, ncol = 4)
FeaturePlot(data.combined, features = c("TJP1","CLDN3", "F11R"), ncol =  3 , pt.size = 1)
VlnPlot(data.combined, features = c("MUC1", "MUC3A","MUC3B", "MUC4", "MUC12", "MUC13", "MUC15", "MUC17", "MUC20", "MUC21",
                                   "MUC2", "MUC5AC", "MUC5B", "MUC6",  "MUC19" , 
                                   "MUC7", "MUC8",  "MUC9"), 
        cols = pal.rename, pt.size = -1, ncol = 4)
VlnPlot(data.combined, features = c("MUC1", "MUC3A","MUC12", "MUC13",  "MUC20","MUC2", "MUC5AC"), 
        cols = pal.rename, pt.size = -1, ncol = 4)
FeaturePlot(data.combined, features = c("MUC12","MUC3A", "MUC13"), ncol =  3, pt.size = 1 )


Idents(data.combined) <- "orig.ident"
VlnPlot(data.combined, features = c("TJP1","TJP2","TJP3","HIF1A","CLDN3", "OCLN","F11R","MYH1","ACTB"), 
        cols = fine.color, pt.size = -1, ncol = 4)
FeaturePlot(data.combined, features = "TJP1")
FeaturePlot(data.combined, features = "TJP1")

VlnPlot(data.combined, features = c("MUC1", "MUC3A","MUC3B", "MUC12", "MUC13", "MUC15",  "MUC20", "MUC21","MUC2", "MUC5AC","MUC6"), 
        cols = fine.color, pt.size = -1, ncol = 4)


VlnPlot(data.surface, features = c("FCGBP" , "TFF3","CLCA1", "CLCA2","CLCA3","CLCA4","RETNLB"), 
        cols = fine.color, pt.size = -1, ncol = 4)
FeaturePlot(data.rm, features = "TJP1")
VlnPlot(data.surface, features = c("TOP2A","CDC14B", "CDKN1A","ACTB"),cols = fine.color, pt.size = -1, ncol = 4)
Idents(data.surface) <- "orig.ident"
VlnPlot(data.surface, features = c("NOS1", "IL6",  "IL12A", "TLR4"), 
        cols = fine.color, pt.size = -1, ncol = 4)
Idents(data.rm) <- "orig.ident"
VlnPlot(data.rm, features = c("SLC16A1", "ANXA2", "ANXA5", "ACTB"), 
        cols = fine.color, pt.size = -1, ncol = 2)
Idents(data.rm) <- "orig.ident"
VlnPlot(data.rm, features = c("TMPRSS2" ,"ACTB"), 
        cols = fine.color, pt.size = -1, ncol = 2)
Idents(data.combined) <- "orig.ident"
VlnPlot(data.combined, features = c("CXCL9","CXCL10", "CXCL11", "CXCL12", 
                              "CCL4", "CCL5", "CCL20", 
                              "CCL17", "CCL22", "CCL2", "CCL5", "CXCL5", "CXCL2", "CXCL6"), 
        cols = fine.color, pt.size = -1, ncol = 4)


Idents(data.combined) <- 'orig.ident'
data.combined.avg <- as.data.frame( AverageExpression(data.combined[,]))
breaksList = seq(-1, 1, by = 0.01)
pheatmap(data.combined.avg[c(grep("^IL[0-9]", rownames(data.combined), value = T)), ], 
         scale = "row",  clustering_distance_cols = "euclidean", border_color = NA,
         clustering_method = "complete",show_rownames = T, show_colnames = T,
         fontsize_row = 12, fontsize_col =12, cluster_cols = F,  cluster_rows = T,   annotation_legend = T,
         color = colorRampPalette(c("blue", "white","red"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList, 
         #annotation_colors = list(Group = pal.group ),annotation_col = annotation_col,
         main = "IL featured genes"
)

#top up-regulated in N_LPS
nlps <- rownames( data.combined.avg[(data.combined.avg$N_LPS/rowMeans(data.combined.avg[,c(1,3)]) >1.5 & data.combined.avg$N_LPS>1),] )
nctrl <- rownames( data.combined.avg[(data.combined.avg$N_Vehicle/rowMeans(data.combined.avg[,c(2,3)]) >1.5 & data.combined.avg$N_Vehicle>1),] )
ncom <- rownames( data.combined.avg[(data.combined.avg$N_Combo/rowMeans(data.combined.avg[,c(2,1)]) >1.5 & data.combined.avg$N_Combo>1),] )
breaksList = seq(-2, 2, by = 0.01)
pheatmap(data.combined.avg[c(nctrl, nlps, ncom), ], 
         scale = "row",  clustering_distance_cols = "euclidean", border_color = NA,
         clustering_method = "complete",show_rownames = T, show_colnames = T,
         fontsize_row = 12, fontsize_col =12, cluster_cols = F,  cluster_rows = F,   annotation_legend = T,
         #color = colorRampPalette(c("blue", "white","red"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         #breaks = breaksList, 
         #annotation_colors = list(Group = pal.group ),annotation_col = annotation_col,
         main = "LPS featured genes"
)

nlps <- rownames( data.combined.avg[(data.combined.avg$N_LPS/data.combined.avg$N_Vehicle>1.5 & data.combined.avg$N_LPS>1),] )
#nlps2 <- rownames( data.combined.avg[(data.combined.avg$N_Vehicle/data.combined.avg$N_LPS>1.5 & data.combined.avg$N_LPS>1),] )
tmp <- data.combined.avg[c(nlps), ]
tmp <- tmp[order(c(tmp$N_LPS/ tmp$N_Vehicle)),]
pheatmap(t(tmp), 
         scale = "column",  clustering_distance_cols = "euclidean", border_color = NA,
         clustering_method = "complete",show_rownames = T, show_colnames = T,
         fontsize_row = 12, fontsize_col =12, cluster_cols = F,  cluster_rows = F,   annotation_legend = T,
         #color = colorRampPalette(c("blue", "white","red"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         #breaks = breaksList, 
         #annotation_colors = list(Group = pal.group ),annotation_col = annotation_col,
         main = "LPS featured genes"
)

Idents(data.combined) <- 'new.cluster.ids'
data.rm <- data.combined[,grep("dying", Idents(data.combined), invert = T)]
DimPlot(data.rm, reduction = "umap", label = F, pt.size = 0.5, cols = pal.rename) +
  ggtitle("Removed dying cells")#+ NoLegend()

Idents(data.rm) <- "orig.ident"
VlnPlot(data.rm, features = c("CDH1", "S100A4","ACTA2", "DKK1", "NFKB1", "HK3", "CXCR4"), 
        cols = fine.color, pt.size = -1, ncol = 4)




#cytokine 
cytokine <- c("TGFBR2", "TGFBI","TGFBRAP1","TGFB1","IL33","IL10RB","TGFA","IL1RN","CCL20","CXCL1","CXCL3","CXCL2","CCL28",
              "IL6ST","IL18","IL32","IL4R","CXCL16","CCL15","IL17RA","IL2RG","IL13RA1")
Idents(data.combined) <- "orig.ident"
DotPlot(data.combined, features = cytokine) + RotatedAxis()
DotPlot(data.combined, features = c("CXCL8", "IL10", "TNF")) + RotatedAxis()

cytokine <- c("TGFBR2", "TGFBI","TGFBRAP1","TGFB1","IL33","IL10RB")
DotPlot(data.combined, features = cytokine) + RotatedAxis() + ggtitle("Anti-Inflammatory")
cytokine <- c("TGFA","IL1RN","CCL20","CXCL1","CXCL3","CXCL2","CCL28",
              "IL6ST","IL18","IL32","IL4R","CXCL16","CCL15","IL17RA","IL2RG","IL13RA1")
DotPlot(data.combined, features = cytokine) + RotatedAxis() + ggtitle("Pro-Inflammatory")


#####LPS pathways ######
library(limma)
Assays(data.combined) <- "integrated"

targets <- data.frame(names = data.combined$orig.ident)
cellType <- factor(targets$names)

# use the above to create a design matrix
design <- model.matrix(~0+cellType, data=targets)
colnames(design) <- c(levels(cellType))

# fit the linear model 
eset <- data.combined@assays$RNA@counts[grep("^MT-|^RP[SL]", rownames(data.combined@assays$RNA@counts), invert = T),]
eset <- log2(eset+1)
fit <- lmFit(eset, design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(N_LPS-N_Vehicle,
                            levels=design)
contMatrix
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))
DMPs <- topTable(fit2, num=Inf, coef=1)
head(DMPs)
DMPs[which(DMPs$logFC > 0 & DMPs$adj.P.Val < 0.01),]







Idents(data.combined) <- "orig.ident"
VlnPlot(data.combined, features = c("ACE2","TMPRSS2"), 
        cols = fine.color, pt.size = -1, ncol = 2)
FeaturePlot(data.combined, features = c("ACE2","TMPRSS2"))
DotPlot(data.combined, features = c("ACE2","TMPRSS2")) + RotatedAxis()+ggtitle("ACE2")


Idents(data.combined) <- "orig.ident"
VlnPlot(data.combined, features = c("NLRP6", "REG3G", "PIGR", "LYPD8", "CCL25", 
                              "ATP1B1", "FABP1", "FABP2", "SLC7A7", "SLC7A9",  "SLC2A2"), 
        cols = fine.color, pt.size = -1, ncol = 4)

Idents(data.combined) <- 'orig.ident'
data.combined.avg <- as.data.frame( AverageExpression(data.combined[,]) )
breaksList = seq(-1, 1, by = 0.01)
pheatmap(data.combined.avg[c("PIGR", "ATP1B1", "FABP1", "SLC7A9",  "SLC2A2"), ], 
         scale = "row",  clustering_distance_cols = "euclidean", border_color = NA,
         clustering_method = "complete",show_rownames = T, show_colnames = T,
         fontsize_row = 12, fontsize_col =12, cluster_cols = F,  cluster_rows = F,   annotation_legend = T,
         color = colorRampPalette(c("blue", "white","red"))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList, 
         #annotation_colors = list(Group = pal.group ),annotation_col = annotation_col,
         main = "Mesenchymal featured genes"
)


Idents(data.combined) <- "orig.ident"
DotPlot(data.combined, features = c("CXCL8", "TNF","IL10")) + RotatedAxis()+ggtitle("Cytokines")




###8. GSEA####
library(fgsea)
library(msigdbr)
library(dplyr)
library(presto)

#Normal
gsea.df <- data.frame(pathway=c(), padj=c(), NES=c())
Idents(data.combined) <- "new.cluster.ids"
list.name <- unique(data.combined$new.cluster.ids)
list.name <- gsub("\\+.*$", "", list.name)
for (i in list.name) {
  tmp <- data.combined[,grep(i, data.combined$new.cluster.ids)]
  print(i)
  Idents(tmp) <- 'orig.ident'
  tmp <- tmp[,grep("N_Vehicle|N_LPS", Idents(tmp))]
  pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
  head(pbmc.genes)
  m_df<- msigdbr(species = "Homo sapiens", category = "H") #Hallmark gene sets
  head(m_df)
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  pbmc.genes %>%
    dplyr::filter(group == "N_Vehicle") %>%
    arrange(desc(logFC), desc(auc)) %>%
    head(n = 10)
  cluster.genes<- pbmc.genes %>%
    dplyr::filter(group == "N_Vehicle") %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  ranks<- cluster.genes$auc
  names(ranks) <- cluster.genes$feature
  head(ranks)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj) %>% 
    head()
  tmp <- as.data.frame(fgseaResTidy %>% 
                         arrange(desc(pathway)) )[,c(1,3,5)]
  gsea.df <- rbind(gsea.df, tmp)
  
  tmp <- data.combined[,grep(i, Idents(data.combined))]
  print(i)
  Idents(tmp) <- 'orig.ident'
  tmp <- tmp[,grep("N_Vehicle|N_Combo", Idents(tmp))]
  pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
  head(pbmc.genes)
  m_df<- msigdbr(species = "Homo sapiens", category = "H") #Hallmark gene sets
  head(m_df)
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  pbmc.genes %>%
    dplyr::filter(group == "N_Vehicle") %>%
    arrange(desc(logFC), desc(auc)) %>%
    head(n = 10)
  cluster.genes<- pbmc.genes %>%
    dplyr::filter(group == "N_Vehicle") %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  ranks<- cluster.genes$auc
  names(ranks) <- cluster.genes$feature
  head(ranks)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj) %>% 
    head()
  tmp <- as.data.frame(fgseaResTidy %>% 
                         arrange(desc(pathway)) )[,c(1,3,5)]
  gsea.df <- rbind(gsea.df, tmp)
}
gsea.df$group <- rep(c(rep("LPS",50), rep("Combo", 50)) , 17)
list.name <- as.character( unique(data.combined$new.cluster.ids) )
gsea.df$cell <- list.name[1]
gsea.df$cell[101:200] <- list.name[2]
gsea.df$cell[201:300] <- list.name[3]
gsea.df$cell[301:400] <- list.name[4]
gsea.df$cell[401:500] <- list.name[5]
gsea.df$cell[501:600] <- list.name[6]
gsea.df$cell[601:700] <- list.name[7]
gsea.df$cell[701:800] <- list.name[8]
gsea.df$cell[801:900] <- list.name[9]
gsea.df$cell[901:1000] <- list.name[10]
gsea.df$cell[1001:1100] <- list.name[11]
gsea.df$cell[1101:1200] <- list.name[12]
gsea.df$cell[1201:1300] <- list.name[13]
gsea.df$cell[1301:1400] <- list.name[14]
gsea.df$cell[1401:1500] <- list.name[15]
gsea.df$cell[1501:1600] <- list.name[16]
gsea.df$cell[1601:1700] <- list.name[17]

idx.pathway <- grep("TNFA|TGF|INTERFERON|INFLAMMATORY|STAT", gsea.df$pathway)
gsea.df$group <- gsub("LPS", "LPS vs Ctrl", gsea.df$group)
gsea.df$group <- gsub("Combo", "Combo vs Ctrl", gsea.df$group)
gsea.df$group  <- factor(gsea.df$group, levels = c("LPS vs Ctrl", "Combo vs Ctrl"))
gsea.df$cell <- factor(gsea.df$cell, levels = c("dying","HMGB1+Enterocyte", "LYZ+Enterocyte", "BEST4+Enterocyte", "Enterocyte1", "Enterocyte2", "immature Enterocyte", 
                                                "Goblet1", "Goblet2", "immature Goblet","TA1", "TA2", "TA3", "Progenitor","ISC",  
                                                "EEC", "Tuft" )  )
ggplot(gsea.df[idx.pathway,], aes(x=group, y=pathway)) +
  geom_point(aes(size=-log10(padj), color = NES)) +
  scale_color_gradient2(midpoint=0,  low="#2166AC", mid="white",high="#B2182B", space ="Lab" ) +
  facet_wrap(vars(cell), nrow = 2)+ theme_classic() +
  RotatedAxis() +ggtitle("Normal")



#Adenoma
gsea.df <- data.frame(pathway=c(), padj=c(), NES=c())
Idents(data.combined) <- "new.cluster.ids"
list.name <- unique(data.combined$new.cluster.ids)
list.name <- gsub("\\+.*$", "", list.name)
list.name <- list.name[-16]
for (i in list.name) {
  tmp <- data.combined[,grep(i, data.combined$new.cluster.ids)]
  print(i)
  Idents(tmp) <- 'orig.ident'
  tmp <- tmp[,grep("A_Vehicle|A_LPS", Idents(tmp))]
  pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
  head(pbmc.genes)
  m_df<- msigdbr(species = "Homo sapiens", category = "H") #Hallmark gene sets
  head(m_df)
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  pbmc.genes %>%
    dplyr::filter(group == "A_Vehicle") %>%
    arrange(desc(logFC), desc(auc)) %>%
    head(n = 10)
  cluster.genes<- pbmc.genes %>%
    dplyr::filter(group == "A_Vehicle") %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  ranks<- cluster.genes$auc
  names(ranks) <- cluster.genes$feature
  head(ranks)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj) %>% 
    head()
  tmp <- as.data.frame(fgseaResTidy %>% 
                         arrange(desc(pathway)) )[,c(1,3,5)]
  gsea.df <- rbind(gsea.df, tmp)
  
  tmp <- data.combined[,grep(i, Idents(data.combined))]
  print(i)
  Idents(tmp) <- 'orig.ident'
  tmp <- tmp[,grep("A_Vehicle|A_Combo", Idents(tmp))]
  pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
  head(pbmc.genes)
  m_df<- msigdbr(species = "Homo sapiens", category = "H") #Hallmark gene sets
  head(m_df)
  fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
  pbmc.genes %>%
    dplyr::filter(group == "A_Vehicle") %>%
    arrange(desc(logFC), desc(auc)) %>%
    head(n = 10)
  cluster.genes<- pbmc.genes %>%
    dplyr::filter(group == "A_Vehicle") %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  ranks<- cluster.genes$auc
  names(ranks) <- cluster.genes$feature
  head(ranks)
  fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  fgseaResTidy %>% 
    dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
    arrange(padj) %>% 
    head()
  tmp <- as.data.frame(fgseaResTidy %>% 
                         arrange(desc(pathway)) )[,c(1,3,5)]
  gsea.df <- rbind(gsea.df, tmp)
}
gsea.df$group <- rep(c(rep("LPS",50), rep("Combo", 50)) , 16)
list.name <- as.character( unique(data.combined$new.cluster.ids) )
gsea.df$cell <- list.name[1]
gsea.df$cell[101:200] <- list.name[2]
gsea.df$cell[201:300] <- list.name[3]
gsea.df$cell[301:400] <- list.name[4]
gsea.df$cell[401:500] <- list.name[5]
gsea.df$cell[501:600] <- list.name[6]
gsea.df$cell[601:700] <- list.name[7]
gsea.df$cell[701:800] <- list.name[8]
gsea.df$cell[801:900] <- list.name[9]
gsea.df$cell[901:1000] <- list.name[10]
gsea.df$cell[1001:1100] <- list.name[11]
gsea.df$cell[1101:1200] <- list.name[12]
gsea.df$cell[1201:1300] <- list.name[13]
gsea.df$cell[1301:1400] <- list.name[14]
gsea.df$cell[1401:1500] <- list.name[15]
gsea.df$cell[1501:1600] <- list.name[16]
#gsea.df$cell[1601:1700] <- list.name[17]

idx.pathway <- grep("TNFA|TGF|INTERFERON|INFLAMMATORY|STAT", gsea.df$pathway)
gsea.df$group <- gsub("LPS", "LPS vs Ctrl", gsea.df$group)
gsea.df$group <- gsub("Combo", "Combo vs Ctrl", gsea.df$group)
gsea.df$group  <- factor(gsea.df$group, levels = c("LPS vs Ctrl", "Combo vs Ctrl"))
gsea.df$cell <- factor(gsea.df$cell, levels = c("dying","HMGB1+Enterocyte", "LYZ+Enterocyte", "BEST4+Enterocyte", "Enterocyte1", "Enterocyte2", "immature Enterocyte", 
                                                "Goblet1", "immature Goblet","TA1", "TA2", "TA3", "Progenitor","ISC",  
                                                "EEC", "Tuft" )  )
ggplot(gsea.df[idx.pathway,], aes(x=group, y=pathway)) +
  geom_point(aes(size=-log10(padj), color = NES)) +
  scale_color_gradient2(midpoint=0, low="#2166AC", mid="white",high="#B2182B", space ="Lab" ) +
  facet_wrap(vars(cell), nrow = 2)+ theme_classic() +
  RotatedAxis() +ggtitle("Adenoma")




##Doheatmap 
Idents(data.rm) <- "orig.ident"
tmp <- data.rm[,grep("N_", Idents(data.rm))]
head(m_df)
m_df$gene_symbol[grep("TNFA", m_df$gs_name)]
tmp@assays$RNA@scale.data <-  tmp@assays$integrated@scale.data
DoHeatmap(tmp, features = m_df$gene_symbol[grep("TNFA", m_df$gs_name)]) + NoLegend() +ggtitle("TNFA")
DoHeatmap(tmp, features = m_df$gene_symbol[grep("TGF", m_df$gs_name)]) + NoLegend()+ggtitle("TGF")
DoHeatmap(tmp, features = m_df$gene_symbol[grep("INTERFERON_GAMMA", m_df$gs_name)]) + NoLegend()+ggtitle("INTERFERON_GAMMA")
DoHeatmap(tmp, features = m_df$gene_symbol[grep("INTERFERON_ALPHA", m_df$gs_name)]) + NoLegend()+ggtitle("INTERFERON_ALPHA")
DoHeatmap(tmp, features = m_df$gene_symbol[grep("INFLAMMATORY", m_df$gs_name)]) + NoLegend()+ggtitle("INFLAMMATORY")
DoHeatmap(tmp, features = m_df$gene_symbol[grep("STAT3", m_df$gs_name)]) + NoLegend()+ggtitle("STAT3")
DoHeatmap(tmp, features = m_df$gene_symbol[grep("STAT5", m_df$gs_name)]) + NoLegend()+ggtitle("STAT5")

act <- c("CXCL9","CXCL10", "CXCL11", "CXCL12", 
         "CCL4", "CCL5", "CCL20")
sup <- c("CCL17", "CCL22", "CCL2", "CCL5", "CXCL5", "CXCL2", "CXCL6", 
         "TGFB1", "CSF1", "CSF2")
grep("^IL[0-9]|^CCL[0-9]|^CXCL[0-9]", rownames(data.combined), value = T)
data.rm@assays$RNA@scale.data <-  data.rm@assays$integrated@scale.data
DoHeatmap(data.rm, features = c(grep("^IL[0-9]|^CCL[0-9]|^CXCL[0-9]", rownames(data.rm), value = T)) )  +ggtitle("")
VlnPlot(tmp, features = c(grep("^IL[0-9]|^CCL[0-9]|^CXCL[0-9]|CSF|TGF", rownames(data.rm), value = T)), cols = fine.color, pt.size = 0, ncol = 5)
DotPlot(data.rm, features = c(grep("^IL[0-9]|^CCL[0-9]|^CXCL[0-9]|CSF|TGF", rownames(data.rm), value = T))) + RotatedAxis()


##
tmp <- data.rm
Idents(tmp) <- 'orig.ident'
tmp <- tmp[,grep("A_Ctrl|A_LPS", Idents(tmp))]
pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
head(pbmc.genes)
m_df<- msigdbr(species = "Homo sapiens", category = "H") #Hallmark gene sets
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
pbmc.genes %>%
  dplyr::filter(group == "A_Ctrl") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
cluster.genes<- pbmc.genes %>%
  dplyr::filter(group == "A_Ctrl") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks<- cluster.genes$auc
names(ranks) <- cluster.genes$feature
head(ranks)

grep("^CCL[0-9]", names(ranks))
grep("^CCL[0-9]", names(ranks), value = T)
Idents(data.rm) <- "orig.ident"
VlnPlot(data.rm, features = c("CCL20", "CCL28", "CCL15"), pt.size = 0, cols = fine.color)

grep("^CXCL[0-9]", names(ranks))
grep("^CXCL[0-9]", names(ranks), value = T)
Idents(data.rm) <- "orig.ident"
VlnPlot(data.rm, features = c("CXCL16", "CXCL2" , "CXCL1" ), pt.size = 0, cols = fine.color)

grep("^IL[0-9]", names(ranks))
grep("^IL[0-9]", names(ranks), value = T)
Idents(data.rm) <- "orig.ident"


head(names(ranks),50)
VlnPlot(data.rm, features = c("IL33" ,"IL12A" ,"IL5" , "IL13", "IL20","IL18" ), pt.size = 0, cols = fine.color)
VlnPlot(data.rm, features = c("LRRC75A","ITGA6" ,"CTNNB1" ,"TCF7L2" , "TGFBR2", "NEAT1","IL18", "IL12A" ), pt.size = 0, cols = fine.color)
tail(names(ranks),50)
VlnPlot(data.rm, features = c("FGFBP1" ,"CRB3" ,"MUC1" ,"TFF1" , "TFF3", "TFF2","ACAT2" ), pt.size = 0, cols = fine.color)
FeaturePlot(data.rm, features = c("FGFBP1" ,"MUC1" ,"TFF1" , "TFF3", "TFF2","ACAT2" ))




##
tmp <- data.rm
Idents(tmp) <- 'orig.ident'
tmp <- tmp[,grep("N_Vehicle|N_LPS", Idents(tmp))]
pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
head(pbmc.genes)
m_df<- msigdbr(species = "Homo sapiens", category = "H") #Hallmark gene sets
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
pbmc.genes %>%
  dplyr::filter(group == "N_Vehicle") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
cluster.genes<- pbmc.genes %>%
  dplyr::filter(group == "N_Vehicle") %>%
  arrange(desc(logFC)) %>% 
  dplyr::select(feature, logFC)
ranks.n<- cluster.genes$logFC
names(ranks.n) <- cluster.genes$feature
head(ranks.n)

tmp <- data.rm
Idents(tmp) <- 'orig.ident'
tmp <- tmp[,grep("A_Ctrl|A_LPS", Idents(tmp))]
pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
head(pbmc.genes)
m_df<- msigdbr(species = "Homo sapiens", category = "H") #Hallmark gene sets
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
pbmc.genes %>%
  dplyr::filter(group == "A_Ctrl") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
cluster.genes<- pbmc.genes %>%
  dplyr::filter(group == "A_Ctrl") %>%
  arrange(desc(logFC)) %>% 
  dplyr::select(feature, logFC)
ranks<- cluster.genes$logFC
names(ranks) <- cluster.genes$feature
head(ranks)
VlnPlot(data.rm, features = c("LRRC75A"), pt.size = 0, cols = fine.color)

intersect( tail(names(ranks), 200), head(names(ranks.n), 200) )
VlnPlot(data.rm, features = intersect( tail(names(ranks), 200), head(names(ranks.n), 200) ), pt.size = 0, cols = fine.color)
VlnPlot(data.rm, features =c("CLDN7", "PLCG2"), pt.size = 0, cols = fine.color)






#https://www.nature.com/articles/s41417-021-00303-x
act <- c("CXCL9","CXCL10", "CXCL11", "CXCL12", 
           "CCL4", "CCL5", "CCL20")
sup <- c("CCL17", "CCL22", "CCL2", "CCL5", "CXCL5", "CXCL2", "CXCL6", 
         "TGFB1", "CSF1", "CSF2")
act <- list(act)
sup <- list(sup)
data.score <- AddModuleScore(
  object = data.rm,
  features = act,
  ctrl = 5,
  name = 'Activation'
)
data.score <- AddModuleScore(
  object = data.score,
  features = sup,
  ctrl = 5,
  name = 'Suppression'
)
library(scales)
data.score@meta.data$Activation <- rescale(data.score@meta.data$Activation1)
FeaturePlot(data.score, features = "Activation")
data.score@meta.data$Suppression <- rescale(data.score@meta.data$Suppression1)
Idents(data.score) <- 'orig.ident'
#FeaturePlot(data.score[,grep("Ctrl", Idents(data.score))] , features = "Inflammation_Score1")
#FeaturePlot(data.score[,grep("LPS", Idents(data.score))] , features = "Inflammation_Score1")
#FeaturePlot(data.score[,grep("Combo", Idents(data.score))] , features = "Inflammation_Score1")
VlnPlot(data.score, features = "Activation", cols = fine.color, pt.size = 0)
VlnPlot(data.score, features = "Activation1", cols = fine.color, pt.size = 0)
VlnPlot(data.score, features = "Suppression", cols = fine.color, pt.size = 0)
VlnPlot(data.score, features = "Suppression1", cols = fine.color, pt.size = 0)

act <- c("CXCL9","CXCL10", "CXCL11", "CXCL12", 
         "CCL4", "CCL5", "CCL20")
sup <- c("CCL17", "CCL22", "CCL2", "CCL5", "CXCL5", "CXCL2", "CXCL6", 
         "TGFB1", "CSF1", "CSF2")
FeaturePlot(data.score, features = act, pt.size = 0)
FeaturePlot(data.score, features = sup, pt.size = 0)
FeaturePlot(data.score, features = c("CCL20", "CXCL2","TGFB1", "CSF1"), cols = c("grey", "darkred", "black"))


##Tight junction GSEA
Idents(data.rm) <- 'orig.ident'
tmp <- data.rm[,grep("N_Vehicle|N_LPS", Idents(data.rm))]
pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
head(pbmc.genes)
m_df<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #Hallmark gene sets
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
pbmc.genes %>%
  dplyr::filter(group == "N_Vehicle") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
cluster.genes<- pbmc.genes %>%
  dplyr::filter(group == "N_Vehicle") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks<- cluster.genes$auc
names(ranks) <- cluster.genes$feature
head(ranks)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()
ggplot(fgseaResTidy %>% filter(padj < 0.95) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 30)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways (Adenoma vs Normal)") + 
  theme_minimal()
plotEnrichment(fgsea_sets[["KEGG_MAPK_SIGNALING_PATHWAY"]],
               ranks) + labs(title="KEGG_MAPK_SIGNALING_PATHWAY")
plotEnrichment(fgsea_sets[["KEGG_RNA_DEGRADATION"]],
               ranks) + labs(title="KEGG_RNA_DEGRADATION")
plotEnrichment(fgsea_sets[["KEGG_TIGHT_JUNCTION"]],
               ranks) + labs(title="KEGG_TIGHT_JUNCTION")

Idents(data.rm) <- 'orig.ident'
tmp <- data.rm[,grep("N_Vehicle|N_Combo", Idents(data.rm))]
pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
head(pbmc.genes)
m_df<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #Hallmark gene sets
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
pbmc.genes %>%
  dplyr::filter(group == "N_Vehicle") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
cluster.genes<- pbmc.genes %>%
  dplyr::filter(group == "N_Vehicle") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks<- cluster.genes$auc
names(ranks) <- cluster.genes$feature
head(ranks)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()
ggplot(fgseaResTidy %>% filter(padj < 0.95) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 30)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways (Adenoma vs Normal)") + 
  theme_minimal()
plotEnrichment(fgsea_sets[["KEGG_MAPK_SIGNALING_PATHWAY"]],
               ranks) + labs(title="KEGG_MAPK_SIGNALING_PATHWAY")
plotEnrichment(fgsea_sets[["KEGG_RNA_DEGRADATION"]],
               ranks) + labs(title="KEGG_RNA_DEGRADATION")
plotEnrichment(fgsea_sets[["KEGG_TIGHT_JUNCTION"]],
               ranks) + labs(title="KEGG_TIGHT_JUNCTION")


Idents(data.rm) <- 'orig.ident'
tmp <- data.rm[,grep("A_Ctrl|A_LPS", Idents(data.rm))]
pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
head(pbmc.genes)
m_df<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #Hallmark gene sets
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
pbmc.genes %>%
  dplyr::filter(group == "A_Ctrl") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
cluster.genes<- pbmc.genes %>%
  dplyr::filter(group == "A_Ctrl") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks<- cluster.genes$auc
names(ranks) <- cluster.genes$feature
head(ranks)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()
ggplot(fgseaResTidy %>% filter(padj < 0.95) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 30)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways (Adenoma vs Normal)") + 
  theme_minimal()
plotEnrichment(fgsea_sets[["KEGG_MAPK_SIGNALING_PATHWAY"]],
               ranks) + labs(title="KEGG_MAPK_SIGNALING_PATHWAY")
plotEnrichment(fgsea_sets[["KEGG_RNA_DEGRADATION"]],
               ranks) + labs(title="KEGG_RNA_DEGRADATION")
plotEnrichment(fgsea_sets[["KEGG_TIGHT_JUNCTION"]],
               ranks) + labs(title="KEGG_TIGHT_JUNCTION")



##Normal vs Adenoma
Idents(data.rm) <- 'orig.ident'
tmp <- data.rm[,grep("Ctrl", Idents(data.rm))]
pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
head(pbmc.genes)
m_df<- msigdbr(species = "Homo sapiens", category = "H") #Hallmark gene sets
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
pbmc.genes %>%
  dplyr::filter(group == "N_Vehicle") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
cluster.genes<- pbmc.genes %>%
  dplyr::filter(group == "N_Vehicle") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks<- cluster.genes$auc
names(ranks) <- cluster.genes$feature
head(ranks)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()

ggplot(fgseaResTidy %>% filter(padj < 0.95) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 30)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways (Adenoma vs Normal)") + 
  theme_minimal()
plotEnrichment(fgsea_sets[["HALLMARK_P53_PATHWAY"]],
               ranks) + labs(title="P53")
plotEnrichment(fgsea_sets[["HALLMARK_OXIDATIVE_PHOSPHORYLATION"]],
               ranks) + labs(title="OXIDATIVE PHOSPHORYLATION")
plotEnrichment(fgsea_sets[["HALLMARK_MYC_TARGETS_V1"]],
               ranks) + labs(title="MYC_TARGETS")
plotEnrichment(fgsea_sets[["HALLMARK_G2M_CHECKPOINT"]],
               ranks) + labs(title="G2M_CHECKPOINT")
plotEnrichment(fgsea_sets[["HALLMARK_KRAS_SIGNALING_DN"]],
               ranks) + labs(title="HALLMARK_KRAS_SIGNALING_DN")
plotEnrichment(fgsea_sets[["HALLMARK_KRAS_SIGNALING_UP"]],
               ranks) + labs(title="HALLMARK_KRAS_SIGNALING_UP")

#C2 CP, BIOCARTA
Idents(data.rm) <- 'orig.ident'
tmp <- data.rm[,grep("Ctrl", Idents(data.rm))]
pbmc.genes <- wilcoxauc(tmp, 'orig.ident')
head(pbmc.genes)
m_df<- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #Hallmark gene sets
head(m_df)
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
pbmc.genes %>%
  dplyr::filter(group == "N_Vehicle") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)
cluster.genes<- pbmc.genes %>%
  dplyr::filter(group == "N_Vehicle") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
ranks<- cluster.genes$auc
names(ranks) <- cluster.genes$feature
head(ranks)
fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  head()
grep("ERK", fgseaResTidy$pathway, value = T)
ggplot(fgseaResTidy %>% filter(padj < 0.95) %>% head(n= 50), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 30)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways (Adenoma vs Normal)") + 
  theme_minimal()
plotEnrichment(fgsea_sets[["KEGG_MAPK_SIGNALING_PATHWAY"]],
               ranks) + labs(title="KEGG_MAPK_SIGNALING_PATHWAY")
plotEnrichment(fgsea_sets[["KEGG_COLORECTAL_CANCER"]],
               ranks) + labs(title="KEGG_COLORECTAL_CANCER")
plotEnrichment(fgsea_sets[["KEGG_TIGHT_JUNCTION"]],
               ranks) + labs(title="KEGG_TIGHT_JUNCTION")


VlnPlot(tmp, features = "KRAS", cols = fine.color, pt.size = 0)+ 
  stat_compare_means(method = "t.test", ref.group = "N_Vehicle") 




###9. Cell cycle####
cc.genes
s.genes <-c(  cc.genes$s.genes)
g2m.genes <- c(cc.genes$g2m.genes)
data.combined <- CellCycleScoring(data.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
DimPlot(data.combined, cols = brewer.pal(3,"Dark2"), pt.size = 1)
FeaturePlot(data.combined, features = c("TOP2A","CCNB2","CCND2", "CCND3", "CDC14A", "CDC14B", "CDKN1A"))
FeaturePlot(data.combined, features = c("TOP2A","CCNB2","CDC14B", "CDKN1A", "MCM6", "MKI67"), ncol = 3)

cyclin.genes <- grep("^CCN[ABDE][0-9]$", rownames(data.rm), value = T)
Idents(data.combined) <- "new.cluster.ids"
VlnPlot(data.combined, features = cyclin.genes, pt.size = 0, cols =  pal.rename)
VlnPlot(data.combined, features = c("CCNA2", "CCNB1", "CCND3", "CCNB2", "PCNA", "TOP2A", "MCM6", "MKI67"), pt.size = 0, cols =  pal.rename, ncol = 4)
VlnPlot(data.combined, features = c("CCNA2", "CCNB1", "CCND3", "CCNB2", "PCNA", "TOP2A", "MCM6", "MKI67"), pt.size = 0, cols =  pal.rename, ncol = 4)
Idents(data.combined) <- "orig.ident"
VlnPlot(data.combined, features = c("CCND1", "CCND2"), pt.size = 0, cols =  fine.color)
FeaturePlot(data.combined, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), pt.size = 0, ncol = 2)


data.combined@meta.data$cycling <- "Cycling"
data.combined@meta.data$cycling[grep("G1", data.combined@meta.data$Phase)]  <- "Non-cycling"
Idents(data.combined) <- "cycling"
DimPlot(data.combined, cols = c("gold2","grey"), pt.size = 1) + ggtitle("Cycling Cells")

library(igraph)
library(schex)
library(scater)
library(scran)
library(ggrepel)
tenx_pbmc3k <- make_hexbin(data.rm, nbins = 40, 
                           dimension_reduction = "UMAP", use_dims=c(1,2))
tenx_pbmc3k@assays$RNA <- tenx_pbmc3k@assays$integrated
plot_hexbin_density(tenx_pbmc3k)
plot_hexbin_feature(tenx_pbmc3k, type="data", feature="MKI67", 
                    action="prop_0", xlab="UMAP1", ylab="UMAP2", 
                    title=paste0("Mean of ", "MKI67"))
plot_hexbin_gene(tenx_pbmc3k, type="data", gene="TOP2A", 
                 action="prop_0", xlab="UMAP1", ylab="UMAP2", 
                 title=paste0("Mean of ", "TOP2A"))




###9. Pseudotime ####
library(SeuratWrappers)
library("monocle3")

cds <- as.cell_data_set(data.combined)
cds <- cluster_cells(cds, resolution=2.5e-3)

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p1+p2
#integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
#cds <- as.cell_data_set(integrated.sub)
cds <- learn_graph(cds, use_partition = TRUE, verbose = TRUE)
plot_cells(cds,
           color_cells_by = "cluster",
           label_groups_by_cluster=F,
           label_leaves=F,
           label_branch_points=FALSE)

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 17])) #Define which cluster is the ending point
plot_cells(cds,
           color_cells_by = "pseudotime",cell_size = 1.2,
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=F,
           label_roots = FALSE,
           trajectory_graph_color = "black")
integrated.sub <- as.Seurat(cds, assay = NULL)
FeaturePlot(integrated.sub, features = "monocle3_pseudotime")+ggtitle("ColoTime")



###10. Adenoma vs Normal #####
#https://www.nature.com/articles/s41467-020-17186-5
VlnPlot(data.combined, features = c("APC", "KRAS", "BRAF", "NRAS"), pt.size = 0, cols = fine.color)
FeaturePlot(data.combined, features = c("APC", "KRAS", "BRAF", "NRAS"))
FeaturePlot(data.rm, features = c("APC", "KRAS", "BRAF", "NRAS"))


###11. Cell prioritization  #####
#devtools::install_github("neurorestore/Augur")
#https://github.com/neurorestore/Augur

load(file = "./34all_organoid_seurat_renaming.rds")
Idents(data.combined) <- 'new.cluster.ids'
data.rm <- data.combined[,grep("dying", Idents(data.combined), invert = T)]

library(Augur)
Idents(data.rm) <- "orig.ident"
#Normal 
augur1 = calculate_auc(data.rm[,grep("N_Vehicle|N_LPS", Idents(data.rm))], n_threads = 16, label_col = "orig.ident", cell_type_col = "new.cluster.ids")
augur1$AUC
augur2 = calculate_auc(data.rm[,grep("N_Vehicle|N_Combo", Idents(data.rm))], n_threads = 16, label_col = "orig.ident", cell_type_col = "new.cluster.ids")
augur2$AUC
augur3 = calculate_auc(data.rm[,grep("N_Combo|N_LPS", Idents(data.rm))], n_threads = 16, label_col = "orig.ident", cell_type_col = "new.cluster.ids")
augur3$AUC

tmp <- as.data.frame(augur1$AUC)
rownames(tmp) <- tmp$cell_type
tmp2 <- as.data.frame(augur2$AUC)
rownames(tmp2) <- tmp2$cell_type
augur.df <- merge(tmp, tmp2, by=0, all=TRUE)
rownames(augur.df) <- augur.df$Row.names
augur.df <- augur.df[,c(3,5)]
tmp <- as.data.frame(augur3$AUC)
rownames(tmp) <- tmp$cell_type
augur.df <- merge(augur.df, tmp, by=0, all=TRUE)
rownames(augur.df) <- augur.df$Row.names
augur.df <- augur.df[,c(2,3,5)]
augur.df$cell_type <- rownames(augur.df)
augur.df$cell_type <- factor(augur.df$cell_type, levels = c( "HMGB1+Enterocyte", "LYZ+Enterocyte", "Enterocyte1", "Enterocyte2", "immature Enterocyte", 
                                                             "Goblet1", "immature Goblet","TA1", "TA2", "TA3", "Progenitor","ISC",  
                                                             "EEC", "Tuft" ) )
colnames(augur.df)[1:3] <- c("auc1", "auc2", "auc3")

ggplot(augur.df[complete.cases(augur.df),], aes(x=cell_type, y=auc1)) +
  geom_segment( aes(x=cell_type, xend=cell_type, y=0.5, yend=auc1), color="#E0C1F6") +
  geom_point( color="#4E0087", size=4, alpha=0.6) +
  theme_light() + ylim(0.5,1)+
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                        # positioning the plot title
    legend.title = element_text(face = "bold") 
  )+ 
  labs(x = "Cell Type",                                              # labelling x axis
       y = "AUC",                                        # labeling y axis
       title = "LPS vs Vehicle",        # title
       fill = "Cluster")

ggplot(augur.df[complete.cases(augur.df),], aes(x=cell_type, y=auc2)) +
  geom_segment( aes(x=cell_type, xend=cell_type, y=0.5, yend=auc2), color="#F08C8C") +
  geom_point( color="#CB0404", size=4, alpha=0.6) +
  theme_light() + ylim(0.5,1)+
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                        # positioning the plot title
    legend.title = element_text(face = "bold") 
  )+ 
  labs(x = "Cell Type",                                              # labelling x axis
       y = "AUC",                                        # labeling y axis
       title = "Combo vs Vehicle",        # title
       fill = "Cluster")


ggplot(augur.df[complete.cases(augur.df),], aes(x=cell_type, y=auc3)) +
  geom_segment( aes(x=cell_type, xend=cell_type, y=0.5, yend=auc3), color="#8DF2E7") +
  geom_point( color="#07AD9C", size=4, alpha=0.6) +
  theme_light() + ylim(0.5,1)+
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                        # positioning the plot title
    legend.title = element_text(face = "bold") 
  )+ 
  labs(x = "Cell Type",                                              # labelling x axis
       y = "AUC",                                        # labeling y axis
       title = "Combo vs LPS",        # title
       fill = "Cluster")

#TLR4
tmp <- data.rm[,grep("N_Vehicle", data.rm$orig.ident)]
Idents(tmp) <- "new.cluster.ids"
VlnPlot(tmp, features = c("TLR4", "NRG1"))
FeaturePlot(tmp, features = c("TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "NRG1"), pt.size = 2, cols = brewer.pal(8, "YlOrRd")[c(1,8)], ncol = 4)
FeaturePlot(tmp, features = c("TLR2", "TLR3", "TLR4", "TLR5"), pt.size = 2, cols = brewer.pal(8, "YlOrRd")[c(1,8)], ncol = 2)

 
VlnPlot(data.combined, features = c("TLR1", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7"))



#Adenoma
Idents(data.rm) <- data.rm$orig.ident
augur1 = calculate_auc(data.rm[,grep("A_Vehicle|A_LPS", Idents(data.rm))], n_threads = 16, label_col = "orig.ident", cell_type_col = "new.cluster.ids")
augur1$AUC
augur2 = calculate_auc(data.rm[,grep("A_Vehicle|A_Combo", Idents(data.rm))], n_threads = 16, label_col = "orig.ident", cell_type_col = "new.cluster.ids")
augur2$AUC
augur3 = calculate_auc(data.rm[,grep("A_Combo|A_LPS", Idents(data.rm))], n_threads = 16, label_col = "orig.ident", cell_type_col = "new.cluster.ids")
augur3$AUC

tmp <- as.data.frame(augur1$AUC)
rownames(tmp) <- tmp$cell_type
tmp2 <- as.data.frame(augur2$AUC)
rownames(tmp2) <- tmp2$cell_type
augur.df <- merge(tmp, tmp2, by=0, all=TRUE)
rownames(augur.df) <- augur.df$Row.names
augur.df <- augur.df[,c(3,5)]
tmp <- as.data.frame(augur3$AUC)
rownames(tmp) <- tmp$cell_type
augur.df <- merge(augur.df, tmp, by=0, all=TRUE)
rownames(augur.df) <- augur.df$Row.names
augur.df <- augur.df[,c(2,3,5)]
augur.df$cell_type <- rownames(augur.df)
augur.df$cell_type <- factor(augur.df$cell_type, levels = c( "BEST4+Enterocyte", "LYZ+Enterocyte", "Enterocyte1", "Enterocyte2", 
                                                             "Goblet1", "immature Goblet","TA1",  "Progenitor","ISC",  
                                                             "EEC", "Tuft" ) )

colnames(augur.df)[1:3] <- c("auc1", "auc2", "auc3")
ggplot(augur.df[complete.cases(augur.df),], aes(x=cell_type, y=auc1)) +
  geom_segment( aes(x=cell_type, xend=cell_type, y=0.5, yend=auc1), color="#E0C1F6") +
  geom_point( color="#4E0087", size=4, alpha=0.6) +
  theme_light() + ylim(0.5,1)+
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                        # positioning the plot title
    legend.title = element_text(face = "bold") 
  )+ 
  labs(x = "Cell Type",                                              # labelling x axis
       y = "AUC",                                        # labeling y axis
       title = "LPS vs Vehicle",        # title
       fill = "Cluster")

ggplot(augur.df[complete.cases(augur.df),], aes(x=cell_type, y=auc2)) +
  geom_segment( aes(x=cell_type, xend=cell_type, y=0.5, yend=auc2), color="#F08C8C") +
  geom_point( color="#CB0404", size=4, alpha=0.6) +
  theme_light() + ylim(0.5,1)+
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                        # positioning the plot title
    legend.title = element_text(face = "bold") 
  )+ 
  labs(x = "Cell Type",                                              # labelling x axis
       y = "AUC",                                        # labeling y axis
       title = "Combo vs Vehicle",        # title
       fill = "Cluster")


ggplot(augur.df[complete.cases(augur.df),], aes(x=cell_type, y=auc3)) +
  geom_segment( aes(x=cell_type, xend=cell_type, y=0.5, yend=auc3), color="#8DF2E7") +
  geom_point( color="#07AD9C", size=4, alpha=0.6) +
  theme_light() + ylim(0.5,1)+
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(), 
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                        # positioning the plot title
    legend.title = element_text(face = "bold") 
  )+ 
  labs(x = "Cell Type",                                              # labelling x axis
       y = "AUC",                                        # labeling y axis
       title = "Combo vs LPS",        # title
       fill = "Cluster")




###12. #####
VlnPlot(data.rm, features = "CFTR", pt.size = 0) 
data.combined.avg["CFTR",]





###Grant app , normal only####
load(file = "./34all_organoid_seurat_renaming.rds")
normal <- data.combined[,grep("N_", data.combined$orig.ident)]

DimPlot(normal, cols = pal.rename, pt.size = 2, label = T) + NoLegend()


idx.cluster <- 17
group.df <- data.frame(Percentage =rep(0, 6*idx.cluster), 
                       group=c(rep("N_Vehicle", idx.cluster), rep("N_LPS", idx.cluster), rep("N_Combo", idx.cluster), 
                               rep("A_Vehicle", idx.cluster), rep("A_LPS", idx.cluster), rep("A_Combo", idx.cluster)))
Idents(data.combined) <- data.combined[["orig.ident"]] 
tmp <- data.combined[,grep("N_Vehicle", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[1:idx.cluster] <- table(Idents(tmp))/ncol(tmp)
tmp <- data.combined[,grep("N_LPS", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[(idx.cluster*1+1):(idx.cluster*2)] <- table(Idents(tmp))/ncol(tmp)
tmp <- data.combined[,grep("N_Combo", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[(idx.cluster*2+1):(idx.cluster*3)] <- table(Idents(tmp))/ncol(tmp)
tmp <- data.combined[,grep("A_Vehicle", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[c((idx.cluster*3+1):59, 61:(idx.cluster*4))] <- table(Idents(tmp))/ncol(tmp)
tmp <- data.combined[,grep("A_LPS", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[(idx.cluster*4+1):(idx.cluster*5)] <- table(Idents(tmp))/ncol(tmp)
tmp <- data.combined[,grep("A_Combo", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
table(Idents(tmp))/ncol(tmp)
group.df$Percentage[c((idx.cluster*5+1):93,95:(idx.cluster*6))] <- table(Idents(tmp))/ncol(tmp)

tmp <- data.combined[,grep("N_Vehicle", Idents(data.combined))]
Idents(tmp) <- tmp[["new.cluster.ids"]] 
group.df$Cluster <- rep(c(levels(Idents(tmp))), 6)
group.df$Cluster <- factor(group.df$Cluster, levels = c("dying","HMGB1+Enterocyte", "LYZ+Enterocyte", "BEST4+Enterocyte", "Enterocyte1", "Enterocyte2", "immature Enterocyte", 
                                                        "Goblet1", "Goblet2", "immature Goblet","TA1", "TA2", "TA3", "Progenitor","ISC",  
                                                        "EEC", "Tuft" ))
group.df$Percentage <- round(group.df$Percentage, 4)
group.df$group <- factor(group.df$group, levels = c("N_Vehicle", "N_LPS", "N_Combo","A_Vehicle", "A_LPS", "A_Combo"))

#https://r-charts.com/colors/
group.df[grep("N_", group.df$group),] %>%
  ggplot(data = ., mapping = aes(x = group, y = Percentage, fill = Cluster)) +
  geom_col() +
  geom_text(mapping = aes(label = paste0(Percentage*100, "%") ),              # converting the values to percent
            size = 5,                                             # size of the font
            position = position_stack(vjust = 0.5)) +             # positioning in the middle
  #scale_fill_brewer(palette = "Accent") +                           # coloring the plot
  scale_fill_manual(values = pal.rename)+
  #facet_grid(.~sex) +
  labs(x = "Group",                                              # labelling x axis
       y = "Percentage",                                        # labeling y axis
       title = "Percentage (All)",        # title
       fill = "Cluster") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 90,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  ) #+ 
