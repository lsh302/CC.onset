library(Seurat) 
library(harmony) 
dir_name <- list.files("GSE197461/") 
scRNAlist <- list()
for(i in 1:length(dir_name)){
  counts <- Read10X(data.dir = paste("GSE197461/",dir_name[i],sep =""))
  scRNAlist[[i]] <- CreateSeuratObject(counts,project = dir_name[i],min.cells = 3, min.features=200)
}
dir_name1 <- list.files("BSST1035/") 
scRNAlist1 <- list()
for(i in 1:length(dir_name1)){
  counts <- read.csv(paste("BSST1035/",dir_name1[i],sep =""),row.names=1)
  scRNAlist1[[i]] <- CreateSeuratObject(counts,project =sub("\\.csv","",dir_name1[i]),min.cells = 3, min.features=200)
  scRNAlist1[[i]]@meta.data$orig.ident <- paste0("BSST1035_",scRNAlist1[[i]]@meta.data$orig.ident)
  Idents(scRNAlist1[[i]]) <- scRNAlist1[[i]]@meta.data$orig.ident
}
scRNAlist <- c(scRNAlist,scRNAlist1)
for(i in 1:length(scRNAlist)){
  sc <- scRNAlist[[i]] 
  sc[["mt_percent"]] <- PercentageFeatureSet(sc,pattern ="^MT-") 
  HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ") 
  HB_m <- match(HB_genes,rownames(sc@assays$RNA)) 
  HB_genes <- rownames(sc@assays$RNA)[HB_m]
  HB_genes <- HB_genes[!is.na(HB_genes)] 
  sc[["HB_percent"]] <- PercentageFeatureSet(sc,features=HB_genes) 
  scRNAlist[[i]] <- sc 
  rm(sc)
}
scRNAlist <- lapply(X = scRNAlist,FUN = function(x){ 
  x<- subset(x,
            subset = nFeature_RNA > 200 & nFeature_RNA <5000 &
            mt_percent < 20 & 
            HB_percent < 3 &
            nCount_RNA < quantile(nCount_RNA,0.97) & 
            nCount_RNA > 1000)})

scRNAlist <- merge(x=scRNAlist[[1]],y=scRNAlist[-1]) 
scRNAlist <- NormalizeData(scRNAlist) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures =2000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 30,verbose = T) 
scRNA_harmony <- RunHarmony(scRNAlist,group.by.vars = "orig.ident") 
pc.num=25 
reso.num=0.2
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction ="harmony", dims = 1:pc.num) %>% FindClusters(resolution = reso.num)
scRNA_harmony <- RunTSNE(scRNA_harmony,reduction = "harmony",dims = 1:pc.num) 
markers <- FindAllMarkers(object =scRNA_harmony,test.use="wilcox",
                          only.pos =TRUE,
                          logfc.threshold =0.25,
                          min.pct = 0.25) 
group_info <- read.csv("Group_info.csv",check.names = F)
for(i in 1:nrow(group_info)){
  scRNA_harmony@meta.data[which(scRNA_harmony@meta.data$orig.ident == group_info$Sample_ID[i]),"Group"] <- group_info$Group_ID[i]} 
celltype=data.frame(ClusterID=0:16,celltype='unkown')
celltype[celltype$ClusterID %in% c(0,1,3,6,11,14,16),2]='Epithelial cell' 
celltype[celltype$ClusterID %in% c(12),2]='Endothelial cell' 
celltype[celltype$ClusterID %in% c(5,9),2]='Fibroblast' 
celltype[celltype$ClusterID %in% c(2,4,15),2]='NK/T cell' 
celltype[celltype$ClusterID %in% c(10),2]='B cell' 
celltype[celltype$ClusterID %in% c(8),2]='Plasma cell' 
celltype[celltype$ClusterID %in% c(7),2]='Myeloid cell' 
celltype[celltype$ClusterID %in% c(13),2]='Mast cell' 
sce.in=scRNA_harmony
sce.in@meta.data$celltype = NA 
for(i in 1:nrow(celltype)){
  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),"celltype"] <- celltype$celltype[i]} 
genes_to_check = c("KRT18","EPCAM","CDKN2A",#Epithelial
                   "PECAM1","VWF","ENG",#Endothelial 
                   "DCN","COL1A1","COL3A1",#Fibroblast
                   "CD3D","CD3E","NKG7",#NK/T
                   "CD14","CD163","CSF1R",#Myeloid
                   "MS4A1","CD19","CD79A",#	B cell
                   "MZB1","IGHG1","TNFRSF17",#Plasma
                   "TPSAB1","MS4A2","CPA3" #Mast
)
sce.in@meta.data$celltype <- factor(sce.in@meta.data$celltype,levels = c("Epithelial cell","Endothelial cell","Fibroblast",
                                                                         "NK/T cell","Myeloid cell","B cell","Plasma cell","Mast cell"))
pMarker <- DotPlot(sce.in,features = genes_to_check, assay="RNA",group.by = "celltype")
sce=sce.in
plot1 = DimPlot(sce,reduction = "tsne",group.by="seurat_clusters",label = F,repel = T, shuffle = T ) +
  ggtitle("Cluster")
plot2 = DimPlot(sce,reduction = "tsne",group.by="celltype",label = T, shuffle = T) +
  scale_color_npg()+NoLegend()+
  ggtitle("Celltype")
plot3 = DimPlot(sce,reduction = "tsne",group.by="Group",label = F,repel = T, shuffle = T ) +
  scale_color_npg()+
  ggtitle("Group")
p.features <- FeaturePlot(sce, features = c("CDKN2A", "PECAM1", "COL3A1", "CD3D","CD14","MS4A1","MZB1","TPSAB1"), cols = c("grey", "red"),
                         reduction="tsne",min.cutoff = "q10",ncol = 2)

scRNA_harmony <- sce
stat<-table(scRNA_harmony$sample.ident,scRNA_harmony$celltype) 
stat<-as.data.frame(stat)
colnames(stat)<-c('sample','celltype','Freq') 
p1<-ggplot(data = stat,aes(x = sample,y = Freq, fill = celltype))+
  geom_bar(stat = 'identity',position = 'fill')+
  labs(x = "",y = "Ratio")+
  mytheme+
  scale_fill_npg()
library(ggalluvial)
Ratio <- scRNA_harmony@meta.data %>%
  group_by(Group, celltype) %>% 
  summarise(n=n()) %>% 
  mutate(relative_freq = n/sum(n))
p <- ggplot(Ratio, aes(x =Group, y= relative_freq, fill = celltype,  
                       stratum=celltype, alluvium=celltype)) + 
  geom_col(width = 0.5)+
  labs(x='',y = 'Ratio')+
  mytheme+
  scale_fill_npg()
scRNA_harmony@meta.data$Group <- factor(scRNA_harmony@meta.data$Group,levels = c("EOCC","LOCC"))
Idents(scRNA_harmony) <- "Group"
table(scRNA_harmony@meta.data$Group,scRNA_harmony@meta.data$celltype)
deg_all=FindMarkers(scRNA_harmony, ident.1 = "EOCC", ident.2 = "LOCC", group.by="Group", min.pct = 0.05)
volcano.data <- deg_all 
volcano.data <- volcano.data[,c(2,5)]
colnames(volcano.data) <- c("log2FoldChange","padj")
volcano.data$change <- as.factor(ifelse(volcano.data$padj<0.05 & abs(volcano.data$log2FoldChange)>log2(2),  
                                        ifelse(volcano.data$log2FoldChange>log2(2), "Up", "Down"), "NoDiff"))
volcano.data$name <- rownames(volcano.data)
up <- volcano.data[volcano.data$log2FoldChange>1 & volcano.data$padj<0.05,]
down <- volcano.data[volcano.data$log2FoldChange< -1 & volcano.data$padj<0.05,]
up.top <- rownames(up[order(up$padj,-up$log2FoldChange),])[1:6]
down.top <- rownames(down[order(down$padj,down$log2FoldChange),])[1:6]
top_genes <- c(up.top,down.top)
top_genes_data <- volcano.data[top_genes,]
library(ggrepel)
valcano <- ggplot(data=volcano.data, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  ggtitle("EOCC vs LOCC") + 
  scale_color_manual(name="", values=c("#ae3137", "#497aa2", "grey"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-log2(2), log2(2)), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.05), lty=2, col="gray", lwd=0.5)
library(clusterProfiler)
degs <- rownames(filter(deg_all,p_val_adj<0.05))
library(org.Hs.eg.db)
entrezid_all = mapIds(x = org.Hs.eg.db,  
                      keys = degs, 
                      keytype = "SYMBOL", 
                      column = "ENTREZID") 
entrezid_all = na.omit(entrezid_all) 
entrezid_all = data.frame(entrezid_all) 
KEGG_enrich = enrichKEGG(gene = entrezid_all[,1], 
                         keyType = "kegg",
                         organism= "human",  
                         pAdjustMethod = "fdr", 
                         pvalueCutoff = 1, 
                         qvalueCutoff =1) 
KEGG_enrich  = data.frame(KEGG_enrich)
split_columns1 <- strsplit(as.character(KEGG_enrich$GeneRatio), "/")
df_split1 <- data.frame(do.call(rbind, split_columns1))
KEGG_enrich$GeneRatio <- as.numeric(df_split1[,1])/as.numeric(df_split1[,2])
KEGG_enrich <- filter(KEGG_enrich,p.adjust < 0.05) 
kegg_enrich$Description <- factor(kegg_enrich$Description,levels=kegg_enrich$Description)
p.kegg <- ggplot(kegg_enrich,aes(y=Description,x=GeneRatio))+
  geom_point(aes(size=Count,color=p.adjust))+ 
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="GeneRatio",y="",title="EOCC vs LOCC")+
  theme_bw()