library(Seurat) 
library(harmony) 
Cells.sub<- subset(sce@meta.data,celltype=="NK/T cell")
scRNAsub <- subset(sce,cells=row.names(Cells.sub))
scRNAsub <- FindVariableFeatures(scRNAsub,selection.method = "vst",nfeatures=2000) 
scRNAsub <- ScaleData(scRNAsub) 
scRNAsub <- RunPCA(scRNAsub,features = VariableFeatures(scRNAsub)) 
scRNA_harmony <- RunHarmony(scRNAsub,group.by.var="orig.ident")
pc.num=30 
reso.num=0.4 
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction ="harmony", dims = 1:pc.num) %>% FindClusters(resolution = reso.num)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony",dims = 1:pc.num) 
plot1 = DimPlot(scRNA_harmony,reduction = "tsne",group.by="seurat_clusters",label = F,repel = T, shuffle = T ) +
  ggtitle("Cluster")
plot2 = DimPlot(scRNA_harmony,reduction = "tsne",group.by="Group",label = F,repel = T, shuffle = T ) +
  ggtitle("Group")
plot3 = DimPlot(scRNA_harmony,reduction = "tsne",group.by="sample.ident",label = F) +
  ggtitle("Sample")
markers <- FindAllMarkers(object =scRNA_harmony,test.use="wilcox",only.pos =TRUE,logfc.threshold =0.25,min.pct = 0.25) 
celltype=data.frame(ClusterID=0:10,celltype='unkown')
celltype[celltype$ClusterID %in% c(1,2,5,8),2]='CD4+T' 
celltype[celltype$ClusterID %in% c(0,4,6,7,10),2]='CD8+T' 
celltype[celltype$ClusterID %in% c(3,9),2]='NK'
sce.in=scRNA_harmony
sce.in@meta.data$subcelltype0 = NA 
for(i in 1:nrow(celltype)){
  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),"subcelltype0"] <- celltype$celltype[i]} 
sce.in@meta.data$subcelltype0 <- factor(sce.in@meta.data$subcelltype0,levels=c('CD4+T','CD8+T','NK'))
celltype=data.frame(ClusterID=0:10,celltype='unkown')
celltype[celltype$ClusterID %in% c(5),2]='CD4+T_reg' 
celltype[celltype$ClusterID %in% c(1,2,8),2]='CD4+T_naive' 
celltype[celltype$ClusterID %in% c(0,4),2]='CD8+T_cyto'
celltype[celltype$ClusterID %in% c(7),2]='CD8+T_pro' 
celltype[celltype$ClusterID %in% c(6),2]='CD8+T_ex'
celltype[celltype$ClusterID %in% c(10),2]='CD8+T_cm' 
celltype[celltype$ClusterID %in% c(3),2]='NK_CD16-' 
celltype[celltype$ClusterID %in% c(9),2]='NK_CD16+' 
sce.in@meta.data$subcelltype = NA 
for(i in 1:nrow(celltype)){
  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),"subcelltype"] <- celltype$celltype[i]} 
sce.in@meta.data$subcelltype <- factor(sce.in@meta.data$subcelltype,
                                        levels=c('CD4+T_naive','CD4+T_reg',"CD8+T_cm",'CD8+T_cyto',"CD8+T_pro","CD8+T_ex",'NK_CD16-','NK_CD16+'))
genes_to_check = c("CD2","CD3D","CD3E",
                   "CD4", 
                   "CD8A","CD8B",
                   "KLRD1","FCGR3A","KLRC1",
                   "FOXP3","IL2RA",
                   "CCR7","IL7R", "KLF2",
                   "NKG7", "GNLY", "GZMH","GZMB",
                   "MKI67", "TOP2A", "CENPF","STMN1",
                   "LAYN", "LAG3","PDCD1","CTLA4")
pMarker <- DotPlot(sce.in,features = genes_to_check, assay="RNA",group.by = "subcelltype")
sce=sce.in
p.tsne1=DimPlot(sce,reduction = "tsne",group.by = "seurat_clusters",repel = T, shuffle = T,label =F,label.size = 4, pt.size=1)+
  ggtitle("Cluster")
p.tsne3=DimPlot(sce,reduction = "tsne",group.by = "subcelltype",repel = T, shuffle = T,label =F,label.size = 4, pt.size=1)+
  ggtitle("Celltype")+
  scale_color_npg()
p.tsne4=DimPlot(sce,reduction = "tsne",group.by = "Group",repel = T, shuffle = T,label =F,label.size = 4, pt.size=1)+
  ggtitle("Group")+
  scale_color_npg()
p.feature1 <- FeaturePlot(sce, features = c("CD4"), cols = c("grey", "red"),pt.size = 0.5,min.cutoff = "q10",reduction = "tsne",ncol = 1)
p.feature2 <- FeaturePlot(sce, features = c("CD8A"), cols = c("grey", "red"),pt.size = 0.5,min.cutoff = "q10",reduction = "tsne",ncol = 1)
p.feature3 <- FeaturePlot(sce, features = c("KLRC1"), cols = c("grey", "red"),pt.size = 0.5,min.cutoff = "q10",reduction = "tsne",ncol = 1)
scRNA_harmony <- sce
library(ggalluvial)
Ratio <- scRNA_harmony@meta.data %>%
  group_by(Group, subcelltype) %>% 
  summarise(n=n()) %>% 
  mutate(relative_freq = n/sum(n))
ggplot(Ratio, aes(x =Group, y= relative_freq, fill = subcelltype, stratum=subcelltype, alluvium=subcelltype)) 
Cells.sub<- subset(sce@meta.data,subcelltype0=="CD8+T") #CD4+T,CD8+T,NK
scRNA_harmony <- subset(sce,cells=row.names(Cells.sub))
scRNA_harmony@meta.data$Group <- factor(scRNA_harmony@meta.data$Group,levels = c("EOCC","LOCC"))
Idents(scRNA_harmony) <- "Group"
deg_all=FindMarkers(scRNA_harmony, ident.1 = "EOCC", ident.2 = "LOCC", group.by="Group", min.pct = 0.05)
degs <- rownames(filter(deg_all,p_val_adj<0.05)) 
checkpioints <- read.table("genelist_immunecheckpoints.txt")
checkpioints <- checkpioints[,1]
checkpioints <- intersect(checkpioints,degs)
chemokines <- read.table("genelist_chemokines.txt")
chemokines <- chemokinesANDreceptors[,1]
chemokines <- intersect(chemokinesANDreceptors,degs)
hla_genes <- grep("^HLA", rownames(deg_all), value = TRUE)
hla_genes
hla_genes <- intersect(hla_genes,degs)
ils <- read.table("genelist_interleukins.txt")
ils <- ils[,1]
ils <- intersect(ils,degs)
p.checkpioints <- VlnPlot(scRNA_harmony, features = checkpioints, pt.size = 0) + NoLegend()
p.chemokines <- VlnPlot(scRNA_harmony, features = chemokines, pt.size = 0) + NoLegend()
p.HLA <- VlnPlot(scRNA_harmony, features = hla_genes, pt.size = 0) + NoLegend()
p.ils <- VlnPlot(scRNA_harmony, features = ils, pt.size = 0) + NoLegend()
gene_sets_all <- read.csv("Genesets.cd8t.csv")
for (i in c(1:4)){
  gene_sets <- gene_sets_all[,i,drop=F]
  gene_sets <- gene_sets[gene_sets[,1] %in% rownames(scRNA_harmony),,drop=F]
  genes <- as.data.frame(gene_sets)
  genes <- as.list(gene_sets)
  scRNA_harmony <- AddModuleScore(scRNA_harmony,features = genes,name = colnames(gene_sets_all)[i])
  colnames(scRNA_harmony@meta.data)
  my_comparisons <- list(c('EOCC',"LOCC"))
  library(ggpubr)
  p.AddModuleScore <- ggviolin(scRNA_harmony@meta.data, x = "Group", y = paste0(colnames(gene_sets_all)[i],"1"),fill = 'Group')+
    geom_boxplot(width = 0.15, position = position_dodge(0.9), outlier.shape = NA, color = "black")+
    stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
    theme_test()+
    NoLegend() + 
    labs(x = '',y="Scores")+
    ggtitle(colnames(gene_sets_all)[i])
  p.AddModuleScore
  plotname <- paste("plot_", i, sep = "")
  assign(plotname, p.AddModuleScore)
}
Idents(scRNAsub) <- "subcelltype"
markers <- FindAllMarkers(object = scRNAsub,test.use="wilcox",
                          only.pos = F, 
                          logfc.threshold =0.25) 
all.markers =markers %>% dplyr::select(gene,everything()) %>% subset(p_val<0.05) 
mycolor <-  pal_npg()(9)
jjVolcano(diffData = all.markers,
          log2FC.cutoff = 0.25,
          size =3.5,
          fontface = 'italic',
          tile.col = mycolor,
          topGeneN =5)
expr <- AverageExpression(scRNAsub, assays = "RNA", layer = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  
expr <- as.matrix(expr)
library(msigdbr)
human_KEGG = msigdbr(species = "Homo sapiens",
                     category = "H") 
human_KEGG_Set = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
library(GSVA)
params <- gsvaParam(expr = expr, 
                    geneSets = human_KEGG_Set, 
                    kcdf="Gaussian")
gsva.kegg <- gsva(params)
library(pheatmap)
gsva.kegg1 <- as.data.frame(t(scale(t(gsva.kegg))))
pheatmap(gsva.kegg1, 
             show_colnames = T, 
             angle_col = "45",
             fontsize_row = 14, 
             fontsize_col = 18,
             cluster_row = F,
             cluster_col = F,
             color = colorRampPalette(c("#4B5B95", "white", "#CD322C"))(50))