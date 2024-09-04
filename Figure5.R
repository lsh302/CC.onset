library(Seurat) 
library(harmony) 
Cells.sub<- subset(sce@meta.data,celltype=="Fibroblast")
scRNAsub <- subset(sce,cells=row.names(Cells.sub))
scRNAsub <- FindVariableFeatures(scRNAsub,selection.method = "vst",nfeatures=2000) 
scRNAsub <- ScaleData(scRNAsub) 
scRNAsub <- RunPCA(scRNAsub,features = VariableFeatures(scRNAsub)) 
scRNA_harmony <- RunHarmony(scRNAsub,group.by.var="orig.ident")
pc.num=10 
reso.num=0.1 
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction ="harmony", dims = 1:pc.num) %>% FindClusters(resolution = reso.num)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony",dims = 1:pc.num) 
markers <- FindAllMarkers(object =scRNA_harmony,test.use="wilcox",
                          only.pos =TRUE,
                          logfc.threshold =0.25,
                          min.pct = 0.25) 
celltype=data.frame(ClusterID=0:5,celltype='unkown')
celltype[celltype$ClusterID %in% c(0,2),2]='iCAF' 
celltype[celltype$ClusterID %in% c(1),2]='myCAF' 
celltype[celltype$ClusterID %in% c(3),2]='vCAF' 
celltype[celltype$ClusterID %in% c(4,5),2]='apCAF' 
sce.in=scRNA_harmony
sce.in@meta.data$subcelltype = NA 
for(i in 1:nrow(celltype)){
  sce.in@meta.data[which(sce.in@meta.data$seurat_clusters == celltype$ClusterID[i]),"subcelltype"] <- celltype$celltype[i]} 
sce.in@meta.data$subcelltype <- factor(sce.in@meta.data$subcelltype,levels=c('apCAF','iCAF','myCAF','vCAF'))
genes_to_check =  c("DCN","COL1A1","COL3A1",
                    "IL6", "CXCL2","CCL2",
                    "INHBA","MMP11",
                    "MCAM","RGS5","MYH11",
                     "HLA-DRB1","HLA-DRA","HLA-DPA1","HLA-DPB1")
pMarker <- DotPlot(sce.in,features = genes_to_check, assay="RNA",group.by = "subcelltype") 
sce=sce.in
p.tsne1=DimPlot(sce,reduction = "tsne",group.by = "seurat_clusters",repel = T, shuffle = T,label =F,label.size = 4, pt.size=1.5)+
  ggtitle("Cluster")
p.tsne2=DimPlot(sce,reduction = "tsne",group.by = "subcelltype",repel = T, shuffle = T,label =F,label.size = 4, pt.size=1.5)+
  ggtitle("Celltype")+
  scale_color_npg()
p.tsne3=DimPlot(sce,reduction = "tsne",group.by = "Group",repel = T, shuffle = T,label =F,label.size = 4, pt.size=1.5)+
  ggtitle("Group")+
  scale_color_npg()
p.feature <- FeaturePlot(sce, features = c("HLA-DRA","CXCL2","INHBA","MCAM"), cols = c("grey", "red"),
                         reduction = "tsne",pt.size=1.2,
                         min.cutoff = "q10",ncol = 2)
scRNA_harmony <- sce
library(ggalluvial)
Ratio <- scRNA_harmony@meta.data %>%
  group_by(Group, subcelltype) %>% 
  summarise(n=n()) %>% 
  mutate(relative_freq = n/sum(n))
convert_to_percentage <- function(decimal) {
  percentage <- decimal * 100
  formatted_percentage <- sprintf("%.2f%%", percentage)
  return(formatted_percentage)
}
Ratio$percent <- convert_to_percentage(Ratio$relative_freq)
Ratio$labs <- paste0(Ratio$subcelltype, " (", Ratio$percent, ")")
library(ggpubr)
Ratio1 <- filter(Ratio,Group == "EOCC")
Ratio2 <- filter(Ratio,Group == "LOCC")
p1 <- ggpie(Ratio1,"n",
      label = "labs", 
      lab.pos = "in",  
      fill = "subcelltype",                            
      palette = color)+
  ggtitle("EOCC")+theme(legend.position = "right",plot.title = element_text(hjust = 0.5))
p2 <- ggpie(Ratio2,"n",
            label = "labs",  
            lab.pos = "in",  
            fill = "subcelltype",                            
            palette = color)+
  ggtitle("LOCC")+theme(legend.position = "right",plot.title = element_text(hjust = 0.5))
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
library(monocle)
scRNAsub <- sc
exp <- scRNAsub@assays[["RNA"]]@counts
fdata <- data.frame(gene_short_name = row.names(scRNAsub), row.names = row.names(scRNAsub))
pdata <- scRNAsub@meta.data
fd <- new("AnnotatedDataFrame", data = fdata) 
pd <- new("AnnotatedDataFrame", data = pdata)
CDS <- newCellDataSet(cellData = exp,phenoData = pd,featureData = fd)
CDS <- estimateSizeFactors(CDS)
CDS <- estimateDispersions(CDS)
CDS <- detectGenes(CDS, min_expr = 0.1)
expressed_genes <-  row.names(subset(fData(CDS),num_cells_expressed >= 20));length(expressed_genes)
clustering_DEG_genes <-differentialGeneTest(CDS[expressed_genes,],fullModelFormulaStr = '~subcelltype')
ordering_genes <-row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:2000]
CDS <- setOrderingFilter(CDS,ordering_genes = ordering_genes)
plot_ordering_genes(CDS)
CDS <- reduceDimension(CDS, method = 'DDRTree')
CDS <-orderCells(CDS)
plot_cell_trajectory(CDS, color_by = "Pseudotime")
plot_cell_trajectory(CDS, color_by = "State")
plot_cell_trajectory(CDS, color_by = "subcelltype")
plot_cell_trajectory(CDS, color_by = "Group")
diff_test_res=differentialGeneTest(CDS[expressed_genes[1:1000],],fullModelFormulaStr = "~sm.ns(Pseudotime)",cores=4)
genes <- c("CCL2","CXCL1","CXCL2","CXCL3","CXCL12")
plot_genes_in_pseudotime(CDS[genes,], color_by="Group", ncol = 1) 
library(CytoTRACE)
monocle_meta <- data.frame(t(CDS@reducedDimS), 
                           CDS$Pseudotime, 
                           CDS$State, 
                           CDS$Group)
colnames(monocle_meta) <- c("C1", "C2", "Pseudotime", "State", "group")
phenot1 <- monocle_meta$group
phenot1 <- as.character(phenot1)
names(phenot1) <- rownames(monocle_meta)
emb_monocle <- monocle_meta[,1:2]
mat<-as.matrix(scRNAsub@assays$RNA@counts)
results <- CytoTRACE(mat = mat)
plotCytoTRACE(results, phenotype = phenot1, emb = emb_monocle)
scores <- as.data.frame(results[["CytoTRACE"]])
colnames(scores) <- c("CytoTRACE_scores")
all(names(phenot1)==rownames(scores))
scores$group <- phenot1
scores$group <- factor(scores$group,levels=c("EOCC","LOCC"))
my_comparisons <- list(c('EOCC',"LOCC"))
ggboxplot(scores, x = "group", y = "CytoTRACE_scores",fill = 'group',outlier.shape = NA)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") 