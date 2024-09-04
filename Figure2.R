library(Seurat) 
library(harmony)
Cells.sub<- subset(sce@meta.data,celltype=="Epithelial cell")
scRNAsub <- subset(sce,cells=row.names(Cells.sub))
scRNAsub <- FindVariableFeatures(scRNAsub,selection.method = "vst",nfeatures=2000) 
scRNAsub <- ScaleData(scRNAsub) 
scRNAsub <- RunPCA(scRNAsub,features = VariableFeatures(scRNAsub)) 
scRNA_harmony <- RunHarmony(scRNAsub,group.by.var="orig.ident")
pc.num=15 
reso.num=0.2 
scRNA_harmony <- FindNeighbors(scRNA_harmony,reduction ="harmony", dims = 1:pc.num) %>% FindClusters(resolution = reso.num)
scRNA_harmony <- RunTSNE(scRNA_harmony, reduction = "harmony",dims = 1:pc.num) 
plot1 = DimPlot(scRNA_harmony,reduction = "tsne",group.by="seurat_clusters",label = F,repel = T, shuffle = T ) +
  ggtitle("Cluster")
plot2 = DimPlot(scRNA_harmony,reduction = "tsne",group.by="Group",label = F,repel = T, shuffle = T ) +
  ggtitle("Group")
plot3 = DimPlot(scRNA_harmony,reduction = "tsne",group.by="sample.ident",label = F) +
  ggtitle("Sample")
stat<-table(scRNA_harmony$sample.ident,scRNA_harmony$seurat_clusters)
stat<-as.data.frame(stat)
colnames(stat)<-c('sample','celltype','Freq')
p1<-ggplot(data = stat,aes(x = sample,y = Freq, fill = celltype))+
  geom_bar(stat = 'identity',position = 'fill')+
  labs(x = "",y = "Ratio")+
  scale_fill_npg()
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
p.checkpioints <- VlnPlot(scRNA_harmony, features = checkpioints, pt.size = 0) 
p.chemokines <- VlnPlot(scRNA_harmony, features = chemokines, pt.size = 0) 
p.HLA <- VlnPlot(scRNA_harmony, features = hla_genes, pt.size = 0) 
seurat_object <- scRNA_harmony
library(clusterProfiler)
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
p.kegg <- ggplot(kegg_enrich,
             aes(y=Description,x=GeneRatio))+
  geom_point(aes(size=Count,color=p.adjust))+ 
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(p.adjust,size="Count"), 
       x="GeneRatio",y="",title="EOCC vs LOCC")+
  theme_bw()
check.list <- c(checkpioints,chemokines,hla_genes)
exp <- read.csv("TCGA_TPM_CESC_expdataTumor_304.csv", header=T, check.names=F,row.names = 1)
exp_marker <- exp[check.list,]
exp_marker <- as.data.frame(t(exp_marker))
dat<-read.csv("TCGA_clinicalData_CESC.csv",header = T,row.names = 1,check.names = F)
sameSample=intersect(row.names(dat), row.names(exp_marker))
dat=dat[sameSample,,drop=F]
epi=exp_marker[sameSample,,drop=F]
all(rownames(dat)==rownames(exp_marker)) 
rt=cbind(exp_marker, dat)
library(survival)
library(survminer)
range <- c(1:length(colnames(exp_marker))) 
for (j in range){
  i <- colnames(rt)[j] 
  rt1 <- rt[,c(j,(lengthe(colnames(exp_marker))+1):(lengthe(colnames(exp_marker))+8))] 
  colnames(rt1) <- c("cluster","OS","OS.time","DSS","DSS.time","DFI","DFI.time","PFI","PFI.time")
  rt1$group = ifelse(rt1$cluster >=median(rt1$cluster),'high','low')
  #OS
  survdata = Surv(time = rt1$OS.time,           
                  event = rt1$OS == 1) 
  KMfit <- survfit(survdata ~ rt1$group)  
  p1 <- ggsurvplot(KMfit,                       
                   data = rt1,               
                   pval = TRUE,                 
                   surv.median.line = "hv",     
                   conf.int = FALSE,
                   risk.table = FALSE,           
                   xlab = "Follow up time(months)",  
                   ylab = "OS",
                   break.x.by = 12,             
                   xlim = c(0,120),
                   title=i,
  )
  #DSS
  survdata = Surv(time = rt1$DSS.time,            
                  event = rt1$DSS == 1) 
  KMfit <- survfit(survdata ~ rt1$group)   
  p2 <- ggsurvplot(KMfit,                       
                   data = rt1,             
                   pval = TRUE,               
                   surv.median.line = "hv",   
                   conf.int = FALSE,
                   risk.table = FALSE,         
                   xlab = "Follow up time(months)",  
                   ylab = "CSS",
                   break.x.by = 12,             
                   xlim = c(0,120),
                   title=i,
  )
  #PFI
  survdata = Surv(time = rt1$PFI.time,            
                  event = rt1$PFI == 1) 
  KMfit <- survfit(survdata ~ rt1$group) 
  p3 <- ggsurvplot(KMfit,                      
                   data = rt1,               
                   pval = TRUE,                 
                   surv.median.line = "hv",    
                   conf.int = FALSE,
                   risk.table = FALSE,          
                   xlab = "Follow up time(months)", 
                   ylab = "PFI",
                   break.x.by = 12,             
                   xlim = c(0,120),
                   title=i,
  )
  plot_grid(p1$plot+theme1,p2$plot+theme1,p3$plot+theme1) 
  ggsave(filename=paste0(i,"_survival_ROCcutoff.pdf"),height = 8,width = 12)
}
scRNAsub <- scRNA_harmony
Idents(scRNAsub) <- "seurat_clusters"
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
          topGeneN =3)
expr <- AverageExpression(scRNAsub, assays = "RNA", layer = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  
expr <- as.matrix(expr)
library(msigdbr)
human_KEGG = msigdbr(species = "Homo sapiens",category = "H") 
human_KEGG_Set = human_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
library(GSVA)
params <- gsvaParam(expr = expr, 
                    geneSets = human_KEGG_Set, 
                    kcdf="Gaussian")
gsva.kegg <- gsva(params)
library(pheatmap)
gsva.kegg1 <- as.data.frame(t(scale(t(gsva.kegg))))
p<- pheatmap(gsva.kegg1, 
             show_colnames = T, 
             angle_col = "45",
             fontsize_row = 12, 
             fontsize_col = 14,
             cluster_row = F,
             cluster_col = F,
             color = colorRampPalette(c("#4B5B95", "white", "#CD322C"))(50))
filelist<-paste0("#CancerSEA/",dir("#CancerSEA/")) 
CancerSEA_Set<-list()
for(i in 1:length(filelist)){   
  CancerSEA_Set[[i]]<-read.table(filelist[i],header = T,sep = "\t")[,2]
}
names(CancerSEA_Set) <-  str_extract(filelist, "(?<=/).*?(?=\\.)") 
params <- gsvaParam(expr = expr, 
                    geneSets = CancerSEA_Set, 
                    kcdf="Gaussian")
gsva.CancerSEA <- gsva(params)
gsva.CancerSEA1 <- as.data.frame(t(scale(t(gsva.CancerSEA))))
p<- pheatmap(gsva.CancerSEA1, 
             show_colnames = T, 
             angle_col = "45", 
             fontsize_row = 8.5, 
             fontsize_col = 8.5,
             cluster_row = F,
             cluster_col = F,
             color = colorRampPalette(c("#4B5B95", "white", "#CD322C"))(50))
library(infercnv)
Idents(sce) <- "celltype"
table(sce@meta.data$celltype)
seurat_object <- subset(sce,idents = c("Epithelial cell",'Endothelial cell','Fibroblast'))
seurat_object <- seurat_object[, sample(1:ncol(seurat_object),round(ncol(seurat_object)/5))] 
seurat_object@meta.data$celltype <- as.character(seurat_object@meta.data$celltype) 
counts <- GetAssayData(seurat_object, slot = 'counts') 
anno <- subset(seurat_object@meta.data, select='celltype') 
gencode <- read_tsv("inferCNV_geneLocate.txt", col_names = c("gene", "chr", "start", "end"))
gencode <- gencode[!duplicated(gencode$gene),]
common_genes <- intersect(gencode$gene, rownames(counts)) 
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts[common_genes,],
                                    annotations_file = anno,
                                    delim="\t",
                                    gene_order_file = "input/#inferCNV_geneLocate.txt",
                                    min_max_counts_per_cell = c(100, +Inf),
                                    ref_group_names = c('Endothelial cell','Fibroblast'))
dir.create("inferCNV")
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1, 
                             out_dir = "inferCNV/", 
                             HMM = F, 
                             denoise = T, 
                             num_threads = 8, 
                             write_expr_matrix=T) 
library(decoupleR)
data <- scRNA_harmony
Idents(data) <- "Group" 
data <- data[, sample(1:ncol(data),round(ncol(data)/5))] 
net <- get_collectri(organism='human', split_complexes=FALSE)
mat <- as.matrix(data@assays$RNA@data)
acts <- run_ulm(mat=mat, net=net, .source='source', .target='target',
                .mor='mor', minsize = 5)
data[['tfsulm']] <- acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition',
              values_from = 'score') %>%
  column_to_rownames('source') %>%
  Seurat::CreateAssayObject(.)
DefaultAssay(object = data) <- "tfsulm"
data <- ScaleData(data)
data@assays$tfsulm@data <- data@assays$tfsulm@scale.data
p.nfkb <- (FeaturePlot(data, features = c("NFKB"),reduction = "tsne",pt.size = 1) &
         scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red'))+
  ggtitle('NFKB activity')
p2.irf1 <- (FeaturePlot(data, features = c("IRF1"),reduction = "tsne",pt.size = 1) &
              scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red'))+
  ggtitle('IRF1 activity')
p2.CIITA <- (FeaturePlot(data, features = c("CIITA"),reduction = "tsne",pt.size = 1) &
              scale_colour_gradient2(low = 'blue', mid = 'white', high = 'red'))+
  ggtitle('CIITA activity')
df <- t(as.matrix(data@assays$tfsulm@data)) %>%
  as.data.frame() %>%
  mutate(cluster = Idents(data)) %>%
  pivot_longer(cols = -cluster, names_to = "source", values_to = "score") %>%
  group_by(cluster, source) %>%
  summarise(mean = mean(score))
tfs <- c("NFKB","AP1","CIITA",
         "STAT1","STAT3","STAT6",
         "IRF1","IRF3","IRF4","IRF5","IRF7",
         "TBX21","GATA3","FOXO1","EOMES","BCL6",
         "PRDM1","HIF1A","MYC","RUNX1","RUNX3","TP53","CTCF",
         "NFATC1","NFATC2","NFATC3","NFATC4")
top_acts_mat <- df %>%
  filter(source %in% tfs) %>%
  pivot_wider(id_cols = 'cluster', names_from = 'source',
              values_from = 'mean') %>%
  column_to_rownames('cluster') %>%
  as.matrix()
palette_length = 50
my_color = colorRampPalette(c("#2067AE", "white","#B11F2B"))(palette_length)
my_breaks <- c(seq(-0.5, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 0.5, length.out=floor(palette_length/2)))
pheatmap(top_acts_mat, border_color = "lightgrey", color=my_color, breaks = my_breaks,
        main="CollecTRI",cluster_rows = FALSE,cluster_cols = FALSE)