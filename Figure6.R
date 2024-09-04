library(CellChat)
library(Seurat)
eocc.object <- subset(sce,Group=="EOCC")
eocc.data.input <- GetAssayData(eocc.object, assay = "RNA", slot = "data")
eocc.meta = eocc.object@meta.data[,c("celltype", "Group")]
eocc.cellchat <- createCellChat(object = eocc.data.input)  
eocc.cellchat <- addMeta(eocc.cellchat, meta = eocc.meta)
eocc.cellchat <- setIdent(eocc.cellchat, ident.use = "celltype")  
eocc.cellchat@DB <- CellChatDB.human 
eocc.cellchat <- subsetData(eocc.cellchat) 
future::plan("multisession", workers = 4) 
eocc.cellchat <- identifyOverExpressedGenes(eocc.cellchat)
eocc.cellchat <- identifyOverExpressedInteractions(eocc.cellchat)
eocc.cellchat <- computeCommunProb(eocc.cellchat) 
eocc.cellchat <- filterCommunication(eocc.cellchat, min.cells = 10)
eocc.cellchat <- computeCommunProbPathway(eocc.cellchat)
eocc.cellchat <- aggregateNet(eocc.cellchat)
eocc.cellchat <- netAnalysis_computeCentrality(eocc.cellchat, slot.name = "netP") 
locc.object <- subset(sce,Group=="LOCC")
locc.data.input <- GetAssayData(locc.object, assay = "RNA", slot = "data")
locc.meta = locc.object@meta.data[,c("celltype", "Group")]
locc.cellchat <- createCellChat(object = locc.data.input)
locc.cellchat <- addMeta(locc.cellchat, meta = locc.meta)
locc.cellchat <- setIdent(locc.cellchat, ident.use = "celltype")
locc.cellchat@DB <- CellChatDB.human 
locc.cellchat <- subsetData(locc.cellchat) 
future::plan("multisession", workers = 4) 
locc.cellchat <- identifyOverExpressedGenes(locc.cellchat)
locc.cellchat <- identifyOverExpressedInteractions(locc.cellchat)
locc.cellchat <- computeCommunProb(locc.cellchat) 
locc.cellchat <- filterCommunication(locc.cellchat, min.cells = 10)
locc.cellchat <- computeCommunProbPathway(locc.cellchat)
locc.cellchat <- aggregateNet(locc.cellchat)
locc.cellchat <- netAnalysis_computeCentrality(locc.cellchat, slot.name = "netP") 
object.list <- list(EOCC = eocc.cellchat, LOCC = locc.cellchat) 
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
compareInteractions(cellchat, show.legend = F, group = c(1,2))
compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
netVisual_heatmap(cellchat)
netVisual_heatmap(cellchat, measure = "weight")
rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE) 
rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE) 
pathway.union <- c("MHC-I","MHC-II","CCL","CXCL","COMPLEMENT","IFN-II","TGFb","TNF",
                   "IL16","CD6","CD70","CD96","CD99","CD137","CSF","PVR","SPP1","ICAM")
netAnalysis_signalingRole_heatmap(object.list[[1]],
                                        pattern = "outgoing", 
                                        signaling = pathway.union,
                                        title = names(object.list)[1],
                                        width = 5,
                                        height = 6)
netAnalysis_signalingRole_heatmap(object.list[[2]],
                                        pattern = "outgoing", 
                                        signaling = pathway.union,
                                        title = names(object.list)[2],
                                        width = 5,
                                        height = 6)
netAnalysis_signalingRole_heatmap(object.list[[1]],
                                        pattern = "incoming", 
                                        signaling = pathway.union,
                                        title = names(object.list)[1],
                                        width = 5, height = 6,
                                        color.heatmap = "GnBu")
netAnalysis_signalingRole_heatmap(object.list[[2]],
                                        pattern = "incoming", 
                                        signaling = pathway.union,
                                        title = names(object.list)[2],
                                        width = 5, height = 6,
                                        color.heatmap = "GnBu")
netVisual_bubble(cellchat,
                 sources.use = 1,
                 targets.use = c(4,5),
                 comparison = c(1, 2),
                 angle.x = 45)
netVisual_bubble(cellchat,
                 sources.use = c(4,5),
                 targets.use = 1,
                 comparison = c(1, 2),
                 angle.x = 45)
netVisual_bubble(cellchat,
                 sources.use = c(3),
                 targets.use = c(4,5),
                 comparison = c(1, 2),
                 angle.x = 45)
netVisual_bubble(cellchat,
                sources.use = c(4,5),
                targets.use = c(3),
                comparison = c(1, 2),
                angle.x = 45)
netVisual_bubble(cellchat,
                 sources.use = c(4,5),
                 targets.use = c(4,5),
                 comparison = c(1, 2),
                 angle.x = 45)