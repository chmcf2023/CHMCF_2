library(DoubletFinder)

c(length(Cells(mole1)), length(Cells(normal1)), length(Cells(mole2)), length(Cells(normal2)))

doublet.rate <- c(0.08, 0.08, 0.055, 0.06)

seurdata <- c(mole1, normal1, mole2, normal2)

doub <- mclapply(c(mole1, normal1, mole2, normal2), function(seur) {
	gc(reset = T)
    seur <- subset(seur, subset = nFeature_RNA > 500 & percent.mt < 25)
    seur <- NormalizeData(seur, verbose = FALSE)
    seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    seur <- ScaleData(seur, verbose = F)
    seur <- RunPCA(seur, npcs = 30, verbose = F)
    
    use.dims <- 10
    print("start tsne")
    seur <- RunTSNE(seur, reduction = "pca", dims = 1:use.dims)
    seur <- FindNeighbors(seur, reduction = "pca", dims = 1:use.dims)
    seur <- FindClusters(seur, resolution = 0.4)
    gc(reset = T)
    sweep.data <- paramSweep_v3(seur, PCs = 1:use.dims)
    sweep.stats <- summarizeSweep(sweep.data, GT = F)
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
    homotypic.prop <- modelHomotypic(seur$seurat_clusters)   
    nExp_poi <- round(doublet.rate[i]*ncol(seur)) 
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    gc(reset = T)
    doubletFinder_v3(seur, PCs = 1:use.dims, pN = 0.25, pK = pK_bcmvn, 
                     nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
}) 

DF.classification <- lapply(1:4, function(i) {
    is.dbl = doub[[i]]@meta.data[,length(doub[[i]]@meta.data)]
    names(is.dbl) <- paste0(rownames(doub[[i]]@meta.data), "_", i)
    is.dbl
})

DF.classification <- unlist(DF.classification)


