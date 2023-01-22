library(dplyr)
library(Seurat)
library(ggplot2)

create_seurat <- function(filename,project,feature,mt){
    control1 <- Read10X(data.dir = filename)
    C1data <- CreateSeuratObject(counts = control1 , project = project , min.cells = 3 , min.features = 200) 
    C1data[["percent.mt"]] <- PercentageFeatureSet(C1data, pattern = "^MT-")
    C1data[["percent.rb"]] <- PercentageFeatureSet(C1data, pattern = "^RP[LS]")
    C1data <- subset(C1data,subset = nFeature_RNA > feature & percent.mt < mt)
}


mole1 <- create_seurat(filename = 'data/mole1/filtered_feature_bc_matrix/','mole1',500,100)
normal1 <- create_seurat(filename = 'data/normal1/filtered_feature_bc_matrix/','normal1',500,100)
mole2 <- create_seurat(filename = 'data/mole2/filtered_feature_bc_matrix/','mole2',500,100)
normal2 <- create_seurat(filename = 'data/normal2/filtered_feature_bc_matrix/','normal2',500,100)

alldata <- merge(x = mole1, y = c(normal1, mole2, normal2))
VlnPlot(alldata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), ncol = 4, pt.size = 0.01)
alldata <- subset(alldata, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & percent.mt < 25)
DF.classification <- readRDS("runtime/env/DF.classification.rds")
alldata <- AddMetaData(alldata, DF.classification, col.name = "DF.classification")
singletdata <- subset(alldata, subset = DF.classification == "Singlet")

vil.list <- SplitObject(singletdata, split.by = "orig.ident")

for (i in 1:length(vil.list)) {
    vil.list[[i]] <- NormalizeData(vil.list[[i]], verbose = FALSE)
    vil.list[[i]] <- FindVariableFeatures(vil.list[[i]], selection.method = "vst", 
                                          nfeatures = 2000, verbose = FALSE)
}

vil.anchors <- FindIntegrationAnchors(object.list = vil.list, dims = 1:30)
vil.combined <- IntegrateData(anchorset = vil.anchors, dims = 1:30)
DefaultAssay(vil.combined) <- "integrated"

vil.combined <- ScaleData(vil.combined, verbose = FALSE)
vil.combined <- RunPCA(vil.combined, npcs = 30, verbose = FALSE)

use.dims <- 15

vil.combined <- RunTSNE(vil.combined, reduction = "pca", dims = 1:use.dims)
vil.combined <- FindNeighbors(vil.combined, reduction = "pca", dims = 1:use.dims)
vil.combined <- FindClusters(vil.combined, resolution = 0.6)

VlnPlot(vil.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.01)

celltype.cluster <- list(
    EVT = c(0,10,1,2,5,3,8,12), 
    STB = c(13,14),    
    CTB = c(4,7,9),          
    cCTB = c(11),       
    HC = c(15),    
    MAC = c(6,16),
    FB = c(18),
    DL = c(17)       
)
new.cluster.ids <- names(sort(unlist(celltype.cluster)))
vil.combined$celltype <- new.cluster.ids[vil.combined$seurat_clusters]
vil.combined$celltype.1 <- stringr::str_split(vil.combined$celltype, "[0-9]+$", simplify = T)[, 1]
vil.combined$celltype.1 <- factor(vil.combined$celltype.1, levels = c("CTB", "EVT", "STB", "EVTP", "HC", "FB", "DL"))
vil.combined$sample <- ifelse(vil.combined$orig.ident %in% c("mole1", "mole2"), "AnCHM", "Normal")
vil.combined$sample <- factor(vil.combined$sample, levels = c("Normal", "AnCHM"))
vil.combined$case <- ifelse(vil.combined$orig.ident %in% c("mole1", "normal1"), "Case1", "Case2")
DimPlot(vil.combined, reduction = "tsne", group.by = "celltype.1", label = T)
#find marker
vil.markers <- FindAllMarkers(vil.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

vil.top10 <- vil.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)


vil.clus.expr.avg <-  AverageExpression(vil.combined)
vil.celltype.expr.avg <- AverageExpression(vil.combined, group.by = "celltype.1")

gene <- c("THY1","COL1A1","COL1A2", "ACTA2","DLK1","EGFL6", 
          "CD14","AIF1","CD163", "FCGRT", "S100A8",
          "ITGA2",
          "CGA","ERVFRD-1","CYP19A1","CGB3",
          "LAIR2", "HLA-G","NOTUM", "MMP2", "HTRA4", "DIO2", "ERBB2", "HAPLN3", "ITGA5",
		  "PAGE4","DUSP9","PARP1","PEG10","SLC27A2",
		  "GATA3", "PERP", "KRT7")
{
    showdata <- vil.clus.expr.avg$RNA
    showdata <- showdata[rownames(showdata) %in% gene,]
    showdata <- showdata[match(gene, rownames(showdata)),]
    showdata <- showdata[, unlist(celltype.cluster) + 1]
}
pheatmap::pheatmap(showdata, scale = "row", cluster_rows = F, cluster_cols = F, angle_col = 45,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),border_color = "white")

{
    showdata <- vil.celltype.expr.avg$RNA
    showdata <- showdata[rownames(showdata) %in% gene,]
    showdata <- showdata[match(gene, rownames(showdata)),]
    showdata <- showdata[, levels(vil.combined$celltype.1)[1:6]]
}
pheatmap::pheatmap(showdata, scale = "row", cluster_rows = F, cluster_cols = F, angle_col = 45,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),border_color = "white")

Seurat::DoHeatmap(vil.combined, features = gene, group.by = 'celltype.1')
plot <- DoHeatmap(object = vil.combined, features = gene, group.by = 'celltype.1') + NoLegend()
ggplot2::ggsave(filename = "feature.pdf", plot = plot, device = "pdf", dpi = 300, width = 10, height = 8) # can add additional parameters similar to png 

prop.dataset.count <- as.data.frame(table(vil.combined$orig.ident, vil.combined$celltype.1)) %>% 
    group_by(Var1) %>% 
    group_map(function(df, celltype) {
        print(celltype)
        vec <- df$Freq
        names(vec) <- df$Var2
        vec
    }, .keep = T)
names(prop.dataset.count) <- c("mole1", "mole2", "normal1", "normal2")

###### CTB #######
vil.ctb <- subset(vil.combined, subset = celltype %in% c("CTB2", "CTB3", "EVTP"))
vil.ctb.reclu <- NormalizeData(vil.ctb, verbose = FALSE)
vil.ctb.reclu <- FindVariableFeatures(vil.ctb.reclu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
DefaultAssay(vil.ctb.reclu) <- "integrated"
vil.ctb.reclu <- ScaleData(vil.ctb.reclu, verbose = FALSE)
vil.ctb.reclu <- RunPCA(vil.ctb.reclu, npcs = 30, verbose = FALSE)

use.dims <- 6

vil.ctb.reclu <- RunTSNE(vil.ctb.reclu, reduction = "pca", dims = 1:use.dims)
vil.ctb.reclu <- FindNeighbors(vil.ctb.reclu, reduction = "pca", dims = 1:use.dims)
vil.ctb.reclu <- FindClusters(vil.ctb.reclu, resolution = 0.1)
DefaultAssay(vil.ctb.reclu) <- "RNA"

Idents(vil.ctb.reclu) <- "seurat_clusters"

vil.ctb.reclu$samplename <- paste(vil.ctb.reclu$case, vil.ctb.reclu$sample, sep = " ")
vil.ctb.reclu$subtype <- paste0("CTB", as.numeric(vil.ctb.reclu$seurat_clusters))
vil.ctb.reclu$subtype <- factor(vil.ctb.reclu$subtype, levels = paste0("CTB", 1:2))
vil.ctb.reclu$sample <- factor(vil.ctb.reclu$sample, levels = c("Normal", "AnCHM"))
vil.ctb.reclu$subtype.by.sample <- factor(paste0(vil.ctb.reclu$subtype, vil.ctb.reclu$sample),
										  levels = expand_grid(
										  	levels(vil.ctb.reclu$subtype), 
										  	levels(vil.ctb.reclu$sample)) %>% apply(1, function(x) paste0(x[1], x[2])))


library(clusterProfiler)
library(org.Hs.eg.db)

ctb.reclu.marker <- FindAllMarkers(vil.ctb.reclu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
ctbGO.ctb.reclu.marker <- mclapply(levels(ctb.reclu.marker$cluster), function(c) {
    g <- ctb.reclu.marker[ctb.reclu.marker$cluster == c, "gene"]
    # g <- g[!g %in% c("CDKN1C", "PHLDA2")]
    tmp <- bitr(g, fromType = "SYMBOL",
                toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) 
    enrichGO(gene = tmp$ENTREZID,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             pAdjustMethod = "BH",
             ont = "BP",
             pvalueCutoff = 0.5,
             qvalueCutoff = 0.5,
             readable = T
    )
})

{i <- 2
	dotplot(ctbGO.ctb.reclu.marker[[i]], font.size = 8)
}
vil.ctb.reclu.expr.avg <- AverageExpression(vil.ctb.reclu, group.by = c("subtype"))
gene <- ctb.reclu.marker.top10$gene
{
    showdata <- vil.ctb.reclu.expr.avg$RNA
    showdata <- showdata[rownames(showdata) %in% gene,]
    showdata <- showdata[match(gene, rownames(showdata)),]
}
pheatmap::pheatmap(showdata, scale = "row", cluster_rows = F, cluster_cols = F, angle_col = 45,
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),border_color = "white")

library(GSEABase)
ctb.gsea.res <- lapply(0:1, function(x) {
	submarker = ctb.reclu.marker[ctb.reclu.marker$cluster == x & ctb.reclu.marker$p_val_adj < 0.05, ]
	genelist = submarker[order(submarker$avg_log2FC, decreasing = T), "avg_log2FC"] %>% 
		setNames(submarker[order(submarker$avg_log2FC, decreasing = T), "gene"])
	GSEA(genelist, 
		 TERM2GENE = read.gmt("~/share/msigdb/c5.go.bp.v7.5.1.symbols.gmt"),
		 minGSSize = 10, pvalueCutoff = 0.99, pAdjustMethod = "BH")
})
library(DOSE)
library(GOSemSim)
GO.key <- keys(org.Hs.eg.db, keytype = "GO")
metab.GS <- sapply(c(GLYCOLYSIS = "GO:0006096", TCA_Cycle = "GO:0006099", OXPHOS = "GO:0006119", 
					 HYPOXIA = "GO:0001666", Electron = "GO:0022900",
					 Nucleotide = "GO:0009117", FattyAcid = "GO:0019395", 
					  Glutamine = "GO:0006541", Pentose = "GO:0006098"), function(i) {
	unique(AnnotationDbi::select(org.Hs.eg.db, keys = i, columns = c("ENTREZID","SYMBOL","ONTOLOGY"), keytype="GO")$SYMBOL)
}, USE.NAMES = T, simplify = F)

minexprpct <- function(object, features, min.pct = 10, group.by = NA) {
	dat = FetchData(object, vars = c(features), slot = "counts")
	if (is.na(group.by)) {
		cnum = nrow(dat)
		return(colnames(dat[, colSums(dat > 0) / cnum * 100 > min.pct]))
	}
	group = object@meta.data[, group.by]
	names(group) <- rownames(object@meta.data)
	littleexpr <- Reduce(intersect, lapply(unique(group), function(x) {
		dat.1 = dat[names(group)[group == x], ]
		cnum = nrow(dat.1)
		colnames(dat.1[, colSums(dat.1 > 0) / cnum * 100 <= min.pct])
	}))
	return(setdiff(colnames(dat), littleexpr))
}

minexprpct(object = vil.ctb.reclu, features = metab.GS$GLYCOLYSIS, min.pct = 10, group.by = "subtype.by.sample")

library(GSVA)
library(GSEABase)
library(limma)

ctb.reclu.bulk <- do.call(cbind, mclapply(levels(vil.ctb.reclu$subtype.by.sample), function(x) {     # pseudobulk
	subtype = subset(vil.ctb.reclu, subset = subtype.by.sample == x)
	Idents(subtype) <- "subtype.by.sample"
	lapply(1:floor(length(Cells(subtype)) / 10), function(i) {
		set.seed(i)
		rdmcell = sample(Cells(subtype), size = 50)
		AverageExpression(subset(subtype, cells = rdmcell))$RNA
	}) %>% as.data.frame() %>% data.table::setnames(paste0(x, "_", 1:floor(length(Cells(subtype)) / 10)))
}))

gs.c2.kegg <- getGmt("~/projects/share/msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
gs.c5.go <- getGmt("~/projects/share/msigdb/c5.go.bp.v7.5.1.symbols.gmt")
gs.c6 <-  getGmt("~/projects/share/msigdb/c6.all.v7.5.1.symbols.gmt")
gs.h <-  getGmt("~/projects/share/msigdb/h.all.v7.5.1.symbols.gmt")

ctb.reclu.gsva.c2.kegg <- GSVA::gsva(expr=as.matrix(ctb.reclu.bulk), gset.idx.list=gs.c2.kegg, kcdf="Gaussian", parallel.sz=30, method = "gsva")
ctb.reclu.gsva.c5.go <- GSVA::gsva(expr=as.matrix(ctb.reclu.bulk), gset.idx.list=gs.c5.go, kcdf="Gaussian", parallel.sz=30, method = "gsva")
ctb.reclu.gsva.c6 <- GSVA::gsva(expr=as.matrix(ctb.reclu.bulk), gset.idx.list=gs.c6, kcdf="Gaussian", parallel.sz=30, method = "gsva")
ctb.reclu.gsva.h <- GSVA::gsva(expr=as.matrix(ctb.reclu.bulk), gset.idx.list=gs.h, kcdf="Gaussian", parallel.sz=30, method = "gsva")


ctb.reclu.gsva <- ctb.reclu.gsva.c5.go
grouP <- stringr::str_split(colnames(ctb.reclu.bulk), "_", simplify = T)[,1] %>% as.factor()
desigN <- model.matrix(~ grouP + 0)

rownames(desigN) <- colnames(ctb.reclu.bulk)
colnames(desigN) <- unique(grouP)
comparE <- makeContrasts(CTB2AnCHM - CTB2Normal, levels=desigN)
fiT <- lmFit(ctb.reclu.gsva, desigN)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
gsvaDiff <- topTable(fiT3, coef=1, number=Inf, adjust.method="BH", p.value = 0.05)
gsvaDiff$pathway <- rownames(gsvaDiff)

pheatmap::pheatmap(ctb.reclu.gsva[gsvaDiff[gsvaDiff$adj.P.Val < 0.01 & (gsvaDiff$logFC > 0.5 | gsvaDiff$logFC < -0.5), "pathway"], grep("CTB2", colnames(ctb.reclu.gsva))], show_colnames = F)
pheatmap::pheatmap(ctb.reclu.gsva[gsvaDiff[gsvaDiff$adj.P.Val < 0.01 & (gsvaDiff$logFC > 1 | gsvaDiff$logFC < -1), "pathway"], grep("CTB2", colnames(ctb.reclu.gsva))], 
				   show_colnames = F, annotation_col = annotation, fontsize_row = 8)

annotation = data.frame(sample = str_split(string =grep("CTB2", colnames(ctb.reclu.gsva), value = T), pattern = "_", simplify = T)[,1] )
rownames(annotation) <- grep("CTB2", colnames(ctb.reclu.gsva), value = T)


ctb.prop.dataset.count <- as.data.frame(table(vil.ctb.reclu$orig.ident, vil.ctb.reclu$subtype)) %>% 
    group_by(Var1) %>% 
    group_map(function(df, celltype) {
        print(celltype)
        vec <- df$Freq 
        names(vec) <- df$Var2 
        vec
    }, .keep = T)
names(ctb.prop.dataset.count) <- c("C1M", "C2M", "C1N", "C2N")
ctb.prop.dataset.count <- as.data.frame(t(as.data.frame(ctb.prop.dataset.count)))

chisqres.1 <- chisq.test(ctb.prop.dataset.count[c(1,3),], correct = F)
chisqres.2 <- chisq.test(ctb.prop.dataset.count[c(2,4),], correct = F)

chisqres.prop.1 <- as.data.frame(t(sapply(1:ncol(ctb.prop.dataset.count), function(i) {
    dat = data.frame(ifelse(ncol(ctb.prop.dataset.count) > 2, 
    						rowSums(ctb.prop.dataset.count[, -i]), 
    						ctb.prop.dataset.count[, -i]), 
    				 ctb.prop.dataset.count[, i])
    
    colnames(dat) <- c(paste0("non-ctb", i), paste0("ctb", i))
    dat = rbind(colSums(dat[1,]), colSums(dat[3,]))
    rownames(dat) = c("AnCHM", "Normal")
    res = chisq.test(dat, correct = T)
    c(res$statistic, pval = res$p.value)
})))
chisqres.prop.2 <- as.data.frame(t(sapply(1:ncol(ctb.prop.dataset.count), function(i) {
    dat = data.frame(ifelse(ncol(ctb.prop.dataset.count) > 2, 
    						rowSums(ctb.prop.dataset.count[, -i]), 
    						ctb.prop.dataset.count[, -i]), 
    				 ctb.prop.dataset.count[, i])
    
    colnames(dat) <- c(paste0("non-ctb", i), paste0("ctb", i))
    dat = rbind(colSums(dat[2,]), colSums(dat[4,]))
    rownames(dat) = c("AnCHM", "Normal")
    
    res = chisq.test(dat, correct = T)  
    c(res$statistic, pval = res$p.value)
})))
chisqres.prop.1$adj.pval <- p.adjust(chisqres.prop.1$pval, method = "bonferroni")
chisqres.prop.2$adj.pval <- p.adjust(chisqres.prop.2$pval, method = "bonferroni")
chisqres.prop.1$cluster <- paste0("CTB", 1:ncol(ctb.prop.dataset.count))
chisqres.prop.2$cluster <- paste0("CTB", 1:ncol(ctb.prop.dataset.count))

library(scProportionTest)
prop_test <- sc_utils(vil.ctb.reclu)
prop_permu.1 <- permutation_test(
    prop_test, cluster_identity = "subtype",
    sample_1 = "normal1", sample_2 = "mole1",
    sample_identity = "orig.ident", n_permutations = 1000
)
prop_permu.2 <- permutation_test(
    prop_test, cluster_identity = "subtype",
    sample_1 = "normal2", sample_2 = "mole2",
    sample_identity = "orig.ident", n_permutations = 1000
)


chisqres.prop.1 <- cbind(chisqres.prop.1, 
   prop_permu.1@results$permutation[, c("obs_log2FD", "boot_mean_log2FD", "boot_CI_2.5" ,"boot_CI_97.5")])
chisqres.prop.2 <- cbind(chisqres.prop.2, 
   prop_permu.2@results$permutation[, c("obs_log2FD", "boot_mean_log2FD", "boot_CI_2.5" ,"boot_CI_97.5")])

log2FC_threshold = log2(1.5)
chisqres.prop <- rbind(chisqres.prop.1, chisqres.prop.2)
chisqres.prop$case <- c(rep("1", ncol(ctb.prop.dataset.count)), rep("2", ncol(ctb.prop.dataset.count)))
chisqres.prop$fc <- ifelse(chisqres.prop$adj.pval < 0.05 & chisqres.prop$obs_log2FD > log2FC_threshold, "UP", 
                     ifelse(chisqres.prop$adj.pval < 0.05 & chisqres.prop$obs_log2FD < -log2FC_threshold, "DOWN", "~"))
chisqres.prop$sig <- ifelse(chisqres.prop$adj.pval < 0.05, "sig", "n.s.")
chisqres.prop$neg.log.p <- -log2(chisqres.prop$adj.pval)


######## STB ###############
vil.stb <- subset(vil.combined, subset = celltype.1 == "STB")

vil.stb.reclu <- NormalizeData(vil.stb, verbose = FALSE)
vil.stb.reclu <- FindVariableFeatures(vil.stb.reclu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
DefaultAssay(vil.stb.reclu) <- "integrated"
vil.stb.reclu <- ScaleData(vil.stb.reclu, verbose = FALSE)
vil.stb.reclu <- RunPCA(vil.stb.reclu, npcs = 30, verbose = FALSE)
use.dims <- 5
vil.stb.reclu <- RunTSNE(vil.stb.reclu, reduction = "pca", dims = 1:use.dims)
vil.stb.reclu <- FindNeighbors(vil.stb.reclu, reduction = "pca", dims = 1:use.dims)
vil.stb.reclu <- FindClusters(vil.stb.reclu, resolution = 0.2)
DefaultAssay(vil.stb.reclu) <- "RNA"
vil.stb.reclu$subtype <- factor(paste0("STB", c(3,1,2,4))[vil.stb.reclu$seurat_clusters], levels = paste0("STB", 1:4))

vil.stb.reclu.dirty <- vil.stb.reclu
vil.stb.reclu <- subset(vil.stb.reclu, subset = subtype != "STB4")
Idents(vil.stb.reclu) <- "subtype"
Idents(vil.stb.reclu) <- "seurat_clusters"
stb.reclu.markers <- FindAllMarkers(vil.stb.reclu, only.pos = T, min.pct = 0.2, logfc.threshold = 0.5)
stb.marker.top10 <- stb.reclu.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

stbGO.stb.reclu.marker <- lapply(paste0("STB", 1:2), function(c) {
	enrichGO(bitr(
		stb.reclu.markers[stb.reclu.markers$cluster == c, ]$gene, 
		fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    ont = "BP",
    pvalueCutoff = 0.5,
    qvalueCutoff = 0.5,
    readable = T
)})


vil.stb1.MvN$log.p.adj <- -log10(vil.stb1.MvN$p_val_adj)
threshold = c(p = 0.05, fc = 0.6)
vil.stb1.MvN$sig <- ifelse(
	vil.stb1.MvN$p_val_adj < threshold["p"] & vil.stb1.MvN$avg_log2FC > threshold["fc"], 
	"Up", ifelse(
		vil.stb1.MvN$p_val_adj < threshold["p"] & vil.stb1.MvN$avg_log2FC < -1*threshold["fc"], 
		"Down","No Sig"))



library(GSVA)
library(GSEABase)
library(limma)


stb.reclu.bulk <- do.call(cbind, mclapply(levels(vil.stb.reclu$subtype.by.sample), function(x) {     # pseudobulk
	subtype = subset(vil.stb.reclu, subset = subtype.by.sample == x)
	Idents(subtype) <- "subtype.by.sample"
	lapply(1:floor(length(Cells(subtype)) / 10), function(i) {
		set.seed(i)
		rdmcell = sample(Cells(subtype), size = 50)
		AverageExpression(subset(subtype, cells = rdmcell))$RNA
	}) %>% as.data.frame() %>% data.table::setnames(paste0(x, "_", 1:floor(length(Cells(subtype)) / 10)))
}))

gs.c2.kegg <- getGmt("~/projects/share/msigdb/c2.cp.kegg.v7.5.1.symbols.gmt")
gs.c5.go <- getGmt("~/projects/share/msigdb/c5.go.bp.v7.5.1.symbols.gmt")
gs.c6 <-  getGmt("~/projects/share/msigdb/c6.all.v7.5.1.symbols.gmt")
gs.h <-  getGmt("~/projects/share/msigdb/h.all.v7.5.1.symbols.gmt")

ctb.reclu.gsva.c2.kegg <- GSVA::gsva(expr=as.matrix(stb.reclu.bulk), gset.idx.list=gs.c2.kegg, kcdf="Gaussian", parallel.sz=30, method = "gsva")
stb.reclu.gsva.c5.go <- GSVA::gsva(expr=as.matrix(stb.reclu.bulk), gset.idx.list=gs.c5.go, kcdf="Gaussian", parallel.sz=30, method = "gsva")
stb.reclu.gsva.c6 <- GSVA::gsva(expr=as.matrix(stb.reclu.bulk), gset.idx.list=gs.c6, kcdf="Gaussian", parallel.sz=30, method = "gsva")
stb.reclu.gsva.h <- GSVA::gsva(expr=as.matrix(stb.reclu.bulk), gset.idx.list=gs.h, kcdf="Gaussian", parallel.sz=30, method = "gsva")


stb.reclu.gsva <- stb.reclu.gsva.c5.go
grouP <- stringr::str_split(colnames(stb.reclu.bulk), "_", simplify = T)[,1] %>% as.factor()
desigN <- model.matrix(~ grouP + 0)

rownames(desigN) <- colnames(stb.reclu.bulk)
colnames(desigN) <- unique(grouP)
comparE <- makeContrasts(STB2Normal - STB2AnCHM, levels=desigN)
fiT <- lmFit(stb.reclu.gsva, desigN)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
gsvaDiff <- topTable(fiT3, coef=1, number=Inf, adjust.method="BH", p.value = 0.05)
gsvaDiff$pathway <- rownames(gsvaDiff)


####### EVT #############
vil.evt <- subset(vil.combined, subset = celltype.1 %in% c( "EVT"))
vil.evt.reclu <- NormalizeData(vil.evt, verbose = FALSE)
vil.evt.reclu <- FindVariableFeatures(vil.evt.reclu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
DefaultAssay(vil.evt.reclu) <- "integrated"
vil.evt.reclu <- ScaleData(vil.evt.reclu, verbose = FALSE)
vil.evt.reclu <- RunPCA(vil.evt.reclu, npcs = 30, verbose = FALSE)

use.dims <- 10

vil.evt.reclu <- RunTSNE(vil.evt.reclu, reduction = "pca", dims = 1:use.dims)
vil.evt.reclu <- FindNeighbors(vil.evt.reclu, reduction = "pca", dims = 1:use.dims)
vil.evt.reclu <- FindClusters(vil.evt.reclu, resolution = 0.1)
DefaultAssay(vil.evt.reclu) <- "RNA"


vil.evt.reclu$subtype <- paste0("EVT", c(2, 1, 4, 3, 5))[vil.evt.reclu$seurat_clusters]
VlnPlot(vil.evt.reclu, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.01, group.by = 'subtype')
VlnPlot(vil.evt.reclu, features = c("HLA-A", "HLA-B"), ncol = 2, pt.size = 0.01, group.by = 'subtype')
vil.evt.reclu.dirty <- vil.evt.reclu
vil.evt.reclu <- subset(vil.evt.reclu.dirty, subset = subtype %in% c("EVT1", "EVT2", "EVT3"))
vil.evt.reclu$seurat_clusters <- factor(vil.evt.reclu$seurat_clusters, levels = 0:2)
Idents(vil.evt.reclu) <- "seurat_clusters"

vil.evt.reclu$samplename <- paste(vil.evt.reclu$case, vil.evt.reclu$sample, sep = " ")

Idents(vil.evt.reclu) <- "subtype"
evt.reclu.marker <- FindAllMarkers(vil.evt.reclu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
Idents(vil.evt.reclu) <- "seurat_clusters"
evtGO.evt.reclu.marker <- mclapply(paste0("EVT", 1:3), function(c) {
    g <- evt.reclu.marker[evt.reclu.marker$cluster == c, "gene"]
    tmp <- bitr(g, fromType = "SYMBOL",
                toType = c("ENTREZID"), OrgDb = org.Hs.eg.db) 
    enrichGO(gene = tmp$ENTREZID,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             pAdjustMethod = "BH",
             ont = "BP",
             pvalueCutoff = 0.5,
             qvalueCutoff = 0.5,
             readable = T
    )
})

{i <- 1
	dotplot(evtGO.evt.reclu.marker[[i]], font.size = 8)
}


vil.evt.reclu.expr.avg <-  AverageExpression(vil.evt.reclu, group.by = c("subtype", "sample"))


gene <- evt.reclu.marker.top10$gene
gene <- c("AOC1", "PAPPA", "TIMP3", "PRG2", "SERPINE1", "FN1", "CSH1",
		  "LAIR2", "HTRA4", "MMP2","CLIC3","MMP12",
		  "RRM2", "TOP2A", "MKI67", "HMGB1", "CENPF")
gene <- c( "RRM2", "MKI67","HMGB2", "CENPF", "TOP2A", 
		   "MMP2", "MMP12", "LAIR2", "HTRA4", "TAC3","PRG2", "JAM2", "SERPINE1", "FN1", "FBN1")
{
    showdata <- vil.evt.reclu.expr.avg$RNA
    showdata <- showdata[rownames(showdata) %in% gene,]
    showdata <- showdata[match(gene, rownames(showdata)),]
    # showdata <- showdata[, colnames(showdata) != "2_AnCHM"]
}
pheatmap::pheatmap(showdata, scale = "row", cluster_rows = F, cluster_cols = F, angle_col = 45, 
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),border_color = "white")


library(clusterProfiler)
library(GSVA)
evt.expr <- as.data.frame(vil.evt.reclu@assays$RNA@data)

s2e <- bitr(rownames(evt.expr), 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
evt.expr <- evt.expr[s2e$SYMBOL, ]
rownames(evt.expr) <- s2e$ENTREZID

c2.kegg.gmt <- read.gmt("runtime/ref/c2.cp.kegg.v7.5.1.entrez.gmt")
c2.kegg.list = split(c2.kegg.gmt$gene, c2.kegg.gmt$term)
evt.expr <- as.matrix(evt.expr)
evt.gsva.kegg.res <- gsva(evt.expr, c2.kegg.list, kcdf="Gaussian", method = "gsva", parallel.sz = 20)

vil.evt.reclu.kegg <- vil.evt.reclu
vil.evt.reclu.kegg[["KEGG"]] <- CreateAssayObject(evt.gsva.kegg.res)
# vil.evt.reclu.kegg <- AddMetaData(vil.evt.reclu, data.frame(t(evt.gsva.kegg.res), stringsAsFactors = F))
FeaturePlot(vil.evt.reclu.kegg, features = "KEGG_RIBOSOME", reduction = 'tsne')
FeaturePlot(vil.evt.reclu.kegg, features = "KEGG_ASTHMA", reduction = 'tsne')

library(limma)
calc.gsva.DE <- function(exprSet, meta, compare = NULL){
  allDiff <- list()
  design <- model.matrix(~factor(meta))
  colnames(design) <- levels(factor(meta))
  rownames(design) <- colnames(exprSet)
  
  fit <- lmFit(exprSet, design)
  
  if(length(unique(meta)) == 2){
    if(is.null(compare)) {
      stop("There are 2 Groups, Please set compare value")
    }
    contrast.matrix <- makeContrasts(contrasts = compare, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix) 
    fit2 <- eBayes(fit2)
    tempOutput = topTable(fit2,adjust='fdr', coef = 1, number=Inf)
    allDiff[[compare]] = na.omit(tempOutput)
  } else if(length(unique(meta)) > 2) {
    for(g in colnames(design)) {
      fm = ""
      for(gother in colnames(design)[which(!colnames(design) %in% g)]){
        fm = paste0(fm, "+", gother)
      }
      fm = paste0(g, "VsOthers = ", g, "-(", substring(fm, 2), ")/", ncol(design)-1)
      contrast.matrix <- makeContrasts(contrasts = fm, levels = design)
      fit2 <- contrasts.fit(fit, contrast.matrix) 
      fit2 <- eBayes(fit2)
      allDiff[[g]] = topTable(fit2, adjust = 'fdr', coef = 1, number = Inf)
    }
  } else {
    stop("error only have one group")
  }
  return(allDiff)
}

gsva.DE <- calc.gsva.DE(
 	exprSet = evt.gsva.kegg.res,
 	meta = paste0("newEVT.", vil.evt.reclu$seurat_clusters)
)

idiff <- gsva.DE[[1]]
idiff$PATH_NAME <- rownames(idiff)
idiff.up <- idiff[idiff$logFC > 0,]
idiff.up$ifupdown <- "up"
idiff.down <- idiff[idiff$logFC < 0,]
idiff.down$ifupdown <- "down"
idifftop <- rbind(
	idiff.up[order(idiff.up$t, decreasing = T)[1:10],],
	idiff.down[order(idiff.down$t, decreasing = F)[1:10],]
)
Padj_threshold <- 0.05
df$group <- ifelse(idifftop$logFC > 0 & idifftop$adj.P.Val < Padj_threshold, "up", 
				   ifelse(idifftop$logFC < 0 & idifftop$adj.P.Val < Padj_threshold, "down", "non-sig"))
idifftop$hjust = ifelse(idifftop$t > 0, 1, 0)
idifftop$nudge_y = ifelse(idifftop$t > 0, -0.1, 0.1)
idifftop <- idifftop[order(idifftop$t),]
idifftop$PATH_NAME <- factor(idifftop$PATH_NAME, levels = idifftop$PATH_NAME)
limt <- max(abs(idifftop$t))
ggplot(idifftop, aes(PATH_NAME, t, fill = ifupdown)) + 
    geom_bar(stat = 'identity',alpha = 0.7) + 
    scale_fill_manual(breaks=c("down","noSig","up"), values = c("#008020","grey","#08519C"))+
    geom_text(data = idifftop, aes(label = PATH_NAME, y = nudge_y),
              nudge_x =0, nudge_y = 0, hjust = idifftop$hjust, size = 2.5)+
    labs(x = "KEGG pathways", y = "t value of GSVA score", title = "EVT Cluster 0") +
    scale_y_continuous(limits=c(-limt,limt)) +
    coord_flip() + 
    theme_bw() +
    theme(panel.grid =element_blank()) +
    theme(panel.border = element_rect(size = 0.6)
          #panel.border = element_blank()
    ) +
    theme(plot.title = element_text(hjust = 0.5,size = 12),
          axis.text.y = element_blank(),
          axis.title = element_text(hjust = 0.5,size = 12),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = limt
    )

DotPlot(vil.evt.reclu, features = c("EIF4E", "FOS", "MTA3", "ID1", "STAT1", "ASCL2", "GCM1"),  dot.scale = 8, 
        group.by = "subtype.by.sample", ) + RotatedAxis() + NoAxes() + 
	theme(
		text = element_text(size = 16),
		# legend.position="bottom", 
		panel.border = element_rect(fill=NA,color="black", linetype="solid"),
		axis.ticks.x = element_line(colour = "black"),
		axis.ticks.y = element_line(colour = "black"),
		axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
		axis.text.y = element_text(),
		axis.title.x = element_blank(),
		axis.title.y = element_blank()
	) 

DotPlot(vil.evt.reclu, features = c("ID1", "STAT1", "ASCL2", "GCM1"), cols = c("grey", "firebrick"), dot.scale = 8, 
        group.by = "subtype.by.sample") + RotatedAxis() + NoAxes() + 
	theme(
		text = element_text(size = 16),
		# legend.position="bottom", 
		panel.border = element_rect(fill=NA,color="black", linetype="solid"),
		axis.ticks.x = element_line(colour = "black"),
		axis.ticks.y = element_line(colour = "black"),
		axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
		axis.text.y = element_text(),
		axis.title.x = element_blank(),
		axis.title.y = element_blank()
	) 
gene <- c("EIF4E", "FOS", "MTA3", "ID1", "STAT1", "ASCL2", "GCM1")

{
    showdata <- vil.evt.reclu.expr.avg$RNA
    showdata <- showdata[rownames(showdata) %in% gene,]
    showdata <- showdata[match(gene, rownames(showdata)),]
    # showdata <- showdata[, colnames(showdata) != "2_AnCHM"]
}
pheatmap::pheatmap(showdata, scale = "row", cluster_rows = F, cluster_cols = F, angle_col = 45, 
                   color = colorRampPalette(c("navy", "white", "firebrick3"))(50),border_color = NA)


###### pseudotime ############
monocle.vil <- merge(vil.evt.reclu, list(vil.ctb.reclu, vil.stb.reclu))
monocle.vil <- merge(vil.evt.reclu, list(vil.ctb.reclu, vil.cctb))
monocle.vil <- merge(vil.ctb.reclu, list(vil.stb.reclu))
################## monocle3 ##########################################
library(monocle3)
data3 <- GetAssayData(monocle.vil, assay = 'integrated', slot = 'data')
fData3 <- data.frame(gene_short_name = rownames(data3), row.names = row.names(data3))
cds <- new_cell_data_set(data3, cell_metadata = monocle.vil@meta.data, gene_metadata = fData3)
cds <- preprocess_cds(cds, num_dim = 30, norm_method = "none")
# plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds, cores = 1, reduction_method = "UMAP", umap.n_neighbors = 30L, umap.min_dist = 0.2)
plot_cells(cds, color_cells_by="subtype")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(monocle.vil, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

cds <- cluster_cells(cds)
 
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "cluster")
cds <- learn_graph(cds, close_loop = T,
				   learn_graph_control = list(
				   	euclidean_distance_ratio = 1.2,
				   	# geodesic_distance_ratio = 0.4,
				   	# orthogonal_proj_tip = T,
				   	minimal_branch_len = 20
				   	),
				   use_partition = F)
plot_cells(cds, reduction_method="UMAP", color_cells_by="subtype", 
		   label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE,
		   group_label_size = 4, cell_size = 0.4) + 
	scale_color_manual(values = c("#90aacb", "#ff9e6f", "#00bfc4", "#f8766d", "#a3a500", "#00bf7d"))
	# scale_color_manual(values = c("#f8766d", "#00bfc4", "#fc8d62", "#66c2a5"))
plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
cds <- order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime",
		   # genes = sfeature('Kif5c'),
           label_cell_groups=F,
           label_leaves=F,
           label_branch_points=F,
		   cell_size = 0.4,
		   graph_label_size=.0)

plot_cells(cds, color_cells_by = "sample",
		   # genes = sfeature('Kif5c'),
           label_cell_groups=T,
           label_leaves=F,
           label_branch_points=F,
		   cell_size = 0.8,
		   graph_label_size=.0)

cds.pseudotime <- cds@principal_graph_aux@listData$UMAP$pseudotime

