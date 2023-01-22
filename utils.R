## calculate number of PCA components to explain certain variance
npc <- function(seuratobj, percent) {
    stdev = seuratobj@reductions$pca@stdev
    var = stdev ^ 2
    total = sum(var)
    i = 0
    while (i < length(var)) {
        i = i + 1
        pct = sum(var[1:i]) / total
        if (pct >= percent) {
            return (i)
        }
    }
    return(length(var))
}
## calculate variance explained percent of given number of PCA components
pc.explain <- function(seuratobj, npcs){
    stdev = seuratobj@reductions$pca@stdev
    var = stdev ^ 2
    total = sum(var)
    return(sum(var[1:npcs]) / total)
}



prop <- function(seuratdata, split.by, prop.by, stat = "percent", plot = T, return.plot = F,
				 orient = "v", plot.title = NULL, plot.legend.title = NULL, show.axis.title = F,
				 plot.colors = NULL) {
	prop.bar <- as.data.frame(table(seuratdata@meta.data[, split.by], seuratdata@meta.data[, prop.by])) %>% 
    	group_by(Var1) %>% 
    	group_map(function(df, key) {
    		if (stat == "percent") setNames(df$Freq / sum(df$Freq), df$Var2) 
    		else if (stat == "count") setNames(df$Freq, df$Var2) 
    	}) %>% as.data.frame(col.names = rownames(table(seuratdata@meta.data[, split.by], seuratdata@meta.data[, prop.by])), check.names = F)

	if (plot) {
		colvar <- c("value", "variable")
		x.title <- NULL
		if (orient == "v" & show.axis.title) x.title = split.by
		if (orient == "h") colvar <- rev(colvar)
		if (is.null(plot.legend.title)) plot.legend.title = prop.by
		
		prop.melt <- prop.bar
		prop.melt[, plot.legend.title] <- factor(rownames(prop.bar), levels = rownames(prop.bar))

		p <- ggplot(reshape2::melt(prop.melt), aes_string(y = colvar[1], x = colvar[2], fill = plot.legend.title)) + 
    		labs(y = NULL, x = x.title, title = plot.title) +
    		geom_bar(stat = "identity", color = "white", width = 0.95, position = position_stack(reverse = TRUE)) + 
			# scale_x_discrete(limits = levels(plot$celltype.1)) +
			theme(
				plot.margin = unit(c(0, 0.5, 0, 0.5), "cm"),
				panel.border = element_blank(),
				panel.background = element_blank(),
				axis.line.x = element_blank(),
				axis.ticks = element_blank(),
				axis.text.y = element_text(),
				axis.text.x = element_text(angle = 45, hjust =1, vjust = 1),
				axis.title.y = element_blank(),
				axis.title.x = element_text(size = 11),
				legend.position= "top"
			)
		if (orient == "v") p <- p + scale_y_continuous(expand = c(0,0))
		else if (orient == "h") p <- p + scale_x_continuous(expand = c(0,0))
		if (!is.null(plot.colors)) p <- p + scale_fill_manual(values = plot.colors)
		#+ guides(fill=guide_legend(nrow=1,bycol=F))
		show(p)
		if (return.plot) return(p)
	}
	return (prop.bar)
}

# DRTest <- function(count, property) {
#     count.DR <- DirichletReg::DR_data(count)
#     count$property <- property
#     fit <- DirichletReg::DirichReg(count.DR ~ property, count)
#     u <- summary(fit)
#     print(u$coef.mat)
#     pvals <- u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4] 
#     v <- names(pvals)
#     pvals <- matrix(pvals, ncol=length(u$varnames))
#     rownames(pvals) <- gsub('condition', '', v[1:nrow(pvals)])
#     colnames(pvals) <- u$varnames
#     fit$pvals <- pvals
#     fit
# }



getGOfromGene <- function(gene) {
    tmp <- bitr(gene, fromType = "SYMBOL",
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
}

GOgene <- function(id) {
	if (!startsWith(id, "GO:")) {
		id = paste0("GO", id)
	}
	unique(AnnotationDbi::select(org.Hs.eg.db, keys = id, 
							 columns = c("ENTREZID","SYMBOL","ONTOLOGY"), keytype="GO")$SYMBOL)
}
