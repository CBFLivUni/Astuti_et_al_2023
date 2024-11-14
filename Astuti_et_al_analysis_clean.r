# Astuti et al scRNA-seq analysis ----
# Author: Dr Kim Clarke
# Date: July 2022

#' You may see small variations in the results if you do not use the same RNG seed. 
#' On rare occasions some systems may produce different results using the same seed.
#' Use the document outline panel in RStudio to find a specific section of the code.

# Libraries ----

library("tidyverse")
library("Seurat")
library("Matrix")
library("ggplot2")
library("patchwork")
library("dplyr")
library("tidyr")
library("ggrepel")
library("clusterProfiler")
library("RColorBrewer")
library("ggrepel")
library("org.Mm.eg.db")
library("extrafont")
library("DESeq2")
library("pheatmap")
library("GSVA")
library("parallel")
library("EnsDb.Hsapiens.v86")

# DATA FILTERING AND PREPARATION ----

rm(list = ls())

# Additional functions ----

source("code/r/additional_seurat_functions_1_loading_normalisation.r")

# Parameters ----

# remove mitochondrial and ribosomal subunit genes?

REMOVE_MITO_AND_RIBO = F

# basic filtering
MIN_GENES_LOAD = 200
MIN_CELLS_LOAD = 20

# clustering resolution
CLUSTER_RESOLUTION = 0.5

# number of principle components to use
PC_TO_USE_N = 20

# random seed where needed
RANDOM_SEED = 1234

# Sample identifiers ----

# batch 1

sample_ids = c("s1","s2","s3","s4")
sample_labels = c("Proximal","Distal","EarlyTumour","Naive")

# Data paths ----

# Cellranger batch 1
data_path = "data/cellranger_output_from_dhammond/mtx/"


# Load data and create Seurat objects ----

seurat_objs = create_seurat_objects_from_mtx(file_path = data_path, 
                                             sample_ids = sample_ids, 
                                             MIN_GENES_LOAD = MIN_GENES_LOAD, 
                                             MIN_CELLS_LOAD = MIN_CELLS_LOAD)


for(i in 1:length(sample_ids)){
  seurat_objs[[sample_ids[i]]]$orig.ident = sample_labels[i]
  seurat_objs[[sample_ids[i]]] = SetIdent(object = seurat_objs[[sample_ids[i]]], value = sample_labels[i])
}

names(seurat_objs) = sample_labels

# Adding metadata ----

seurat_objs = add_metadata_to_seurat_objects(obj_list = seurat_objs, mito_pattern = "^mt-")

# Assessing initial filtering ----

# before filtering, look at plots for the different filtering parameters

pdf(file = "combined.nfeat.ncount.pcMito.plots.raw.data.pdf", onefile = T)
for(id in sample_labels){
  p = VlnPlot(object = seurat_objs[[id]],features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), combine = F)
  p = lapply(p, function(x)return(x + theme(legend.position = "none")))
  
  # median nFeatures + 2*SD
  p[[1]] = p[[1]] + geom_hline(color = "red", yintercept = median(FetchData(seurat_objs[[id]], "nFeature_RNA")[,1]) + 2 * sd(FetchData(seurat_objs[[id]], "nFeature_RNA")[,1]))
  
  # median nFeatures - 1.8*SD
  p[[1]] = p[[1]] + geom_hline(color = "red", yintercept = median(FetchData(seurat_objs[[id]], "nFeature_RNA")[,1]) - 1.8 * sd(FetchData(seurat_objs[[id]], "nFeature_RNA")[,1]))
  
  # nFeatures = 200
  p[[1]] = p[[1]] + geom_hline(color = "blue", yintercept = 200)
  
  # median percent.mito + 2*SD
  p[[3]] = p[[3]] + geom_hline(color = "red", yintercept = median(FetchData(seurat_objs[[id]], "percent.mito")[,1]) + 2 * sd(FetchData(seurat_objs[[id]], "percent.mito")[,1]))
  
  # percent.mito = 0.1
  p[[3]] = p[[3]] + geom_hline(color = "blue", yintercept = 0.1)
  
  print(p[[1]] + p[[2]] + p[[3]])
  
  #print(p)
}
dev.off()

# Initial filtering ----

seurat_objs = lapply(seurat_objs, function(seurat_obj){
  
  # basic filtering
  MAX_MITO = 0.1
  MIN_GENES_FILTER = 200
  MAX_GENES_FILTER = 6000
  
  seurat_obj = seurat_obj[,FetchData(seurat_obj,"percent.mito") < MAX_MITO]
  seurat_obj = seurat_obj[,FetchData(seurat_obj,"nFeature_RNA") > MIN_GENES_FILTER & FetchData(seurat_obj,"nFeature_RNA") < MAX_GENES_FILTER]
  
  print(paste0("MAX_MITO: ", MAX_MITO))
  print(paste0("MIN_GENES_FILTER: ", MIN_GENES_FILTER))
  print(paste0("MAX_GENES_FILTER: ", MAX_GENES_FILTER))
  print(paste0("Dimensions of seurat object after filtering - ",FetchData(seurat_obj, "orig.ident")[1,]))
  print(dim(seurat_obj))
  
  return(seurat_obj)
  
})

# Manual Filtering ----

# how much of the counts from each cell are coming from ribosomal subunits

seurat_objs = lapply(seurat_objs, function(obj){
  
  feats = grep("^RPL[0-9XYAL]+$|^RPS[0-9XYAL]+$", 
               rownames(obj), 
               value = T, 
               ignore.case = T)
  obj[["ribosomal_pct"]] = PercentageFeatureSet(obj, features = feats)
  return(obj)
  
})

# we might also want to remove ribosomal and mitochondrial genes, as is common in the field

if(REMOVE_MITO_AND_RIBO){
  
  seurat_objs = lapply(seurat_objs, function(obj){
    
    ribo_genes = grep("^RPL[0-9XYAL]+$|^RPS[0-9XYAL]+$", rownames(obj), value = T, ignore.case = T)
    mito_genes = grep("^MT-", rownames(obj), value = T, ignore.case = T)
    
    to_remove = c(ribo_genes, mito_genes)
    
    print(paste0("removing ", length(to_remove), " genes out of ", nrow(obj)))
    print(paste0("lowest number of genes in a cell before filtering is ", min(obj$nFeature_RNA)))
    
    obj = obj[!rownames(obj) %in% to_remove,]
    
    print(paste0("lowest number of genes in a cell is now ", min(obj$nFeature_RNA)))
    
    return(obj)
    
  })
  
}

# Merge seurat objects ----

seurat_merge = merge(x = seurat_objs[[1]], y = seurat_objs[-c(1)])


# LogNormalize and Scaling ----

seurat_merge_lognorm = NormalizeData(seurat_merge, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_merge_lognorm = ScaleData(object = seurat_merge_lognorm)
seurat_merge_lognorm = FindVariableFeatures(object = seurat_merge_lognorm, selection.method = "vst", nfeatures = 2000)

# remove Adgre1 and Cd68 negative cells ----

seurat_merge_lognorm = PercentageFeatureSet(seurat_merge_lognorm, features = c("Adgre1"), col.name = "Adgre1")
seurat_merge_lognorm = PercentageFeatureSet(seurat_merge_lognorm, features = c("Cd68"), col.name = "Cd68")
seurat_merge_filt = seurat_merge_lognorm[,seurat_merge_lognorm$Adgre1 > 0 & seurat_merge_lognorm$Cd68 > 0]


PC_TO_USE_N = 30 # Try with 30 PCs in the filtered data
seurat_merge_filt <- RunPCA(object = seurat_merge_filt, verbose = F, npcs = PC_TO_USE_N, seed.use = RANDOM_SEED)
ElbowPlot(seurat_merge_filt, ndims = PC_TO_USE_N)
PC_TO_USE_N = 20 # Reduce to 20 based on ElbowPlot

seurat_merge_filt@misc$dims = 1:PC_TO_USE_N

seurat_merge_filt = RunUMAP(seurat_merge_filt, reduction = "pca", dims = seurat_merge_filt@misc$dims, seed.use = RANDOM_SEED)
seurat_merge_filt <- FindNeighbors(object = seurat_merge_filt, 
                                   dims = seurat_merge_filt@misc$dims)
seurat_merge_filt <- FindClusters(object = seurat_merge_filt, 
                                  resolution = 0.44, 
                                  random.seed = RANDOM_SEED)

# create metadata column to split clusters by sample

seurat_merge_filt[["sample_clusters"]] = paste(seurat_merge_filt$orig.ident, seurat_merge_filt$seurat_clusters)

# relevel the sample name factor 

seurat_merge_lognorm$orig.ident = factor(seurat_merge_lognorm$orig.ident, levels = c("Naive","EarlyTumour","Proximal","Distal"))
seurat_merge_filt$orig.ident = factor(seurat_merge_filt$orig.ident, levels = c("Naive","EarlyTumour","Proximal","Distal"))

# Adding sample information ----

# Adgre1 and Cd68 filtered data
condition = rep("naive", ncol(seurat_merge_filt))
condition[seurat_merge_filt$orig.ident %in% c("Proximal","Distal")] = "late"
condition[seurat_merge_filt$orig.ident %in% c("EarlyTumour")] = "early"
seurat_merge_filt[["condition"]] = condition

# Creating dataset without naive clusters ----

seurat_merge_filt_nonnaive = seurat_merge_filt[, !seurat_merge_filt$seurat_clusters %in% c(0, 5, 9)]
seurat_merge_filt_nonnaive$seurat_clusters = droplevels(seurat_merge_filt_nonnaive$seurat_clusters)

# DIFFERENTIAL EXPRESSION ----

# additional functions

source("code/r/getFunctionalProfile.r")
adj_pval_threshold = 0.01

# functional analysis outdir
outdir = "filtered_data_markers"

# Kupffer / MoM cell labels ----

kupffer_cell_clusters = c(0, 1, 4, 5, 6, 8, 9) 
mom_cell_clusters = c(2, 3, 7)

kupffer_or_mom = rep("undefined", ncol(seurat_merge_filt))
kupffer_or_mom[ seurat_merge_filt$seurat_clusters %in% mom_cell_clusters ] = "mom"
kupffer_or_mom[ seurat_merge_filt$seurat_clusters %in% kupffer_cell_clusters ] = "kupffer"
seurat_merge_filt[["kupffer_or_mom"]] = kupffer_or_mom



# Marker detection for each cluster vs the rest of the sample ----

outdir = "filtered_data_markers"

# All clusters

unique_clusters = as.numeric(levels(seurat_merge_filt$seurat_clusters))
adj_pval_threshold = 0.01
marker_results_list = list(all_clusters = list(), non_naive = list())

for(cl in unique_clusters){
  #for(cl in c(12)){
  
  print("Analysing cluster:")
  print(cl)
  
  m = FindMarkers(object = seurat_merge_filt, ident.1 = cl)
  m = m[m$p_val_adj < adj_pval_threshold,]
  marker_results_list$all_clusters[[paste0("cluster_", cl)]] = m
  
  # write table for IPA
  
  fh = paste0("marker.table.for.IPA.with_naive.cluster-",  cl, ".txt")
  
  write.table(m[,c("avg_logFC","p_val_adj")], file = file.path(outdir,fh), sep="\t", quote=F, col.names = NA)
  
}


# COLOUR PALETTES AND THEME ----


# clusters / heatmap column bars
all_clusters_pal = scales::hue_pal()(10)
names(all_clusters_pal) = 0:9

# naive, early mets, late mets
time_point_pal = c(Naive = "#117733",`Early Mets` = "#AA4499", `Late Mets` = "#0072B2")

# pmam and dmam
pmam_pal = c(Proximal = "#332288", Distal = "#CC6677")

# heatmap colour scale
# heatmap_pal = c()
heatmap_pal = c("#AAEEFF","#CC6677","#661100") # lightblue / maroon 
# heatmap_pal = c("#456DA2","#FAF6B6","#CA3D38") # blue yellow red
# heatmap_pal = c("black","purple","yellow")
  
# ggplot theme

general_theme = theme(text = element_text(size = 12, family = "ArialMT"))

umap_theme = theme(plot.title = element_text(hjust = 0.5), 
                   legend.position = "none", 
                   axis.text = element_blank(), 
                   axis.title = element_text(size = 12), 
                   text = element_text(size = 12, family = "ArialMT"),
                   axis.ticks = element_blank())



# RESULTS USING ALL CLUSTERS ----

# Cluster composition ----

df = data.frame(apply(t(table(seurat_merge_filt$seurat_clusters, seurat_merge_filt$condition)), 2, as.numeric))
# colnames(df) = paste0("Cluster", gsub("X","",colnames(df)))
df  = df %>% rowwise() %>% mutate(size = sum(across(everything()))) %>% ungroup()
df$sample = c("Early","Late","Naive")

df_pct = df %>% mutate(across(!contains(c("size", "sample")), function(x)x/size))
df_pct_long = pivot_longer(df_pct, cols = !contains(c("size", "sample")))
df_pct_long$sample = factor(df_pct_long$sample, levels = c("Naive","Early","Late"))
df_pct_long = df_pct_long %>% mutate(cluster = gsub("X", "", name)) %>% dplyr::select(-name)
df_pct_long$cluster = factor(df_pct_long$cluster, levels = rev(sort(unique(df_pct_long$cluster))))

ggplot(df_pct_long, aes(x = sample, y = value*100, fill = cluster, color = I("black"))) + 
  geom_bar(position = "stack", stat = "identity") + 
  general_theme +
  ylab("Composition (%)") +
  scale_fill_manual(values = all_clusters_pal) + 
  coord_flip() +
  guides(fill = guide_legend(reverse = T))


# Heatmap of degs ----

marker_top_frac = sapply(marker_results_list$all_clusters, function(x)rownames(x[order(x$avg_logFC, decreasing = T),])[1:20])
marker_top_frac = pivot_longer(data.frame(marker_top_frac, stringsAsFactors = F),cols = everything()) %>% arrange(name)
heatmap_order = tapply(marker_top_frac$name, marker_top_frac$value, function(x)paste(x, collapse = ";"))
heatmap_order = data.frame(gene = names(heatmap_order), cluster = heatmap_order, stringsAsFactors = F) %>% arrange(cluster)
heatmap_order$hm_splits = sapply(strsplit(heatmap_order$cluster, ";"), function(x)x[1])

yintercepts = cumsum(rev(table(heatmap_order$hm_splits))) + 0.5
xintercepts = cumsum(table(seurat_merge_filt$seurat_clusters)) + seq(8, (24*9)+8, length.out = 10) # manual offset to make the lines fit the Seurat heatmap

h1 = DoHeatmap(object = seurat_merge_filt, features = heatmap_order$gene, 
               group.by = "seurat_clusters", 
               group.colors = all_clusters_pal,
               group.bar.height = 0.075, 
               size = 4, 
               slot = "scale.data", 
               disp.max = 2, 
               disp.min = 0, 
               raster = F, combine = T, draw.lines = T) + 
  scale_fill_gradientn(colors = c(heatmap_pal[1], heatmap_pal[2], heatmap_pal[3])) + 
  theme(legend.position = "none", axis.text.y = element_blank(), text = element_text(size = 9)) +
  geom_hline(yintercept = yintercepts, color = "white", size = 0.7) +
  geom_vline(xintercept = xintercepts, color = "white", size = 0.7)

h1

# Heatmap of average expression of degs ----

marker_top_frac = sapply(marker_results_list$all_clusters, function(x)rownames(x[order(x$avg_logFC, decreasing = T),])[1:20])
marker_top_frac = pivot_longer(data.frame(marker_top_frac, stringsAsFactors = F),cols = everything()) %>% arrange(name)

feats = unique(marker_top_frac$value)

feats_avg = t(apply(seurat_merge_filt@assays$RNA@data[feats,], 1, function(x)tapply(x, seurat_merge_filt$seurat_clusters, mean)))
pheatmap_col_annotation = data.frame(cluster = factor(0:9), row.names = 0:9)
pheatmap_colours = list(cluster = all_clusters_pal)

hm = pheatmap(feats_avg, cluster_cols = F, scale = "row", cluster_rows = F, breaks = seq(-1.5, 1.5, length.out = 101),
         annotation_col = pheatmap_col_annotation, annotation_colors = pheatmap_colours, show_rownames = F)
#larger version with row labels
hm = pheatmap(feats_avg, cluster_cols = F, scale = "row", cluster_rows = F, breaks = seq(-1.5, 1.5, length.out = 101),
              annotation_col = pheatmap_col_annotation, annotation_colors = pheatmap_colours, fontsize_row = 8, 
              filename = "heatmap.pdf", height = 18)


# KC/MoM markers: heatmap of KC vs MoM clusters showing top DEGs, with cluster labels

df = data.frame(custom_marker_lists$KC_vs_MoM)
sel = df %>% arrange(desc(abs(avg_logFC))) %>% head(50) %>% arrange(desc(avg_logFC)) %>% rownames
df_sel = df %>% arrange(desc(abs(avg_logFC))) %>% head(50) %>% arrange(desc(avg_logFC))
seurat_merge_filt[["kupffer_or_mom_with_cluster"]] = paste(seurat_merge_filt$kupffer_or_mom, seurat_merge_filt$seurat_clusters)

pdf(file = "pdf/KC_vs_MoM_heatmap_top50byFC.pdf", width = 10, height = 8)
DoHeatmap(object = seurat_merge_filt, features = sel, group.by = c("kupffer_or_mom_with_cluster"), 
          cells = Cells(seurat_merge_filt)[seurat_merge_filt$kupffer_or_mom %in% c("kupffer", "mom")], 
          disp.min = -2, disp.max = 2, raster = FALSE) # +
  # scale_fill_gradientn(colors = c(heatmap_pal[1], heatmap_pal[2], heatmap_pal[3]))
dev.off()




# UMAP plots: clusters ----

DimPlot(object = seurat_merge_filt, reduction = "umap", label = T, label.size = 5) +
  umap_theme +
  ggtitle("Clusters")

# UMAP plots: timepoint ----

DimPlot(object = seurat_merge_filt, reduction = "umap", order = T, sizes.highlight = 0.05, pt.size = 0.05, 
        cells.highlight = colnames(seurat_merge_filt)[seurat_merge_filt$condition == "naive"], cols.highlight = time_point_pal[1]) + 
  umap_theme + ggtitle("Naive")

DimPlot(object = seurat_merge_filt, reduction = "umap", order = T, sizes.highlight = 0.05, pt.size = 0.05, 
        cells.highlight = colnames(seurat_merge_filt)[seurat_merge_filt$condition == "early"], cols.highlight = time_point_pal[2]) + 
  umap_theme + ggtitle("Early Mets")

DimPlot(object = seurat_merge_filt, reduction = "umap", order = T, sizes.highlight = 0.05, pt.size = 0.05, 
        cells.highlight = colnames(seurat_merge_filt)[seurat_merge_filt$condition == "late"], cols.highlight = time_point_pal[3]) + 
  umap_theme + ggtitle("Late Mets")

# UMAP plots: proximal / distal ----

DimPlot(object = seurat_merge_filt[, seurat_merge_filt$orig.ident %in% c("Proximal", "Distal")], 
        reduction = "umap", 
        group.by = "orig.ident",
        order = F, 
        pt.size = 0.1) +
  scale_colour_manual(values = c(pmam_pal[2], pmam_pal[1]))


# UMAP plots: gene expression ----

# Violin plots ----

# Vsig4, Clec4f, Timd4, Ccr2

feats = c("Vsig4", "Clec4f", "Timd4", "Ccr2")

plt = VlnPlot(object = seurat_merge_filt[,seurat_merge_filt$orig.ident %in% c("Proximal","Distal")], features = feats, pt.size = 0, group.by = "orig.ident", same.y.lims = T, combine = F)
plt_a = plt[[1]] + theme(legend.position = "left", axis.title.x = element_blank()) + 
  scale_color_manual(values = pmam_pal, aesthetics = "fill") + 
  general_theme  
plt_b = lapply(plt[2:length(plt)], function(p)p + theme(axis.line.y = element_blank(), 
                                                     axis.title.y = element_blank(), 
                                                     axis.text.y = element_blank(), 
                                                     axis.ticks.y = element_blank(),
                                                     axis.title.x = element_blank(),
                                                     legend.position = "none") +
                scale_color_manual(values = pmam_pal, aesthetics = "fill") +
                general_theme)
plt_a | plt_b



# ANALYSIS OF MoM CLUSTERS ONLY ----

# Filter the data and re-do the PCA and UMAP

seurat_merge_filt_mom = seurat_merge_filt[,seurat_merge_filt$kupffer_or_mom == "mom"]

seurat_merge_filt_mom = RunPCA(seurat_merge_filt_mom, seed.use = RANDOM_SEED)
# ElbowPlot(seurat_merge_filt_mom, ndims = 50) # 20 PCs is still OK
seurat_merge_filt_mom = RunUMAP(seurat_merge_filt_mom, dims = 1:20, seed.use = RANDOM_SEED)

seurat_merge_filt_mom = FindNeighbors(seurat_merge_filt_mom, dims = seurat_merge_filt_mom@misc$dims)
seurat_merge_filt_mom = FindClusters(object = seurat_merge_filt_mom, resolution = 0.44, random.seed = RANDOM_SEED)

# comparing new clusters with old clusters

g1 = Cells(seurat_merge_filt_mom)[seurat_merge_filt_mom$seurat_clusters == 2]
g2 = Cells(seurat_merge_filt)[seurat_merge_filt$sample_clusters == "EarlyTumour 2"]

length(intersect(g1, g2)) / length(g2)

# combine plots into an initial PDF to send to Yuliana ----

pdf("pdf/mom_reclustering_umap.pdf", width = 4, height = 4)

# UMAP plots

DimPlot(object = seurat_merge_filt_mom, reduction = "umap", label = T, label.size = 7) +
  umap_theme +
  ggtitle("Clusters")

dev.off()

# timepoint 

pdf("pdf/mom_reclustering_umap_samples.pdf", width = 9, height = 3)

u1 = DimPlot(object = seurat_merge_filt_mom, reduction = "umap", order = T, sizes.highlight = 0.05, pt.size = 0.05, 
        cells.highlight = colnames(seurat_merge_filt)[seurat_merge_filt$condition == "naive"], cols.highlight = time_point_pal[1]) + 
  umap_theme + ggtitle("Naive")

u2 = DimPlot(object = seurat_merge_filt_mom, reduction = "umap", order = T, sizes.highlight = 0.05, pt.size = 0.05, 
        cells.highlight = colnames(seurat_merge_filt)[seurat_merge_filt$condition == "early"], cols.highlight = time_point_pal[2]) + 
  umap_theme + ggtitle("Early Mets")

u3 = DimPlot(object = seurat_merge_filt_mom, reduction = "umap", order = T, sizes.highlight = 0.05, pt.size = 0.05, 
        cells.highlight = colnames(seurat_merge_filt)[seurat_merge_filt$condition == "late"], cols.highlight = time_point_pal[3]) + 
  umap_theme + ggtitle("Late Mets")

u1 + u2 + u3

dev.off()

# DEGs

unique_clusters = as.numeric(levels(seurat_merge_filt_mom$seurat_clusters))
adj_pval_threshold = 0.01
marker_results_list$mom = list()

for(cl in unique_clusters){
  
  print("Analysing cluster:")
  print(cl)
  
  m = FindMarkers(object = seurat_merge_filt_mom, ident.1 = cl)
  m = m[m$p_val_adj < adj_pval_threshold,]
  marker_results_list$mom[[paste0("cluster_", cl)]] = m
  
  # write table for IPA
  
  fh = paste0("marker.table.MoM.reclustering.cluster-",  cl, ".txt")
  
  write.table(m[,c("avg_logFC","p_val_adj")], file = file.path(fh), sep="\t", quote=F, col.names = NA)
  
}

marker_top_frac = sapply(marker_results_list$mom, function(x)rownames(x[order(x$avg_logFC, decreasing = T),])[1:20])
marker_top_frac = pivot_longer(data.frame(marker_top_frac, stringsAsFactors = F),cols = everything()) %>% arrange(name)
heatmap_order = tapply(marker_top_frac$name, marker_top_frac$value, function(x)paste(x, collapse = ";"))
heatmap_order = data.frame(gene = names(heatmap_order), cluster = heatmap_order, stringsAsFactors = F) %>% arrange(cluster)
heatmap_order$hm_splits = sapply(strsplit(heatmap_order$cluster, ";"), function(x)x[1])

# yintercepts = cumsum(rev(table(heatmap_order$hm_splits))) + 0.5
# xintercepts = cumsum(table(seurat_merge_filt$seurat_clusters)) + seq(8, (24*9)+8, length.out = 10) # manual offset to make the lines fit the Seurat heatmap

pdf("pdf/mom_reclustering_results_top_degs_heatmap.pdf", height = 14)

DoHeatmap(object = seurat_merge_filt_mom, features = heatmap_order$gene, 
               group.by = "seurat_clusters", 
               # group.colors = all_clusters_pal,
               group.bar.height = 0.075, 
               size = 4, 
               slot = "scale.data", 
               disp.max = 2, 
               disp.min = -2, 
               raster = F, combine = T, draw.lines = T) + 
  # scale_fill_gradientn(colors = c(heatmap_pal[1], heatmap_pal[2], heatmap_pal[3])) + 
  theme(legend.position = "none", text = element_text(size = 10)) # +
  # geom_hline(yintercept = yintercepts, color = "white", size = 0.7) +
  # geom_vline(xintercept = xintercepts, color = "white", size = 0.7)

dev.off()






# ~~~~~~~~ ----
# Human bulk RNA-seq analysis ----

# GSVA analysis of bulk liver mets RNA-seq with custom gene sets ----
# Adapted from code provided by Peter Bailey, University of Glasgow

# Load and process data ----

mets <- read.delim("~/schmid_single_cell/glasgow_mets_rnaseq/merged_gene_counts.txt", sep="\t", row.names=1, stringsAsFactors = FALSE)
mets = mets[,-c(grep("Takara", colnames(mets)))]

groups <- sapply(colnames(mets), function(i) strsplit(i, "_")[[1]][1])

tabFinal <- data.frame(ID=colnames(mets), groups=groups, row.names=colnames(mets))

MJ.raw.counts.rsub.filt <- mets[ rowSums(edgeR::cpm(mets) > 1) >= 4, ]

# Produce a DESeq dataset
MJ.raw.rsub.se <- DESeqDataSetFromMatrix(countData = MJ.raw.counts.rsub.filt, colData = tabFinal, design = ~ groups)

# Perform rlog transformation of the raw count data
MJ.raw.rsub.rlog <- rlog(MJ.raw.rsub.se)
dat.rlog <- assay(MJ.raw.rsub.rlog)

geneNames <- select(EnsDb.Hsapiens.v86, key=rownames(dat.rlog), columns=c("ENTREZID", "SYMBOL"), keytype="GENEID")
geneNames <- geneNames[!duplicated(geneNames$GENEID), ]
geneNames <- geneNames[!is.na(geneNames$SYMBOL), ]
rownames(geneNames) <- geneNames$GENEID

dat.rlog.genes <- dat.rlog
rownames(dat.rlog.genes) <- geneNames[rownames(dat.rlog.genes), "SYMBOL"]

# Subtype enrichment ----

activated_stroma_genes <- c("ZNF469", "VCAN", "THBS2", "SULF1", "SPARC", "SFRP2", "POSTN", "MMP11", "LUM", "ITGA11", "INHBA", "GREM1", "FNDC1", "FN1", "FAP", "CTCHR1", "COMP", "COL5A2", "COL5A1", "COL3A1", "COL1A2", "COL1A1", "COL11A1", "COL10A1", "CDH11")
normal_stroma_genes <- c("VIT", "SYNM", "SCRG1", "RSPO3", "RERGL", "RBPMS2", "PTX3", "PLP1", "OGN", "MYH11", "MEOX2", "LPHN3", "LMOD1", "IGF1", "ID4", "GPM6B", "FABP4", "DES", "CDH19", "ANGPTL7", "ADAMTS1", "ACTG2", "ABCA8")
basal_like <- c("VGLL1", "UCA1", "S100A2", "LY6D", "SPRR3", "SPRR1B", "LEMD1", "KRT15", "CTSL2", "DHRS9", "AREG", "CST6", "SERPINB3", "KRT6C", "KRT6A", "SERPINB4", "FAM83A", "SCEL", "FGFBP1", "KRT7", "KRT17", "GPR87", "TNS4", "SLC2A1", "ANXA8L2")
Classical <- c("BTNL8", "FAM3D", "ATAD4", "AGR3", "CTSE", "LOC400573", "LYZ", "TFF2", "TFF1", "ANXA10", "LGALS4", "PLAG2G10", "CEACAM6", "VSIG2", "TSPAN8", "ST6GALNAC1", "AGR2", "TFF3", "CYP3A7", "MYO1A", "CLRN3", "KRT20", "CDH17", "SPINK4", "REG4")

# immune cell markers from Rooney et al

rooney_markers_list = list(
  b_cells = c("CD79B","BTLA","FCRL3","BANK1","CD79A","BLK","RALGPS2","FCRL1","HVCN1","BACH2"),
  cd4_t_cells = c("FOXP3","C15orf53","IL5","CTLA4","IL32","GPR15","IL4"),
  cd8_t_cells = c("CD8A"),
  macrophages = c("FUCA1","MMP9","LGMN","HS3ST2","TM4SF19","CLEC5A","GPNMB","C11orf45","CD68","CYBB"),
  neutrophils = c("KDM6B","HSD17B11","EVI2B","MNDA","MEGF9","SELL","NLRP12","PADI4","TRANK1","VNN3"),
  nk_cells = c("KLRF1","KLRC1"),
  pdc = c("LILRA4","CLEC4C","PLD4","PHEX","IL3RA","PTCRA","IRF8","IRF7","GZMB","CXCR3"),
  co_stimulation_apc = c("ICOSLG","CD70","TNFSF14","CD40","TNFSF9","TNFSF4","TNFSF15","TNFSF18","TNFSF8","SLAMF1","CD58"),
  co_stimulation_t_cell = c("ICOS","CD28","CD27","TNFSF14","CD40LG","TNFRSF9","TNFRSF4","TNFRSF25","TNFRSF18","TNFRSF8","SLAMF1","CD2","CD226"),
  co_inhibition_apc = c("PDCD1LG2","CD274","C10orf54","LGALS9","PVRL3"),
  co_inhibition_t_cell = c("LAG3","CTLA4","CD274","CD160","BTLA","C10orf54","LAIR1","HAVCR2","CD244","TIGIT"),
  cytolytic_activity = c("GZMA","PRF1"),
  type_i_ifn_response = c("MX1", "TNFSF10", "RSAD2", "IFIT1", "IFIT3", "IFIT2", "IRF7", "DDX4", "MX2", "ISG20"),
  type_ii_ifn_response = c("GPR146", "SELP", "AHR")
)

# heatmap ----

heatmap_order = c(`Macrophage` = "macrophages",
                  `Neutrophils` = "neutrophils",
                  `CD8 T cells` = "cd8_t_cells",
                  `CD4 T cells` = "cd4_t_cells",
                  `NK cells` = "nk_cells",
                  `B cells` = "b_cells",
                  `Co-stimulation, APC` = "co_stimulation_apc",
                  `Co-inhibitory, APC` = "co_inhibition_apc",
                  `Co-stimulation, T cells` = "co_stimulation_t_cell",
                  `Co-inhibitory, T cells` = "co_inhibition_t_cell",
                  `Cytolytic activity` = "cytolytic_activity",
                  `Type I IFN response` = "type_i_ifn_response",
                  `Type II IFN response` = "type_ii_ifn_response")

pdf(file = "pdf/Immune_signature_scores_v2.pdf", width = 10, height = 5)
pheatmap(t(scores_per_sample[heatmap_order,]), 
         scale = "none", 
         breaks = seq(-0.4, 0.4, length.out = 101), 
         cluster_rows = F, cluster_cols = F)
dev.off()


# Cluster 1, 4, 2, 3 signatures in patient 1-5 bulk RNAseq data ---- 
# heatmap of gene expression
# (and dot plot if can show 2 values: score and p-value)

cluster_1234_sigs = list(cluster1 = marker_results_list$non_naive$cluster_1_nonnaive %>% dplyr::filter(avg_logFC > 0, p_val_adj < 0.01) %>% rownames() %>% toupper(),
     cluster2 = marker_results_list$non_naive$cluster_2_nonnaive %>% dplyr::filter(avg_logFC > 0, p_val_adj < 0.01) %>% rownames() %>% toupper(),
     cluster3 = marker_results_list$non_naive$cluster_3_nonnaive %>% dplyr::filter(avg_logFC > 0, p_val_adj < 0.01) %>% rownames() %>% toupper(),
     cluster4 = marker_results_list$non_naive$cluster_4_nonnaive %>% dplyr::filter(avg_logFC > 0, p_val_adj < 0.01) %>% rownames() %>% toupper())

cluster_1234_sigs_table = do.call(rbind, lapply(names(cluster_1234_sigs), function(x)
  data.frame(symbol = toupper(cluster_1234_sigs[[x]])[toupper(cluster_1234_sigs[[x]]) %in% rownames(dat.rlog.genes)],
             cluster = x)
))

cluster_1234_sigs_table = cluster_1234_sigs_table[!duplicated(cluster_1234_sigs_table$symbol),]

mat = dat.rlog.genes[match(cluster_1234_sigs_table$symbol, rownames(dat.rlog.genes)),]
matm = t(apply(mat, 1, function(x)tapply(x, tabFinal$groups, mean)))


pheatmap(matm, cluster_rows = F, cluster_cols = F, scale = "column",
         annotation_row = data.frame(cluster = cluster_1234_sigs_table$cluster, row.names = cluster_1234_sigs_table$symbol), 
         breaks = seq(-2, 2, length.out = 101))

# GSVA

scores  <- GSVA::gsva(expr = dat.rlog.genes, gset.idx.list=cluster_1234_sigs, method="ssgsea", ssgsea.norm = T, 
                      parallel.sz = 1)

df <- scores %>% t() %>% as.data.frame() %>% mutate(description=MJ.raw.rsub.se$groups)
scoresMelt <- reshape2::melt(df)

scoresMelt <- scoresMelt %>% mutate(description=  factor(description, levels = c("WholeTissue.01", "WholeTissue.02",  "WholeTissue.03", "WholeTissue.04", "WholeTissue.05"))) %>% arrange(description)

ord = names(sort(tapply(scoresMelt$value, scoresMelt$variable, median), decreasing = T))

scoresMelt$variable = factor(scoresMelt$variable, levels = ord)

scores_per_sample = sapply(by(t(scores), tabFinal$groups, function(x)colMeans(x)), function(x)x)

# dotplot using geom_point ----

scores_melt = pivot_longer(data = data.frame(signature = rownames(scores_per_sample), scores_per_sample), cols = starts_with("Whole"))
scores_melt = scores_melt %>% mutate(abs_value = abs(value))
names(scores_melt) = c("Signature","Sample","Score","Abs(Score)")
scores_melt$Signature = factor(scores_melt$Signature, levels = c("cluster1", "cluster4", "cluster2", "cluster3"))
scores_melt[["size"]] = ifelse(scores_melt$Score > 2, 2, scores_melt$Score)


dp = ggplot(scores_melt, aes(x = Signature, y = Sample, size = size)) + 
  geom_point(aes(fill = Score), color = "black", pch = 21) + 
  scale_size_continuous(range = c(1,15), limits = c(1.4,2.6)) +
  scale_fill_distiller(palette = "Spectral", limits = c(1, 2.5), oob = scales::squish) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0)) +
  scale_x_discrete(position = "top")
dp

# heatmap ----

scores_hm = data.frame(signature = rownames(scores_per_sample), scores_per_sample) %>% 
  select(starts_with("WholeTissue"))

scores_hm = t(scores_hm)[,c(1,4,2,3)]

pheatmap(scores_hm, cluster_rows = F, cluster_cols = F,
         filename = "pdf/cluster_1234_sig_scores_in_human_heatmap.pdf",
         height = 2, width = 3)








































