############# Functions to help create a list of Seurat objects from a set of single-cell RNAseq cellranger files and perform filtering and normalisation ############

#' Inputs:
#' 
#' file_path = directory containing filtered CellRanger expression matrices in mtx format, feature and cell barcode files in tsv format
#' sample_ids = unique identifier appended to the beginning of each of the files in file.path. Format should be "[sample_id]_[barcodes|features/genes|matrix].[tsv.gz|mtx.gz]"
#' whitelist = list of cell barcodes to keep when filtering a seurat object
#' filter parameters, see below


#' Filters used when loading data:
#' 1) MIN_GENES_LOAD = minimum number of genes detected per cell (used when loading the data)
#' 2) MIN_CELLS_LOAD = minimum number of cells that express a gene (used when loading the data)

create_seurat_objects_from_mtx = function(file_path, sample_ids, MIN_GENES_LOAD, MIN_CELLS_LOAD){
  
  seurat_objs = lapply(sample_ids, function(id){
    
    print(paste0("Creating Seurat Object for:", id))
    
    barcode.path <- paste0(file_path, id, "_barcodes.tsv.gz")
    features.path <- paste0(file_path, id, "_features.tsv.gz")
    matrix.path <- paste0(file_path, id, "_matrix.mtx.gz")
    
    mat <- readMM(file = matrix.path)
    
    feature.names = read.delim(features.path, 
                                 header = FALSE,
                                 stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)
    
    colnames(feature.names) = c("ENSEMBL","SYMBOL","TYPE")
    feature.names$SYMBOL = make.unique(feature.names$SYMBOL)
    rownames(feature.names) = feature.names$SYMBOL
    
    colnames(mat) = barcode.names$V1
    rownames(mat) = feature.names$SYMBOL
    
    seurat_obj = CreateSeuratObject(
      counts = mat,
      min.cells = MIN_CELLS_LOAD ,
      min.features = MIN_GENES_LOAD,
      project = paste0("10X_",id)
    )
    
    seurat_obj@assays$RNA = AddMetaData(object = seurat_obj[["RNA"]], metadata = feature.names)
    
    return(seurat_obj)
    
  })
  
  names(seurat_objs) = sample_ids

  return(seurat_objs)
  
}



# filtering cells in a seurat object using a whitelist of cell barcodes

filter_seurat_cb_whitelist = function(obj_list, whitelist_list){
  
  for(id in 1:length(obj_list)){
    
    whitelist = whitelist_list[[id]]
    suffix = unique(gsub("[ATCG]","",colnames(obj_list[[id]])))
    if(length(suffix) == 1){
      whitelist = paste(whitelist, suffix, sep = "")
      
      sel = colnames(obj_list[[id]]) %in% whitelist
      print(paste0("Object ", id, " cells to keep: ", length(which(sel)), " / ", length(whitelist)))
      
      obj_list[[id]] = obj_list[[id]][,colnames(obj_list[[id]]) %in% whitelist]
    } else {
      print("more than 1 suffix for some reason")
    }
    
  }
  
  return(obj_list)
  
}
  

# adding standard metadata metrics to seurat objects

# metadata:
# percent.mito = percentage of counts coming from genes starting with "MT" i.e. mitochondrially encoded genes

add_metadata_to_seurat_objects = function(obj_list, mito_pattern = "^MT-"){

  print("Adding percent.mito to metadata")
  
  seurat_objs = lapply(obj_list, function(seurat_obj){
    
    symbols = seurat_obj@assays$RNA@meta.features$SYMBOL
    mito.genes =  grep(pattern = mito_pattern,x = symbols)
    
    percent.mito = Matrix::colSums(GetAssayData(object=seurat_obj,slot="counts",assay="RNA")[mito.genes,]) / Matrix::colSums(GetAssayData(object=seurat_obj,slot="counts",assay="RNA"))
    
    # AddMetaData adds columns to object@meta.data, a great place to stash QC stats
    seurat_obj = AddMetaData(object = seurat_obj, metadata = percent.mito, col.name="percent.mito")
    
    return(seurat_obj)
    
  })
  
}
  

sctransform_vignette = function(objs, outdir = "sctransform_vignette"){
  
  if(!dir.exists(outdir)){ dir.create(outdir) }
  
  for(s in names(objs)){
    
    pdf(file = paste0(outdir, "/sctransform.vignette.",s,".pdf"), width = 12)
    
    vignette_data = as.matrix(objs[[s]]@assays$RNA@counts)
    
    gene_attr <- data.frame(mean = rowMeans(vignette_data), 
                            detection_rate = rowMeans(vignette_data > 0),
                            var = apply(vignette_data, 1, var))
    gene_attr$log_mean <- log10(gene_attr$mean)
    gene_attr$log_var <- log10(gene_attr$var)
    rownames(gene_attr) <- rownames(vignette_data)
    cell_attr <- data.frame(n_umi = colSums(vignette_data),
                            n_gene = colSums(vignette_data > 0))
    rownames(cell_attr) <- colnames(vignette_data)
    
    # Plot: log mean expression vs log mean variance
    print("log mean vs log var")
    
    print(ggplot(gene_attr, aes(log_mean, log_var)) + 
            geom_point(alpha=0.3, shape=16) + 
            geom_density_2d(size = 0.3) +
            geom_abline(intercept = 0, slope = 1, color='red') +
            ggtitle("Log mean expression vs log variance") +
            theme_bw())
    
    
    # Plot: log mean expression vs detection rate
    print("detection rate")
    
    # add the expected detection rate under Poisson model
    x = seq(from = -3, to = 2, length.out = 1000)
    poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
    
    print(ggplot(gene_attr, aes(log_mean, detection_rate)) + 
            geom_point(alpha=0.3, shape=16) + 
            geom_line(data=poisson_model, color='red') +
            theme_gray(base_size = 8) +
            ggtitle("Log mean expression vs detection rate") +
            theme_bw())
    
    # Plot: nUMI vs detection rate relationship
    
    print("UMIs vs no of genes")
    
    print(ggplot(cell_attr, aes(n_umi, n_gene)) + 
            geom_point(alpha=0.3, shape=16) + 
            geom_density_2d(size = 0.3) +
            ggtitle("nUMIs vs nGenes") +
            theme_bw())
    
    # VST
    
    print("VST")
    
    set.seed(44)
    vst_out <- sctransform::vst(vignette_data, latent_var = c('log_umi'), return_gene_attr = TRUE, return_cell_attr = TRUE, show_progress = FALSE)
    
    print(sctransform::plot_model_pars(vst_out) + ggtitle("vst model parameters"))
    
    
    plot_model_genes = names(head(sort(rowMeans(objs[[s]]@assays$RNA@counts), decreasing = T))[1:5])
    plot_model_genes_hvgs = objs[[s]]@assays$SCT@var.features[1:5]
    
    sctransform::plot_model(vst_out, vignette_data, plot_model_genes, plot_residual = TRUE, gg_cmds = ggtitle('Highest expressed genes'))
    
    sctransform::plot_model(vst_out, vignette_data, plot_model_genes_hvgs, plot_residual = TRUE, gg_cmds = ggtitle('Highly variable genes'))
    
    print(ggplot(vst_out$gene_attr, aes(residual_mean)) + 
            geom_histogram(binwidth=0.01) +
            ggtitle("Distribution of residuals") +
            theme_bw())
    
    print(ggplot(vst_out$gene_attr, aes(log10(gmean), residual_variance)) + 
            geom_point(alpha=0.3, shape=16) +
            geom_density_2d(size = 0.3) + 
            ggtitle("Residual variance vs log mean expression") + 
            theme_bw())
    
    # ranked residual variance plot
    
    rv = log10(objs[[s]]@assays$SCT@meta.features$sct.residual_variance)
    varfeat = objs[[s]]@assays$SCT@meta.features$sct.variable
    ord = order(rv)
    plot(rv[ord], col = c("black","red")[as.numeric(varfeat[ord])+1], pch = 19, cex = 0.5, main = "Sorted residual variance", ylab = "Residual variance")
    legend(x = "topleft", legend = c("Var feat", "None var feat"), col = c("red","black"), pch = 19)
    
    dev.off()
    
  }
  
}


