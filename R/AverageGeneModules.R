#' AverageGeneModules
#' 
#' Calculate average expression for groups of genes across Seurat clusters and add to seurat_obj@meta.data.
#' 
#' @param seurat_obj Seurat object
#' @param gene_list Named list of genes passed to seurat::AddModuleScore
#' @return seurat_obj
#' @export
AverageGeneModules <- function(seurat_obj, gene_list){
  if(is.null(names(gene_list))){stop('Gene list must contain named arrays')}
  
  # Override previous meta.data columns if present
  if(any(colnames(seurat_obj@meta.data) %in% names(gene_list))){
    seurat_obj@meta.data[, names(gene_list)] <- NULL
  }
  
  seurat_obj <- AddModuleScore(seurat_obj, features = gene_list, name = 'genes')
  colnames(seurat_obj@meta.data) <- replace(colnames(seurat_obj@meta.data), grepl('genes', colnames(seurat_obj@meta.data)), names(gene_list))
  return(seurat_obj)
}

