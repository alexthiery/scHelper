#' Generate individual FeaturePlots for given list of genes
#' 
#' @param seurat_object seurat object
#' @param gene_list array of genes to plot
#' @param plot_path path to plot FeaturePlots to
#' @param zip_files boolean for whether or not to zip plots
#' @export
UMAPGeneList <- function(seurat_object, gene_list, plot_path, zip_files = TRUE){
  dir.create(plot_path, recursive = TRUE)
  
  missing_genes = c()
  for(g in gene_list){
    if(!g %in% rownames(seurat_object)){
      missing_genes = c(g, missing_genes)
      next
    }else{
      print(g)
      pdf(paste0(plot_path, g, "_UMAP.pdf"), height = 5, width = 10)
      plot(gridExtra::grid.arrange(grobs = c(list(DimPlot(seurat_object) +
                                                    ggtitle("Seurat clusters") +
                                                    theme(plot.title = element_text(hjust = 0.5))),
                                             list(FeaturePlot(seurat_object, g))), ncol = 2))
      dev.off()
    }
  }
  cat('\nFollowing genes not expressed in dataset:', missing_genes, '\n\n')

  if(zip_files == TRUE){
    system(paste0("zip -rj ", dirname(plot_path), "/", basename(plot_path), ".zip ", plot_path))
    unlink(plot_path, recursive=TRUE, force=TRUE)
  }
}