#' PCA level comparison
#' 
#' Plot UMAPs using different numbers of principle components.
#' 
#' @param data Seurat object or list of Seurat objects
#' @param PCA_levels array of four max principle components to compare
#' @param cluster_res resolution used for KNN clustering
#' @export
PCALevelComparison <- function(data, PCA_levels = c(10, 15, 20, 30), cluster_res = 0.5){
  if(class(data) == "list") {
    data.stage.PCA <- list()
    for (stage in names(data)){
      for (val in PCA_levels){
        name<-paste0(stage, " PCA ", val)
        data.stage.PCA[[name]] <- FindNeighbors(data[[stage]], dims = 1:val, verbose = FALSE)
        data.stage.PCA[[name]] <- FindClusters(data.stage.PCA[[name]], resolution = cluster_res, verbose = FALSE)
        data.stage.PCA[[name]] <- RunUMAP(object = data.stage.PCA[[name]], dims = 1:val, verbose = FALSE)
      }
    }
    plots <- lapply(seq_along(1:length(data.stage.PCA)), function(dat){
      DimPlot(data.stage.PCA[[dat]]) + ggtitle(paste(names(data.stage.PCA[dat])))
    })
    return(gridExtra::grid.arrange(grobs = plots))
  } else {
    data.PCA <- list()
    for (val in PCA_levels){
      name<-paste0(" PCA ", val)
      data.PCA[[name]] <- FindNeighbors(data, dims = 1:val, verbose = FALSE)
      data.PCA[[name]] <- FindClusters(data.PCA[[name]], resolution = cluster_res, verbose = FALSE)
      data.PCA[[name]] <- RunUMAP(object = data.PCA[[name]], dims = 1:val, verbose = FALSE)
    }
    plots <- lapply(seq_along(1:length(data.PCA)), function(dat){
      DimPlot(data.PCA[[dat]]) + ggtitle(paste(names(data.PCA[dat])))
    })
    return(gridExtra::grid.arrange(grobs = plots))
  }
}
