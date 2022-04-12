#' An S4 class to extend Seurat class to include scHelper specific slots.
#'
#' @slot feature.meta A dataframe containing feature information. At minimum containing gene identifiers and associated gene names
#' @slot gene.modules A list of gene modules obtained from Antler
#' @name scHelper-class
#' @export
scHelper <- setClass(
  "scHelper",
  contains="Seurat",
  slots=c(feature.meta="data.frame",
          gene.modules = "list")
)

# Add show method for scHelper based on seurat default
setMethod("show", "scHelper", getMethod('show', 'Seurat'))

# Convert seurat object to scHelper class
setMethod("scHelper", "scHelper", function(object){
  as(object, "scHelper")
})