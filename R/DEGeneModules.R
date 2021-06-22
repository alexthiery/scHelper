#' Filter gene modules based on differential expression (seurat::FindAllMarkers)
#'
#' @param seurat_object Seurat object
#' @param gene_modules named list of gene modules
#' @param logfc logFC cutoff used for differential expression
#' @param pval adjusted pvalue
#' @param selected_gene_proportion proportion of genes in gene module which are required to pass differential expression test
#' @return gene module list
#' @export
DEGeneModules <- function(seurat_data, gene_modules, logfc = 0.25, pval = 0.001, selected_gene_proportion = 0.5, active_ident = NULL){
  if(!is.null(active_ident)){Idents(object = seurat_data) <- active_ident}
  DE_genes <- FindAllMarkers(seurat_data, only.pos = T, logfc.threshold = logfc) %>% filter(p_val_adj < pval)
  gms <- SubsetGeneModules(gene_modules, selected_genes = DE_genes$gene, keep_mod_ID = T, selected_gene_proportion = selected_gene_proportion)
  return(gms)
}