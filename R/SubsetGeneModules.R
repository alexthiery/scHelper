#' Subset Antler gene modules based on selected gene list provided
#' 
#' @param gms List of gene modules obtained from Antler
#' @param selected_genes array of genes used to filter gene module list
#' @param keep_mod_ID boolean value for whether to keep original gene module ID
#' @param selected_gene_proportion ratio of selected genes which must be present in order for the module to be kept
#' @export
SubsetGeneModules <- function(gms, selected_genes, keep_mod_ID = FALSE, selected_gene_proportion = 0) {
  outlist <- list()
  if(is.null(names(gms))){
    names(gms) <- paste0("GM: ", 1:length(gms))
  } else {}
  
  # this filters gene modules for which the percentage of differentially expressed genes passes a threshold test
  gms = gms[unlist(lapply(gms, function(x) sum(x %in% selected_genes) >= round(length(x)*selected_gene_proportion)))]
  
  for(mod in names(gms)){
    if(any(selected_genes %in% gms[[mod]])){
      gene.match <- stringr::str_c(selected_genes[which(selected_genes %in% gms[[mod]])], collapse = "; ")
      if(keep_mod_ID == T){
        outlist[[mod]] <- gms[[mod]]
      } else {
        outlist[[gene.match]] <- gms[[mod]]
      }
    }
    else{next}
  }
  return(outlist)
}





