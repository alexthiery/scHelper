#' Exports antler peak modules
#' 
#' This function writes out the peak modules calculated by antler to a text file
#'
#' @param antler_object antler object with the peak modules
#' @param publish_dir directory to write the peak modules to
#' @param names_list names of the peak modules to write out, as determined in antler object. Defaults to unbiasedPMs
#' @return .txt files with the peak modules
#' @export
ExportAntlerModules <- function(antler_object, publish_dir, names_list = "unbiasedPMs"){
  for(gm_list in names_list){
    mods = antler_object$gene_modules$lists[[gm_list]]$content
    for (i in seq(length(mods))) {
      modname = base::names(mods)[i]
      if (is.null(modname)) {
        modname = paste0("PM: ", i)
      }
      write(paste0(modname, "; ", paste0(mods[[i]], collapse = ", ")), file = paste0(publish_dir, '/', gm_list, '.txt'), append = TRUE)
    }
  }
}
