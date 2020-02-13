#' GetGatkPipelineTabs
#'
#' Loads the working parts of tables ("ID", allele count 2x columns, "contig" if specified); concatenates (uniting merge) tables from all provided files.
#'
#' @param inFiles A vector of full pathes to files with alelleic counts tables (necessary columns: "ID", pairs of ref and alt counts; optionally, for filtering, "contig")gene names and allelic counts
#' @param nReps A vector of numbers, either: (1) each entry is a number of replicates in the corresponding file; (2) for one file only, each entry is a number of replicates corresponding to particular experiment
#' @param contigs Parameter defining if the resulting table should be filtered by contig column (preserving only rows corresponding to a given vector), default set to empty vector() and no filtering applied
#'
#' @return A concatenated table with allele counts for all replicates, each row corresponds to a feature ("ID")
#' @export
#'
#' @importFrom "utils" "read.delim"
#'
GetGatkPipelineTabs <- function(inFiles, nReps, contigs = vector()){

  options(stringsAsFactors = FALSE)
  nameColumns <- function(nReps)  {
    if(length(nReps) != 1){
      lapply(1:length(nReps), function(i){
        paste0("rep", i, ".", rep(1:nReps[i], each = 2), c("_ref", "_alt"))
      })
    } else {
      list(paste0("rep", rep(1:sum(nReps), each = 2), c("_ref", "_alt")))
    }
  }
  if(length(contigs) != 0){
    # i.e. filtering by contig is needed
    cs_merge <- c("ID", "contig")
    if(length(inFiles) == 1){
      # contig column identification:
      c_contig <- which(names(read.delim(inFiles[1], sep="\t")) == "contig")
      cs <- list(c(1, 2:(2*sum(nReps)+1), c_contig))
      cs_names <- list(c("ID", unlist(nameColumns(nReps)), "contig"))
    } else {
      cs <- lapply(1:length(nReps), function(i){
        # contig column identification:
        c_contig <- which(names(read.delim(inFiles[i], sep="\t")) == "contig")
        c(1, 2:(2*nReps[i]+1), c_contig)
      })
      cs_names <- lapply(nameColumns(nReps), function(x){c("ID", x, "contig")})
    }
  } else {
    # i.e. NO filtering by contig
    cs_merge <- c("ID")
    if(length(inFiles) == 1){
      cs <- list(c(1, 2:(2*sum(nReps)+1)))
      cs_names <- list(c("ID", unlist(nameColumns(nReps))))
    } else {
      cs <- lapply(1:length(nReps), function(i){c(1, 2:(2*nReps[i]+1))})
      cs_names <- lapply(nameColumns(nReps), function(x){c("ID", x)})
    }
  }

  df <- Reduce(function(x,y){merge(x, y, by=cs_merge)},
               lapply(1:length(inFiles), function(i){
                 df0 <- read.delim(inFiles[i], sep="\t")
                 df0 <- df0[, cs[[i]]]
                 names(df0) <- cs_names[[i]]
                 return(df0)
               })
  )
  if(length(contigs) != 0){
    df <- df[df$contig %in% contigs, ]
    df <- df[, -which(names(df)=="contig")]
    df
  }

  return(df)
}
