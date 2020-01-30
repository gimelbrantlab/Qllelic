#' BuildDesign
#'
#' Creates a design matrix for the experiment
#'
#' @param experimentNames Vector with names of the experiments
#' @param techReps Vector with number of technical replicates in each experiment
#' @param corrConst Optional, Vector with correction constants for each experiment
#'
#' @return Dataframe with experiments numbered and numbers of columns
#' @export
#'
#' @examples BuildDesign(c("clone1", "clone2", "clone3"), c(2, 2, 3), c(1.7, 1.8, mean(1.55, 1.6, 1.57)))
#'
BuildDesign <- function(experimentNames, techReps, corrConst=NA){
  rowsSp <- data.frame(matrix(lapply(1:length(techReps), function(x){
              (2*sum(techReps[1:x-1])+1):(2*sum(techReps[1:x]))
            }),
            nrow=length(techReps), byrow=T), stringsAsFactors=FALSE)
  colnames(rowsSp) <- "replicateCols"

  colExp <- data.frame(matrix(lapply(1:length(techReps), function(x){
              (sum(techReps[1:x])-techReps[x]+1):(sum(techReps[1:x]))
            }),
            nrow=length(techReps), byrow=T), stringsAsFactors=FALSE)
  colnames(colExp) <- "replicateNums"

  designMatrix <- cbind(experimentNames, techReps, rowsSp, colExp)
  colnames(designMatrix) <- c("experimentNames", "techReps", "replicateCols", "replicateNums")

  if (!sum(is.na(corrConst))) {
    designMatrix <- cbind(designMatrix, corrConst)
  }
  return(designMatrix)
}
