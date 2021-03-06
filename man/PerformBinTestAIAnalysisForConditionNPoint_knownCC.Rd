% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/PerformBinTestAIAnalysisForConditionNPoint_knownCC.R
\name{PerformBinTestAIAnalysisForConditionNPoint_knownCC}
\alias{PerformBinTestAIAnalysisForConditionNPoint_knownCC}
\title{PerformBinTestAIAnalysisForConditionNPoint_knownCC}
\usage{
PerformBinTestAIAnalysisForConditionNPoint_knownCC(
  inDF,
  vectReps,
  vectRepsCombsCC,
  pt = 0.5,
  Q = 0.95,
  thr = NA,
  thrUP = NA,
  thrType = "each",
  minDifference = NA
)
}
\arguments{
\item{inDF}{Allele counts dataframe: with 2n+1 columns, "ID" and 2n columns with ref & alt counts (rep1_ref, rep1_alt, rep2_ref, rep2_alt, ...)}

\item{vectReps}{A vector (>=2) of replicate numbers for which the analysis should be applied}

\item{vectRepsCombsCC}{A vector of pairwise-computed correction constants for given replicates (QCC=1 is no correction)}

\item{pt}{Optional (default=0.5), a value to compare with}

\item{Q}{Optional (default=0.95), confidence level, quantile}

\item{thr}{Optional (default=NA), threshold on the overall number of counts for a gene to be considered in the analysis}

\item{thrUP}{Optional (default=NA), threshold for max gene coverage (default = NA)}

\item{thrType}{Optional (default = "each", also can be "average" for average coverage on replicates), threshold type}

\item{minDifference}{Optional (default=NA), if specified, one additional column is added to the output (T/F depending if the gene changed AI expression more than minDifference in addition to passing the test)}
}
\value{
A table of gene names, AIs + CIs, p-values for both non-corrected and (BT..) and QCC corrected (BT_CC..) tests, classification into genes demonstrating signifficant difference (TRUE) from point estimate AI and those that don't (FALSE).
}
\description{
Performs Binomial and QCC-corrected binomial tests (with Bonferroni correction) with a given point estimate, for given QCC.
}
