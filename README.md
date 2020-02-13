
# QCumber (v0.3.2)

R tools for quantification of allelic imbalance within one sample and for differential analysis of allelic imbalance between two or more samples. Relies on comparison of allele-specific signal from technical replicate RNA-seq libraries to compensate for technical noise.  
Input: tables of allelic counts, calculated from RNA-seq data using sequencing analysis pipeline such as [ASEReadCounter*](https://github.com/gimelbrantlab/ASEReadCounter_star).

## Installation

To install this package in R:

``` r
devtools::install_github("gimelbrantlab/QCumber")
```

## Usage

Please see [Wiki page](https://github.com/gimelbrantlab/QCumber/wiki) for worked exaples:
* for [QCC calculation and CI estimation](https://github.com/gimelbrantlab/QCumber/wiki/Use-case-1:-One-biological-sample)
* for [AI differential analysis for two samples](https://github.com/gimelbrantlab/QCumber/wiki/Use-case-2:-Differential-AI-analysis)

## Citations

Please cite _"Unexpected variability of allelic imbalance estimates from RNA sequencing", Mendelevich A.*, Vinogradova S.*, Gupta S., Mironov A., Sunyaev S., Gimelbrant A._, if you used our R-package in your work.

## Reporting bugs

Please report bugs to the Github [issues](https://github.com/gimelbrantlab/QCumber/issues) page.

## License

[GNU General Public License v3.0](https://github.com/gimelbrantlab/QCumber/blob/master/LICENSE)



