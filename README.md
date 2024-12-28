
# Qllelic (v0.3.2)

**Qllelic** is a set of R tools for quantification of allele-specific expression. It relies on two or more technical replicate RNA-seq libraries to calculate Quality Correction Constant (QCC) and use it to correct for allelic imbalance overdipserion.
**Qllelic** analysis starts with a table of allelic counts per gene, calculated from RNA-seq data using any analysis pipeline such as **[ASEReadCounter*](https://github.com/gimelbrantlab/ASEReadCounter_star)**.

Mendelevich A.\*, Vinogradova S.\*, Gupta S., Mironov A., Sunyaev S., Gimelbrant A.  _"Replicate sequencing libraries are important for quantification of allelic imbalance"_, Nat Commun 12, 3370 (2021). [https://doi.org/10.1038/s41467-021-23544-8](https://www.nature.com/articles/s41467-021-23544-8)

## Installation

To install current version of this package in R:

``` r
remotes::install_github("gimelbrantlab/Qllelic")
```
Installation takes less than 1 minute. 

To install a previous version: `remotes::install_github("gimelbrantlab/Qllelic@v0.3.1")`, please take into account that it existed with different [name](https://github.com/gimelbrantlab/QCumber), so it would be installed as "QCumber".


## Manual

Full documentation: [pdf document](https://github.com/gimelbrantlab/Qllelic/wiki/documentation/Qllelic_documentation.pdf).

Please see walk-through examples on the [Wiki page](https://github.com/gimelbrantlab/Qllelic/wiki)
* for [QCC calculation and CI estimation](https://github.com/gimelbrantlab/Qllelic/wiki/Use-case-1:-One-biological-sample)
* for [AI differential analysis for two samples](https://github.com/gimelbrantlab/Qllelic/wiki/Use-case-2:-Differential-AI-analysis)

## Reporting bugs

Please report bugs to the Github [issues](https://github.com/gimelbrantlab/Qllelic/issues) page.

## License

[GNU General Public License v3.0](https://github.com/gimelbrantlab/Qllelic/blob/master/LICENSE)



