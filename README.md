# InProgress Overview
This pacakge is designed to take MANTA or SvABA VCF files which have been run without a paired normal and classify each rearrangement as germline or somatic. 

# System Requirements and Dependencies
- This has been tested on the current versions of OSX, Linux, Ubuntu, Windows (as of April 1, 2022)

# Installation
1. Install R-4.0 or newer (https://www.r-project.org/)
2. Install necessary packages
  ```{r}
  install.pacakges('devtools')
  install.packages('parallel')
  install.packages('BiocManager')
  
  library(BiocManager)
  BiocManager::install('GenomicRanges')
  
  library(devtools)
  devtools::install_github('mskilab/gUtils')
  ```
  Documentation for these packages can be found at: \
    [devtools](https://cran.r-project.org/web/packages/devtools/index.html) \
    [parallel](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf) \
    [BiocManager](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html) \
    [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) \
    [gUtils](https://github.com/mskilab/gUtils) 
  
3. Install InProgress
```bash
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE) 
devtools::install_github('acranej/InProgress')
```
