# ssGSEA2.0/PTM-SEA

This repository provides an R package implementation of [ssGSEA2.0](https://github.com/broadinstitute/ssGSEA2.0).
See the [upstream repository](https://github.com/broadinstitute/ssGSEA2.0) for more details about 
publications, ssGSEA 2.0, PTM-SEA, [MSigDB](http://software.broadinstitute.org/gsea/msigdb/), and 
[PTMsigDB](https://proteomics.broadapps.org/ptmsigdb/). 

## Installation

### R >= 4.0

```r
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
devtools::install_github("nicolerg/ssGSEA2")
```

### R 3.6

You must install the GitHub version of `cmapR` **first**:  
```r
if (!require("devtools", quietly = TRUE)){
  install.packages("devtools")
}
install_github("cmap/cmapR", ref="R-3.6")
```

Then install this package:  
```r
devtools::install_github("nicolerg/ssGSEA2")
```

## Example 

In this example, we perform PTM-SEA with example input provided in [examples/](examples). 
```r
library(ssGSEA2)

# Download example input
download.file(url = "https://raw.githubusercontent.com/nicolerg/ssGSEA2/master/example/PI3K_pert_logP_n2x23936.gct",
              destfile = "/tmp/PI3K_pert_logP_n2x23936.gct")

# Download gene set database 
download.file(url = "https://raw.githubusercontent.com/nicolerg/ssGSEA2/master/example/ptm.sig.db.all.flanking.human.v1.8.1.gmt"),
              destfile = "/tmp/ptm.sig.db.all.flanking.human.v1.8.1.gmt")

res = run_ssGSEA2("/tmp/PI3K_pert_logP_n2x23936.gct",
                  output.prefix = "example",
                  gene.set.databases = "/tmp/ptm.sig.db.all.flanking.human.v1.8.1.gmt",
                  output.directory = "/tmp",
                  sample.norm.type = "none", 
                  weight = 0.75, 
                  correl.type = "rank", 
                  statistic = "area.under.RES",
                  output.score.type = "NES", 
                  nperm = 1000, 
                  min.overlap = 5, 
                  max.overlap = 2000,
                  extended.output = TRUE, 
                  global.fdr = FALSE,
                  log.file = "/tmp/run.log")
```

## Notes 

We are aware that the following warnings are seen upon attaching this package:  
```
Warning messages:
1: multiple methods tables found for ‘aperm’ 
2: replacing previous import ‘BiocGenerics::aperm’ by ‘DelayedArray::aperm’ when loading ‘SummarizedExperiment’
```
This could be avoided by requiring `DelayedArray >= 0.24.0` and `R >= 4.2`. 
However, to make this package accessible to users using older versions of R, we opted not to do this. 
See more details [here](https://github.com/cmap/cmapR/issues/70). 

## Citing this work 

### Citing ssGSEA2.0

Krug, K., Mertins, P., Zhang, B., Hornbeck, P., Raju, R., Ahmad, R., . Szucs, M., 
Mundt, F., Forestier, D., Jane-Valbuena, J., Keshishian, H., Gillette, M. A., Tamayo, 
P., Mesirov, J. P., Jaffe, J. D., Carr, S. A., Mani, D. R. (2019). 
**A curated resource for phosphosite-specific signature analysis**, 
Molecular & Cellular Proteomics, 18(3), 576-593. 
http://doi.org/10.1074/mcp.TIR118.000943

### Citing ssGSEA

Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005).
**Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.**
Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550. 
http://doi.org/10.1073/pnas.0506580102

Abazeed, M. E., Adams, D. J., Hurov, K. E., Tamayo, P., Creighton, C. J., Sonkin, D., et al. (2013).
**Integrative Radiogenomic Profiling of Squamous Cell Lung Cancer.** Cancer Research, 73(20), 6289–6298.
http://doi.org/10.1158/0008-5472.CAN-13-1616
