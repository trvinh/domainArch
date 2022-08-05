# domainArch

From an R terminal, PhyloRBF can be installed using *devtools*:

```r
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("trvinh/domainArch", INSTALL_opts = c('--no-lock'), dependencies = TRUE)
```

Then, to run it, enter:

```r
library(PhyloRBF)
runPhyloRBF()
```
