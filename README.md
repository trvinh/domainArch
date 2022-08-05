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

<img width="1292" alt="image" src="https://user-images.githubusercontent.com/19269760/183078583-a973d96b-4a62-49d7-8a79-408521c7add3.png">
