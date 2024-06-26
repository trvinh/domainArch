# domainArch

## Installation and usage
From an R terminal, **domainArch** can be installed using *devtools*:

```r
if (!requireNamespace("devtools"))
    install.packages("devtools")
devtools::install_github("trvinh/domainArch", INSTALL_opts = c('--no-lock'), dependencies = TRUE)
```

Or you can download the source code and install it with the command:
```r
install.packages(path_to_domainArch, repos = NULL, type="source")
```

To run **domainArch**, enter:

```r
library(domainArch)
runDomainArch()
```

## Input

Currently the tool accepts 3 kinds of input:
- A folder containing several **domain files**, each file for a single seed protein together with its orthologs. These files are the `*.domains` file generated by *[fdog.run](https://github.com/BIONF/fDOG/wiki/Input-and-Output-Files#output-files)*. *__Please note:__* *each domain file must be named following this format `<seed/groupID>.domains`*
- A single **concatenated domain file** for multiple seeds and ortholog proteins. This file is the output of *[fdogs.run](https://github.com/BIONF/fDOG/wiki/Use-standalone-version)* or generated by combining multiple single fDOG outputs using *fdog.mergeOutput* function
- A folder containing **annotation files** for several species in JSON format, e.g. *annotation_dir* in [fDOG data](https://github.com/BIONF/fDOG/wiki/Input-and-Output-Files#data-structure). These annotation files can be created using *fas.doAnno* function of [FAS tool](https://github.com/BIONF/FAS/wiki/Annotation)

In the future one can use directly the protein sequences as input. Stay tuned! ;)

Tip: if your concatenated domain file is too large, you should split it into multiple singe domain files by using the function *Split domain file* provided by **domainArch**.


## Screenshot

![image](https://github.com/trvinh/domainArch/assets/19269760/8bdea154-08c1-4648-bd47-55c4507ce3ae)

![image](https://github.com/trvinh/domainArch/assets/19269760/e7f78bda-dac1-4f72-b05e-9a452d6f0818)

![image](https://github.com/trvinh/domainArch/assets/19269760/5977c79c-1125-4f70-a7dd-3c5ab5ea9782)

## Contact

For any bug reports or questions, please [open an issue on GitHub](https://github.com/trvinh/domainArch/issues/new) or be in touch via email tran@bio.uni-frankfurt.de
