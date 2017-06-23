# mutationCheck
Check and QC mutations from BAM

Manually check mutations by scanning the original BAM file, and recalculate supporting reads and the variant allele frequency(VAF). 
Basicall, two steps are performed, 1. do QC on a read; 2. do QC on a mutant bases. This is now under very naive deveploped, and may 
just be used for temporary in-house use.

## Installation
A demo package could be installed directly from github:
```{r}
install.packages ("devtools")
library(devtools)
install_github("Hiuyu/mutationCheck")
```

