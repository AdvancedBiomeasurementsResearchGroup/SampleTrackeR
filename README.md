# SampleTrackeR

This repository described SampleTrackeR, an R script for sample assurance in multiplexed sequencing experiments, based on tagging of samples with synthetic spike-in control mixtures (STMs).

## Prerequisites

`SampleTrackeR` depends on several external R packages that need to be installed and loaded.

```
library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(tibble)
```

The `SampleTrackeR` script is not part of an R package and the simplest way to load the script is by sourcing the code as follows:

```
source("/absolute/path/to/SampleTrackeR.R")
```

Note that the `SampleTrackeR` script sets the current working directory as the path when searching for input files; all input files (see below for details) thus need to be present in the current working directory.
