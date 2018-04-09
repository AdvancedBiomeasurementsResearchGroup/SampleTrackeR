# SampleTrackeR

This repository described SampleTrackeR, an R script for sample assurance in multiplexed sequencing experiments, based on tagging of samples with synthetic spike-in control mixtures (STMs).

## Prerequisites

`SampleTrackeR` depends on several external R packages that need to be installed and loaded.

```
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(tibble)
```

The `SampleTrackeR` script is not part of an R package and the simplest way to load the script is by sourcing the code as follows:

```
source("/absolute/path/to/SampleTrackeR.R")
```

Note that the `SampleTrackeR` script sets the current working directory as the path when searching for input files; all input files (see below for details) thus need to be present in the current working directory.

## Description of input files

### sample_plate_layout

This tab-delimited file describes the experiments, namely plate layout and STM added to each sample / sequencing library. Samples lacking STMs can be added using a mock STM (designated *e.g.*, stm00); the mock STM compoistion should also be present in the stm_compositions file (see below for details).

The following columns and matching names are required.

  + `libID`: name of the sample or sequencing library.

  +  `stmID`: identifier of the STM to added to the sample.

  +  `row`: row identifier (should be integer).

  +  `column`: column identifier (should be integer).

Other columns can be also added (e.g. description in the example table below); these are however ignored and not included in any of the generated output files.

### read_count_table

This tab-delimited file represent a typical OTU read count table. A column with name otuID is mandatory. Addition column names represent sample identifiers as in `sample_plate_layout`.

### stm_compositions

This tab-delimited file represents a long-format table with the composition of the STMs.

The following columns as well as their names are required.

  +  `stmID`: identifier of the STM
  
  +  `controlID`: identifier of the spike-in control

  +  `value`: value indicating wether the spike-in control is present in the STM (1: present and 0: absent)

To allow analysis of samples lacking spike-in controls, a mock STM (*e.g.*, stm00) can be added with all values set to 0.

## General usage, output and terminology

```
out <- SampleTrackeR(sample_plate_layout = "sample_plate_layout.txt",
                     read_count_table = "read_count_table.txt", 
                     stm_compositions = "stm_compositions.txt",
                     read.threshold = 1,
                     fraction.threshold = 1)
```

### Description input arguments

  +  `sample_plate_layout` (required) Name of tab-delimited file of sample layout.
  
  +  `read_count_table` (required) Name of tab-delimited file with read count data, including both synthetic spike-in controls and environmental OTUs.
  
  +  `stm_compositions` (required) Name of tab-delimited file with STM compositions.
  
  +  `read.threshold` (optional) Minimum number of reads for a given spike-in control to be scored as present (default 1).
  
  +  `fraction.threshold` (optional) Minimum fraction spike-in controls for a given STM to be scored as present (default 1).

### Description outputs

The output of `SampleTrackeR` is a list with four different objects. Assuming that the output list is called `out`, the resultant list contains the following:

  +  `out$tab1`, a data frame containing a summary of the sample identification based on majority STMs.

  +  `out$plot1`, a ggplot2 object visualizing of the sample identification.

  +  `out$tab2`, a data frame containing a summary of between-sample carry-over based on minority STMs.

  +  `out$plot2`, a ggplot2 object visualizing the output of the between-sample carry-over.
  
Note that `out$tab2` and `out$plot2` are only generated when minority STMs are identified in at least one of the samples. In case no minority STMs are found, a message will be output stated that between-sample carry-over was not evaluated, and hence assurmed to be minimal based on the lack of minority STMs.

###  Terminology

`majority_STM`, refers to the STM with the highest cumulative read count in a given sample. The majority STM is used to assign sample identity and detect/resolve potential sample swaps.

`minority_STM`, refers to STMs present in a sample, excluding the majority STM. A sample can contain multiple minority STMs and these are used for assessment of between-sample carry-over (that is, cross-contamination).

`distinguishing_controls`, refers to the number of individual spike-in controls that are not shared between two different STMs. For evaluation of between-sample carry-over, the number of distinguishing controls is the number of spike-in controls present in the minority STM but absent in the majority STM.

`percent_carryover`, refers to the estimated amount of sample carry-over between two samples, as quantified based on minority STMs.


