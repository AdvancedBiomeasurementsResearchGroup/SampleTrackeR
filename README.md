## Contents

This repository describes `SampleTrackeR`, an R script for sample assurance in multiplexed sequencing experiments, based on tagging of samples with synthetic spike-in control mixtures (STMs).

The following files are available in this repository.

1. "SampleTrackeR.R", source code.

2. Data files for results presented in Tourlousse *et al.* (see citation), namely "sample_plate_layout1.txt", "sample_plate_layout2.txt", "sample_plate_layout3.txt", "read_count_table1.txt", "read_count_table2.txt", "read_count_table3.txt" and "stm_compositions.txt".

3. "SampleTrackeR.nb.html", an HTML file demonstrating the usage of `SampleTrackeR` using the above data sets.

## Using SampleTrackeR

### Prerequisites

`SampleTrackeR` depends on a small number of external R packages that need to be installed and loaded.

```
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
```

`SampleTrackeR` consists of a single script that can be loaded by sourcing the code.

```
source("/absolute/path/to/SampleTrackeR.R")
```

Note that the script automaticcaly sets the current directory as the path to search for input files (internally, `path <- getwd()`). All input files thus need to be present in the current directory.

### Description of input files

#### sample_plate_layout

This tab-delimited file describes the experiments, namely plate layout and STM added to each sample / sequencing library. Samples lacking STMs can be added using a *mock* STM (designated *e.g.*, stm0); the mock STM composition should also be present in the `stm_compositions` file (see below for details).

The following columns and matching names are required.

  + `libID`: name of the sample or sequencing library.

  + `stmID`: identifier of the STM to added to the sample.

  + `row`: row identifier (should be integer).

  + `column`: column identifier (should be integer).

Other columns can be present (*e.g.*, *description* in the example table below); these are however ignored and not included in any of the generated output files.

| libID | stmID | row | column | description |
| ------|-------|-----|--------| --------|
| lib1 | stm1 | 1 | 1 | soil_stm1 |
| lib2 | stm2 | 2 | 1 | sludge_stm2 |
| lib3 | stm3 | 3 | 1 | feces_stm3 |
| ... | ...  | ... | ... |
| lib4 | stm0 | 6 | 12 | soil |
| lib5 | stm0 | 7 | 12 | sludge |
| lib6 | stm0 | 8 | 12 | feces |

#### read_count_table

This tab-delimited file represent a typical OTU read count table. A column with name `otuID` is mandatory. The other column names represent sample identifiers as in the `sample_plate_layout` file.

| otuID | lib1 | lib2 | ... |
| ------|-------|-----|--------|
| control1 | 101 | 231 | ... |
| control2 | 3  | 10 | ... |
| control3 | 0 | 1 | ... |
| ... | ...  | ... | ... |
| OTU1 | 0  | 5 | ... |
| OTU2 | 6  | 567 | ... |
| ... | ...  | ... | ... |

#### stm_compositions

This tab-delimited file represents a long-format table with the composition of the STMs.

The following columns and matching names are required.

  + `stmID`: identifier of the STM
  
  + `controlID`: identifier of the spike-in control

  + `value`: value indicating wether the spike-in control is present in the STM (1: present and 0: absent)
  
| stmID | controlID | value | 
| ------|-------|-----|
| stm1 | control1 | 1 | 
| stm1 | control2  | 1 | 
| stm1 | control3  | 0 | 
| stm2 | control1 | 1 | 
| stm2 | control2  | 0 | 
| stm2 | control3  | 1 |  
| stm3 | control1 | 0 | 
| stm3 | control2  | 1 | 
| stm3 | control3  | 1 | 
| ... | ...  | ... | 
| stm0 | control1  | 0 | 
| stm0 | control2  | 0 | 
| stm0 | control3  | 0 | 

For the mock STM for samples without added STMs, all values should be set to 0.

### General usage, output and terminology

```
out <- SampleTrackeR(sample_plate_layout = "sample_plate_layout.txt",
                     read_count_table = "read_count_table.txt", 
                     stm_compositions = "stm_compositions.txt",
                     read.threshold = 1,
                     fraction.threshold = 1)
```

#### Description input arguments

  + `sample_plate_layout` (*required*) Name of tab-delimited file of sample layout.
  
  + `read_count_table` (*required*) Name of tab-delimited file with read count data, including both synthetic spike-in controls and sample OTUs.
  
  + `stm_compositions` (*required*) Name of tab-delimited file with STM compositions.
  
  + `read.threshold` (*optional*) Minimum number of reads for a given spike-in control to be scored as present (default 1).
  
  + `fraction.threshold` (*optional*) Minimum fraction of spike-in controls for a given STM to be scored as present in order for the STM to be scored as present (default 1).

#### Description outputs

The output of `SampleTrackeR` is a list with four different objects. Assuming that the output list is called `out`, the resultant list contains the following:

  + `out$tab1`, a data frame containing a summary of sample identification based on majority STMs.

  + `out$plot1`, a ggplot2 object visualizing the output of the sample identification portion of the script.

  + `out$tab2`, a data frame containing a summary of between-sample carry-over based on minority STMs.

  + `out$plot2`, a ggplot2 object visualizing the output of the between-sample carry-over portion of the script.
  
Note that `out$tab2` and `out$plot2` are *not* generated when no minority STMs are identified in any of the samples. In this case, a message will be output stated that between-sample carry-over was not evaluated, and hence assurmed to be minimal based on the lack of minority STMs.

####  Terminology

  + `majority_STM`, refers to the STM with the highest cumulative read count in a given sample. The majority STM is used to assign sample identity and detect/resolve potential sample swaps.

  + `minority_STM`, refers to STMs present in a sample, excluding the majority STM. A sample can contain multiple minority STMs and these are used for assessment of between-sample carry-over (that is, cross-contamination).

  + `distinguishing_controls`, refers to the number of individual spike-in controls that are not shared between two different STMs. For evaluation of between-sample carry-over, the number of distinguishing controls is the number of spike-in controls present in the minority STM but absent in the majority STM.

  + `percent_carryover`, refers to the estimated amount of sample carry-over between two samples, as quantified based on minority STMs.

### Citation

**Dieter M. Tourlousse, Akiko Ohashi, Yuji Sekiguchi.** Sample tracking in microbiome community profiling assays using synthetic 16S rRNA gene spike-in controls. *Manuscript in revision.*

### Contact

Dieter Tourlousse: dieter.tourlousse@aist.go.jp

