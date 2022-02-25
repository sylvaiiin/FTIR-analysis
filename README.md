# FTIR-analysis
R scripts for FTIR analysis. 
This scripts permit selection, analysis and exportation of the FTIR data.

## How to get started ?
### Package installations
Install the following packages if that is not already done
- ggplot2
- dplyr
- reshape
- cowplot
- stringr
- fs

### Run DataTreatment_FTIR.R
Run DataTreatment_FTIR.R to begin the selection and analysis process

### Data type
This script use **.brut** or **.csv** data coming from FTIR.

### Entry parameters
You can change the entry parameters in **DataTreament_FTIR.R** (line 39 to 50) to suit your needs.

## What is the scripts doing ?
### Selection
The first phase is a phase of spectrum selection. You can go through every spectrum that you took at the FTIR. Then you can delete the unusable spectrums from your dataset.

### Exportation
A csv file is created everytime you run the script. You can find your clean dataset  in this file.
**Be careful to rename your file before passing it in the software to not overrun it**

### Analysis
The second phase is a phase of analysis of your dataset. you have 4 options to visualize your data :
- Graph with all spectrums
- Graph with the mean spectrums of population
- Graph with the mean spectrum of one population
- Mean spectrum comparison between two populations using a T-test


Copyright Urbain Aurelie, Caillaud Sylvain & INRAE
