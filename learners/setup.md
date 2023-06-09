---
title: Setup
---


## Data Sets

This lesson uses a number of RNA-Seq datasets downloaded from public functional genomics repositories. The links to the datasets are provided in each episode. Create a new RStudio project where you will keep all of the files for this lesson. In the project directory (where the .Rproj file is), create a subdirectory called data.

Either download each file as you go along and save to you data subdirectory, or alternatively you can download all the files used in the lesson in this 

Download the [data zip file](data/rna-seq-ml-readiness-data.zip) and unzip it to your data subdirectory.

## Software Setup

::::::::::::::::::::::::::::::::::::::: discussion

### R and RStudio

Learners will need updated versions of R and RStudio installed. There are instructions for installing R, RStudio, and additional R packages for all the main operating systems in the [R Ecology Lesson](https://datacarpentry.org/R-ecology-lesson/#Install_R_and_RStudio).

### R Packages

Please install the following R packages. You will need to install the package `BiocManager` to be able to install packages from Bioconductor.

* `tidyverse`
* `here`
* `BiocManager`
* `Biobase`
* `GEOquery`
* `DESeq2`

Executing the following lines of code in the R console will install all of these packages

```{r}

install.packages(c("tidyverse", "here", "BiocManager"))

BiocManager::install(c("Biobase", GEOquery", "DESeq2"))

```

You can check that all packages have been installed using the following command, which will return `character(0)` if all packages have been successfully installed.

```{r}

setdiff(c("tidyverse", "here", "BiocManager", "GEOquery", "DESeq2"),
        rownames(installed.packages()))

```
::::::::::::::::::::::::::::::::::::::: 
