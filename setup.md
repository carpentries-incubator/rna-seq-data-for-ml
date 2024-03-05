---
title: Setup
---

## Summary

This lesson provides a practical guide to sourcing and pre-processing a bulk RNA-Seq dataset for use in a machine learning classification task. The lessons explains the characteristics of a dataset required for this type of analysis, how to search for and download a dataset from each of the main public functional genomics repositories, and then provides guidelines on how to pre-process a dataset to make it machine learning ready, with detailed examples. The lesson finally explains some of the additional data filtering and transformation steps that will improve the performance of machine learning algorithms using RNA-Seq count data.

The lesson is written in the context of a supervised machine learning classification modelling task, where the goal is to construct a model that is able to differentiate two different disease states (e.g. disease vs. healthy control) based on the gene expression profile.

This work was funded by the ELIXIR-UK: FAIR Data Stewardship training UKRI award (MR/V038966/1)

![Elixir-UK](fig/ELIXIR-UK_logo.png){alt='Elixir-UK'}

:::::::::::::::::::::::::::::::::::::: prereq 

This lesson assumes a working knowledge of programming in R. For learners who aren't familiar with R or feel they need a refresher, the [Programming with R](https://swcarpentry.github.io/r-novice-inflammation/index.html) provides a good introduction to both R and working with R studio.

:::::::::::::::::::::::::::::::::::::::::::::

## Data Sets

This lesson uses a number of RNA-Seq datasets downloaded from public functional genomics repositories.

Before the start of the lesson, create a new RStudio project where you will keep all of the files for this lesson. In the project directory (where the .Rproj file is), **create a subdirectory called `data`**.

Links to the relevant datasets and instructions on how to download them provided in each episode.


## Software Setup

::::::::::::::::::::::::::::::::::::::: discussion

### R and RStudio

Learners will need updated versions of R and RStudio installed. There are instructions for installing R, RStudio, and additional R packages for all the main operating systems in the [R Ecology Lesson](https://datacarpentry.org/R-ecology-lesson/#Install_R_and_RStudio).

Please install the following R packages. You will need to install the package `BiocManager` to be able to install packages from Bioconductor.

* `tidyverse`
* `reshape2`
* `scales`
* `BiocManager`
* `Biobase`
* `GEOquery`
* `DESeq2`

Executing the following lines of code in the R console will install all of these packages.

```r

install.packages(c("tidyverse", "reshape2", "scales", "BiocManager"))

BiocManager::install(c("Biobase", "GEOquery", "DESeq2"))

```

You can check that all packages have been installed using the following command, which will return `character(0)` if all packages have been successfully installed.

```r

setdiff(c("tidyverse", "reshape2", "scales", "BiocManager", "Biobase", "GEOquery", "DESeq2"),
        rownames(installed.packages()))

```
::::::::::::::::::::::::::::::::::::::: 
