---
title: "Data Collection: RNA-Seq Datasets for ML Analysis"
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- Where can I find an RNA-Seq dataset for a ML analysis?
- What characteristics make data appropriate for ML?
- What format is RNA-Seq data stored in?
- Do I use 'raw' or 'processed' data?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Recall the main platforms hosting functional genomic datasets
- Recall the features of a gene expression dataset that make it more appropriate for ML
- Explain the difference between raw and processed RNA-Seq data
- Recall the different stages of processing observed in public datasets

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

This lesson provides a practical guide to sourcing and pre-processing a bulk RNA-Seq dataset for use in a machine learning classification task. The lessons explains the characteristics of a dataset required for this type of analysis, how to search for and download a dataset from each of the main public functional genomics repositories, and then provides guidelines on how to pre-process a dataset to make it machine learning ready, with detailed examples. The lesson finally explains some of the additional data filtering and transformation steps that will improve the performance of machine learning algorithms using RNA-Seq count data.

The lesson is written in the context of a supervised machine learning classification modelling task, where the goal is to construct a model that is able to differentiate two different disease states (e.g. disease vs. healthy control) based on the gene expression profile.


## Functional Genomics platforms

There are two major public repositories for sourcing public functional genomics data sets, in particular microarry and next-generation sequencing data. Here we focus on NGS data, in particular RNA-Seq transcriptomics datasets:

1. The [Array Express](https://www.ebi.ac.uk/biostudies/) collection within the Bio Studies database, maintained by the European Molecular Biology Laboratory - European Bioinformatics Institute (EMBL-EBI).  

2. The [Gene Expression Onmibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/), maintained by the National Center for Biotechnology Information, part of the US National Institues of Health. 

Each database contains thousands of RNA-Seq datasets relating to a wide variety of experiments, not all of which would be suitable for analysis with a machine learning approach.


## Dataset Characteristics for supervised learning

The success of a ML/AI model depends heavily on the input data. Finding an appropriate dataset can be a challenge and the data must be selected carefully as a predictive model trained on the unreliable or inappropriate data will produce unreliable predictions.

There is currently no curated "ML/AI ready" datasets that meet the requirements of machine learning analyses within the ArrayExpress or GEO repositories. It is therefore important to conduct careful searches and examine the data sets identified to confirm their appropriateness for a supervised machine learning task. The following characteristics of the data should be considered:

* Data provenance. The source of the data, its likely quality and whether it is recognised by the community, for example having been used in highly regarded publications.
* Dataset size. Ensure the data set captures the full complexity of the underlying distribution, for example, the biological heterogeneity in the sample population, and that there are sufficient samples to create independent train, validation and tests sets within the dataset. Bulk RNA-Seq datasets with only a few samples for each group are unlikely to be sufficient as a training dataset. If examples are grouped in classes, there must be sufficient examples available in each class.
* Metadata Quality. Training data in classification task require accurate and reliable labels. Datasets must have sufficient metadata, ideally connected to a domain-specific or community-specific ontology so that the meaning is widely interpretable and comparable between studies.

Further information on data considerations in supervised machine learning with biological data can be found in the [DOME recommendations](https://www.nature.com/articles/s41592-021-01205-4).

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: instructor

Add instructor notes here

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


## Raw vs. Processed Data

Data sets will contain both raw and processed RNA-Seq data data sets, as well as information about the relationship between the RNA-Seq data and study samples, and the design of the experiment.

* Raw RNA-Seq data consists of the sequences of each read obtained from the sequencer and is stored in FastQ file format. Processed RNA-Seq data is typically stored as a matrix counts, with the names of each sample on one axis, and the unique identifiers of each transcript identified on the other axis. For more information on FastQ file format, and the transformation between FastQ and counts data see [reference previous elixir training/ RDM bites]. Processed RNA-Seq data is stored in a range of formats, typically as csv or txt files.

* Processed data refers to gene abundance data, which may be interger counts of the number of transcripts mapped to a given gene, or a further processed measure of abundance to normalise for sequence depth and/or gene length such as TPM and FPKM.

A quick summary of raw vs. processed data:

Data Type     | Raw | Processed
------------- | -------------
Contains  |  Strings of based pair sequence | Abundance value by gene (integer, float)
File Format  | fastq | .csv, .txt, .tsv
                

The description processed means different things in different studies. It is important therefore to read the study protocol to determine how the data has been processed.

For more information on the fastq file format, see [RDMBites | FASTQ Format](https://www.youtube.com/watch?v=tO2H3zuBouw). For further information on the processing of raw fastq to gene abundance data, see [RDMbites | RNAseq expression data](https://www.youtube.com/watch?v=3Pe9xcGF_Wo).

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1:

What type of data do you need to look for on a gene expression data repository as input to a machine learning modelling analysis?

:::::::::::::::::::::::: solution 

The starting point for a machine learning model is processed RNA-Seq data, in the form of count data. This could be raw count data, or TPM or FPKM pre-processed data. The data will likely need further pre-processing before input to a ML model as discussed in Episode X.

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::


::::::::::::::::::::::::::::::::::::: keypoints 

- Consider dataset provenance, dataset size and label quality when selecting a RNA-Seq dataset for ML
- The two main repositories for RNA-Seq dataset are ArrayExpress and GEO
- Processed data, in the form of raw counts, or further processed counts data is the starting point for analysis

::::::::::::::::::::::::::::::::::::::::::::::::

[r-markdown]: https://rmarkdown.rstudio.com/
