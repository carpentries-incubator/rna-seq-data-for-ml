---
title: "Introduction: Machine Learning Ready RNA-Seq Data"
teaching: 10
exercises: 1
---

:::::::::::::::::::::::::::::::::::::: questions 

- What characteristics of a dataset do I need to consider to make it 'ready' for a machine learning /AI modelling analysis?
- Where can I find a publicly available RNA-Seq dataset suitable for a machine learning analysis?
- What format is RNA-Seq data stored in on public repositories?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Describe the factors that must be considered in readying an RNA-Seq dataset for a machine learning analysis?
- Recall the main platforms hosting functional genomic datasets
- Explain the difference between raw and processed RNA-Seq data, and the different stages of processing for RNA-Seq data

::::::::::::::::::::::::::::::::::::::::::::::::


## Dataset Characteristics for supervised learning

The success of a ML/AI model depends heavily on the input data. Finding an appropriate dataset can be a challenge and the data must be selected carefully as a predictive model trained on the unreliable or inappropriate data will produce unreliable predictions.

There is currently no curated "ML/AI ready" datasets that meet the requirements of machine learning analyses within public functional genomics repositories. It is therefore important to examine a data set carefully to confirm its appropriateness for a supervised machine learning task. The following characteristics of the data should be considered:


Characteristic | Considerations
--- | ---------
**Quality and Provenance** |  Assess source of the data, its likely quality and whether it is recognised by the community, for example having been used in highly regarded publications. Datasets must have sufficient metadata, including accurate and reliable class labels, ideally connected to a domain-specific or community-specific ontology so that the meaning is widely interpretable and comparable between studies.
**Size**| Ensure the data set captures the full complexity of the underlying distribution. If there is significant biological heterogeneity in the sample population, will there be sufficient samples to create independent and representative train, validation and tests sets? Bulk RNA-Seq datasets with only a few samples for each group are unlikely to be sufficient to train a machine learning classifier.
**Format and Integrity** | Machine learning models require data to be machine readable and input in a specific format. Data collected experimentally and data acquired from public repositories will likely need to be reformatted and checked for integrity issues.
**Technical Noise** | RNA-Seq data is likely to contain some technical noise resulting from experiential processes that is unrelated to the underlying biology. This has the potential to bias machine learning algorithms and should be identified and ideally removed.
**Distribution and Scale** | Many machine learning algorithms are sensitive to the distribution and the scale of the input data. RNA-Seq data must be transformed and re-scaled to perform well in a machine learning analysis.

<br>

In this lesson, we will address each of these issues in turn, working through the collection of a RNA-Seq dataset, reformatting, integrity checking, noise elimination, transformation and rescaling to make the data ready for machine learning. Further information on some of the data considerations in supervised machine learning with biological data can be found in the [DOME recommendations](https://www.nature.com/articles/s41592-021-01205-4).


## Functional Genomics platforms

There are two major public repositories for sourcing public functional genomics data sets, in particular microarray and next-generation sequencing data. Here we focus on NGS data, in particular RNA-Seq transcriptomics datasets:

1. The [Array Express](https://www.ebi.ac.uk/biostudies/arrayexpress/studies) collection within the Bio Studies database, maintained by the European Molecular Biology Laboratory - European Bioinformatics Institute (EMBL-EBI).  

2. The [Gene Expression Onmibus (GEO)](https://www.ncbi.nlm.nih.gov/geo/), maintained by the National Center for Biotechnology Information, part of the US National Institues of Health. 

Each database contains thousands of RNA-Seq datasets relating to a wide variety of experiments, not all of which would be suitable for analysis with a machine learning approach.



## Raw vs. Processed Data

Public datasets will contain both raw and processed RNA-Seq data, as well as information about the relationship between the RNA-Seq data and study samples, and the design of the experiment. Here is a brief summary of the difference between raw and processed RNA-Seq data:

<br>

Data Type         | What it is | File Format
--                 | --------     | --
**Raw data**      |  Sequences of each read obtained from the sequencing instrument, along with quality scores for each sequence | FASTQ
**Processed data**|  Matrix of abundance values (e.g. integer counts) for the feature of interest, which may be genes, transcripts, exons or miRNAs | .csv / .txt / .tsv 

For more information on the fastq file format, see [RDMBites | FASTQ Format](https://www.youtube.com/watch?v=tO2H3zuBouw) and for information on the processing of raw fastq to gene abundance data, see [RDMbites | RNAseq expression data](https://www.youtube.com/watch?v=3Pe9xcGF_Wo).

<br>

:::::::::::::::::::::::::::::::::::::: callout

The description "processed" means different things in different studies. Things to be aware of:

1. Data is often generated and made available at different levels:

* **Raw abundance** of sequencing reads for each feature of interest (e.g. integer counts of the number of reads mapped to a given gene)
* **Normalised abundance **, further processed to account for sequence depth and/or gene length such as TPM and FPKM (output using DESeq2, edgeR)
* **Normalised and transformed** adundance with normalisation and further transformation applied such as log2 or variance stabilising transformation (vst) applied

2. The features of interest, say transcripts, may also have been filtered to remove transcripts with read counts below a particular threshold.

3. Transcripts may have been summarised to genes.

4. The data may be given as a matrix of features of interest and samples. Alternatively there may be a separate file for each sample that need to be combined.

**It is important to carefully read the protocols to understand determine how the data has been processed and if is appropriate for your task.**


:::::::::::::::::::::::::::::::::::::: 


::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1:

What type of data do you need to look for in a gene expression data repository as input to a machine learning modelling analysis?

:::::::::::::::::::::::: solution 

The starting point for a machine learning model is processed RNA-Seq data, in the form of count data. This could be raw count data, or TPM or FPKM pre-processed data.

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::


::::::::::::::::::::::::::::::::::::: keypoints 

- Sourcing and appropriate RNA-Seq dataset and preparing it for a machine learning analysis requires you to consider the dataset quality, size, format, noise content, data distribution and scale
- The two main repositories for sourcing RNA-Seq datasets are ArrayExpress and GEO
- Processed data, in the form of raw counts, or further processed counts data is the starting point for machine learning analysis

::::::::::::::::::::::::::::::::::::::::::::::::
