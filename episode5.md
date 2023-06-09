---
title: "Data Readiness: Technical Noise"
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do technical artefacts in RNA-Seq data impact machine learning algorithms?
- How can technical artefacts such as low count genes and outlier read counts be effectively removed from RNA-Seq data prior to analysis?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Describe the main types of technical noise in RNA-Seq data, namely low count genes and influential outlier read counts, and recall why these exist in raw datasets
- Demonstrate how to remove these elements from the data standard R libraries to further improve machine learning readiness of RNA-Seq count data

::::::::::::::::::::::::::::::::::::::::::::::::

## Technial Artefacts in RNA-Seq Data


## Low Count Filtering


## Outlier Read Count Filtering


::::::::::::::::::::::::::::::::::::: keypoints 

- RNA-Seq read counts contain two main sources of technical 'noise' and are unlikely to represent true biological signals from the samples: Low count genes/transcripts and influential outlier read counts.
- Filtering out low read count and influential outlier genes / transcripts removes potentially biasing variables without negatively impacting the performance of downstream machine learning analysis.

:::::::::::::::::::::::::::::::::::::
