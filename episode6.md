---
title: "Data Readiness: Distribution and Scale"
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- Do I need to transform or rescale RNA-Seq data before inputting into a machine learning algorithm?
- What are the most appropriate transformations and how do these depend on the particular machine learning algorithm being employed?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Describe the main characteristics of the distribution of RNA-Seq count data and how this impacts the performance of machine learning algorithms
- Recall the main types of data transformation commonly applied to RNA-Seq count data to account for the skewedness, zero-inflation and heteroskedasticity
- Demonstrate how to apply these transformations in R
- Recall the main types of data rescaling commonly applied to data prior to inputing into machine learning algorithms, and which rescaling are most appropriate / necessary for a range of commonly used algorithms

::::::::::::::::::::::::::::::::::::::::::::::::


:::::::::::::::::::::::::::::::::::::: callout

## Under Development

This episode is still actively being developed

:::::::::::::::::::::::::::::::::::::: 


## RNA-Seq Counts Data: A Skewed Distribution




## Variance Stabilising and rlog Transformation




## Standardisation and Min-Max Scaling



::::::::::::::::::::::::::::::::::::: keypoints 

- RNA-Seq read count data is heavily skewed with a large percentage of zero values. The distribution is heteroskedastic, meaning the variance depends on the mean.
- Standard transformations such as the variance stabilising transformation and rlog transformation are designed to make the distribution of RNA-Seq data more Gaussian, and therefore more appropriately distributed for many machine learning models. These transformations are improvements to a simple log2 transformation in particular for low count values.
- Many machine learning algorithms require predictor variables to be on the same scale. Standardisation (z-score) and min-max scaling are two common techniques to rescale variables to the same scale. Tree based models are less sensitive to scaling. Standardisation is appropriate for linear models such as logistic regression and svm. min-max scaling is most appropirate for neural network based models.

:::::::::::::::::::::::::::::::::::::
