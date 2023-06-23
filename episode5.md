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

:::::::::::::::::::::::::::::::::::::: callout

## Under Development

This episode is still actively being developed

:::::::::::::::::::::::::::::::::::::: 

## Technical Artefacts in RNA-Seq Data

Machine learning classification algorithms are highly sensitive to any feature data characteristic, regardless of scale, that may differ between experimental groups, and will exploit these data characteristic differences when training a model. Given this, it is important to make sure that data input into a machine learning model reflects true biological signal, and not technical artefacts or noise that stems from the experimental process used to generate the data.

There are two important sources of noise inherent in RNA-Seq data that may negatively impact machine learning modelling performance, namely low read counts, and influential outlier read counts.

<br>

## Low read counts

Genes with consistently low read count values across all samples in a dataset may be technical or biological stochastic artefacts such as the detection of a transcript from a gene that is not uniformly active in a heterogeneous cell population or as the result of a transcriptional error. Below some count threshold, genes are unlikely to be representative of true biological differences related to the condition of interest.

We'll briefly investigate the distribution of read counts in the IBD dataset to illustrate this point. Import the cleaned up counts matrix and sample information text files that we prepared in Episode 4.

<br>


```r
require(tidyverse)
```

```{.output}
Loading required package: tidyverse
```

```{.output}
── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.1.2     ✔ readr     2.1.4
✔ forcats   1.0.0     ✔ stringr   1.5.0
✔ ggplot2   3.4.2     ✔ tibble    3.2.1
✔ lubridate 1.9.2     ✔ tidyr     1.3.0
✔ purrr     1.0.1     
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```


```r
samp.info.ibd.sel <- read.table(file="./data/ibd.sample.info.txt", sep="\t", header=T, fill=T, check.names=F)
```

```{.warning}
Warning in file(file, "rt"): cannot open file './data/ibd.sample.info.txt': No
such file or directory
```

```{.error}
Error in file(file, "rt"): cannot open the connection
```

```r
counts.mat.ibd <- read.table(file="./data/counts.mat.ibd.txt", sep='\t', header=T, fill=T, check.names=F)
```

```{.warning}
Warning in file(file, "rt"): cannot open file './data/counts.mat.ibd.txt': No
such file or directory
```

```{.error}
Error in file(file, "rt"): cannot open the connection
```

Run the following code to view the histogram of the maximum count for each gene in the sample. You'll see in the resulting that over 800 genes have no counts for any gene (maximum = 0), and that there are hundreds of genes where the maximum count for the gene across all the sample is below 10. This compares to a median maximum count value of around 250 for the dataset.


```r
data.frame(max_count = apply(counts.mat.ibd, 1, max, na.rm=TRUE)) %>% 
  ggplot(aes(x = max_count)) + 
    geom_histogram(bins = 200) + 
    xlab("Max Counts (log10 scale)") + 
    ylab("Frequency") +
    scale_x_log10(n.breaks = 6, labels = scales::comma)
```

```{.error}
Error in eval(expr, envir, enclos): object 'counts.mat.ibd' not found
```


## Low Count Filtering

Filtering out low count genes has been show to increase the classification performance of machine learning classifiers, and to increase the stability of the set of genes selected by a machine learning algorithm in the context of selecting relevant genes. A simple approach is to filter out all genes where the maximum count over all samples is below a given threshold.

However, in order to be able to compare read counts between samples, we must first adjust ('normalise') the counts to control for differences in sequence depth and sample composition between samples. To achieve this, we will normalise using the median-of-ratios method implemented in the R package `DESeq2`. For more information on the rationale for scaling RNA-Seq counts and a comparison of the different metrics used see [RDMBites | RNAseq expression data](https://www.youtube.com/watch?v=tO2H3zuBouw).


```r
# convert the condition variable to a factor as required by DESeq2 functions
samp.info.ibd.sel[c('condition')] <- lapply(samp.info.ibd.sel[c('condition')], factor)
```

```{.error}
Error in eval(expr, envir, enclos): object 'samp.info.ibd.sel' not found
```

```r
# create DESeq Data Set object from the raw counts matrix, with condition as the experimental factor of interest
dds.ibd <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts.mat.ibd,
    colData = data.frame(samp.info.ibd.sel, row.names = 'sampleID'),
    design = ~ condition)
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'ncol': object 'counts.mat.ibd' not found
```

```r
# calculate the normalised count values using the median-of-ratios method
dds.ibd <- dds.ibd %>% DESeq2::estimateSizeFactors()
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'estimateSizeFactors': object 'dds.ibd' not found
```

```r
# extract the normalised counts
counts.ibd.norm <- DESeq2::counts(dds.ibd, normalized = TRUE)
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'counts': object 'dds.ibd' not found
```

There are a number of methodologies to calculate the appropriate threshold. One widely used approach calculates the threshold that maximises the similarity between the samples, calculating the mean pairwise Jaccard index over all samples in each conditionn of interest, and setting the threshold that maximises this value. The motivation for this measure is that true biological signal will be consistent between samples in the same condition, however low count noise will be randomly distributed. By maximising the similarity, the threshold is likely to be set above the background noise level. The code for this filter is implemented in the R package `HTSFilter` and described in detail in the article by [Rau et al., 2013](10.1093/bioinformatics/btt350).

Calculating the Jaccard index between all pairs of samples in a dataset does not however scale well with the dimensionality of the data. Here we use an alternative approach that achieves a very similar result using the Multiset Jaccard index, which is faster to compute. Don't worry too much about the code, the main point is that we find a filter threshold, and filter out all the genes where the max count value is below the threshold, on the basis that these are likely technical noise.


```r
# create vector of the class labels
condition <- samp.info.ibd.sel$condition
```

```{.error}
Error in eval(expr, envir, enclos): object 'samp.info.ibd.sel' not found
```

```r
# create a sequence of thresholds to test
t.seq <- seq(1, 25, by = 1)

# Function to calculate the Multiset Jaccard Index' over all samples
Jaccard = function(mat){
  row.sums = rowSums(mat)
  inter = sum(row.sums == ncol(mat))
  union = sum(row.sums > 0)
  return(inter/union)
}

# Calculate vector of the minimum Multiset Jaccard Index over the three classes, for each value of t
ms.jac = sapply(t.seq, function(t){
  filt.counts = counts.ibd.norm >=t
  group.jac = sapply(unique(condition), function(cond){
    Jaccard(filt.counts[,condition==cond])
  })
  return(min(group.jac))
})
```

```{.error}
Error in FUN(X[[i]], ...): object 'counts.ibd.norm' not found
```

```r
# calculate the value of the threshold that maximises the multiset Jaccard index

(t.hold <- which.max(ms.jac))
```

```{.error}
Error in eval(expr, envir, enclos): object 'ms.jac' not found
```

```r
# plot the threshold value against the value of the Multiset Jaccard index to visualise

ggplot(data=data.frame(t = t.seq, jacc = ms.jac)) +
            geom_line(aes(x=t, y=jacc)) +
            geom_hline(yintercept = max(ms.jac), lty=2,col='gray') +
            geom_point(aes(x=t.star, y=max(ms.jac)), col="red", cex=6, pch=1) +
            xlab("Low Count Threshold") + 
            ylab("Multiset Jaccard Index")
```

```{.error}
Error in eval(expr, envir, enclos): object 'ms.jac' not found
```


Having determine a threshold, we then filter the raw counts matrix on the rows (genes) that meet the threshold criterion based on the normalised counts. As you can see, over 4,000 genes are removed from the dataset, approximately 20% of the genes, all of which are likely to contain very little biologically meaningful information but which might easily bias a machine learning classifier.


```r
counts.mat.ibd.filtered <-  counts.mat.ibd[which(apply(counts.ibd.norm, 1, function(x){sum(x > t.hold) >= 1})),]
```

```{.error}
Error in eval(expr, envir, enclos): object 'counts.mat.ibd' not found
```

```r
sprintf("Genes filtered: %s; Genes remaining: %s", nrow(counts.mat.ibd)-nrow(counts.mat.ibd.filtered), nrow(counts.mat.ibd.filtered))
```

```{.error}
Error in eval(expr, envir, enclos): object 'counts.mat.ibd' not found
```


## Outlier Reading Counts

A second source of technical noise is outlier read counts. Relatively high read counts occurring in only
a very small number of samples relative to the size of each patient group are unlikely to be representative
of the general population and a result of biological heterogeneity and technical effects.

These very large count values may bias machine learning algorithms if they help to discriminate between examples in a training set, leading to model overfitting. Such influential outliers may be the result of natural variation between individuals or they may have been introduced during sample preparation. cDNA libraries require PCR amplification prior
to sequencing to achieve sufficient sequence depth. This clonal amplification by PCR is stochastic in nature;
different fragments may be amplified with different probabilities. This leaves the possibility of outlier read
counts having resulted from bias introduced by amplification, rather than biological differences between
samples. The presence of these outlier read counts for a particular gene may therefore inflat the observed association between a particular gene and the condition of interest.

Run the following code to view the top 10 values of read counts in the raw counts matrix. Compare this with the histogram above, and the mean read count value which is approximately 500. The largest read count values range from 2MM to over 3MM counts for a gene in a particular sample. These require investigation.


```r
tail(sort(as.matrix(counts.mat.ibd)),10)
```

```{.error}
Error in eval(expr, envir, enclos): object 'counts.mat.ibd' not found
```

```r
sprintf("The mean read count value: %f", mean(as.matrix(counts.mat.ibd)))
```

```{.error}
Error in eval(expr, envir, enclos): object 'counts.mat.ibd' not found
```

## Outlier Read Count Filtering

As with low counts, multiple methods exist to identify influential outlier read counts. A summary is provided by [Parkinson et al.](https://doi.org/10.3389/fgene.2023.1158352). Here we adopt an approach that is built into DESeq2 that identifies potential influential outiers based on their Cooks distance.

Run the following code to create DESeq Data Set object from the filtered raw counts matrix, with condition as the experimental factor of interest


```r
dds.ibd.filt <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts.mat.ibd.filtered,
    colData = data.frame(samp.info.ibd.sel, row.names = 'sampleID'),
    design = ~ condition)
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'ncol': object 'counts.mat.ibd.filtered' not found
```

Run DESeq2 differential expression analysis to calculate Cook's distances matrix


```r
deseq.ibd <-  DESeq2::DESeq(dds.ibd.filt)
```

```{.error}
Error in eval(expr, envir, enclos): object 'dds.ibd.filt' not found
```

```r
cooks.mat <- SummarizedExperiment::assays(deseq.ibd)[["cooks"]]
```

```{.error}
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'assays': object 'deseq.ibd' not found
```


Calculate the cooks outlier threshold based on the F-distribution


```r
cooks.quantile <- 0.95
m <- ncol(deseq.ibd)     # number of samples
```

```{.error}
Error in eval(expr, envir, enclos): object 'deseq.ibd' not found
```

```r
p <- 3                   # number of model parameters (in the three condition case)

h.threshold <- stats::qf(cooks.quantile, p, m - p)
```

```{.error}
Error in eval(expr, envir, enclos): object 'm' not found
```

Filter the counts matrix to eliminate all genes where the cooks distance is over the outlier threshold 


```r
counts.mat.ibd.final <-  counts.mat.ibd.filtered[which(apply(cooks.mat, 1, function(x){(max(x) < h_threshold) >= 1})),]
```

```{.error}
Error in eval(expr, envir, enclos): object 'counts.mat.ibd.filtered' not found
```

```r
sprintf("Genes filtered: %s; Genes remaining: %s", nrow(counts.mat.ibd.filtered)-nrow(counts.mat.ibd.final), nrow(counts.mat.ibd.final))
```

```{.error}
Error in eval(expr, envir, enclos): object 'counts.mat.ibd.filtered' not found
```


::::::::::::::::::::::::::::::::::::: keypoints 

- RNA-Seq read counts contain two main sources of technical 'noise' and are unlikely to represent true biological signals from the samples: Low count genes/transcripts and influential outlier read counts.
- Filtering out low read count and influential outlier genes / transcripts removes potentially biasing variables without negatively impacting the performance of downstream machine learning analysis.

:::::::::::::::::::::::::::::::::::::
