---
title: "Data Readiness: Data Format and Integrity"
teaching: 20
exercises: 30
---
:::::::::::::::::::::::::::::::::::::: questions 

- What format do I require my data in to build a supervised machine learning classifier?
- What are the potential data format and data integrity issues I will encounter with RNA-Seq datasets, including those downloaded from public repositories, that I'll need to address before beginning any analysis?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Describe the format of a dataset required as input to a simple supervised machine learning analysis
- Recall a checklist of common formatting and data integrity issues encountered with sample information (metadata) in RNA-Seq datasets
- Recall a checklist of common formatting and data integrity issues encountered with counts matrix RNA-Seq datasets
- Recall why formatting and data integrity issues present a problem for downstream ML analysis
- Apply the checklist to an unseen dataset downloaded from either the Array Express or GEO platform to identify potential issues
- Construct a reusable R code pipeline implement required changes to a new dataset

::::::::::::::::::::::::::::::::::::::::::::::::


## Required format for machine learning libraries

Data must be correctly formatted for use in a machine learning / AI pipeline. The *garbage in garbage out* principle applies; a machine learning model is only as good as the input data. This episode first discusses the kind of data required as input to a machine learing classification model.  We'll then illustrates the process of formatting an RNA-Seq dataset downloaded from a public repository. On the way, we'll define a checklist of data integrity issues to watch out for.

For a supervised machine learning classification task, we require:

1. **A matrix of values of all of the predictor variables** to be included in a machine learning model for each of the samples. Example predictor variables would be gene abundance in the form of read counts. Additional predictor variables may be obtained from the sample information including demographic data (e.g. `sex`, `age`) and clinical and laboratory measures.

2. **A vector of values for the target variable** that we are trying to predict for each sample. Example target variables for a classification model would be disease state, for example 'infected with tuburculosis' compared with 'healthy control'. The analogy for a regression task would be a continuous variable such as a disease progression or severity score.

We will follow a series of steps to reformat and clean up the sample information and counts matrix file for the IBD dataset collected form ArrayExpress in Episode 2 by applying a checklist of data clean up items. If required, the data can be downloaded again by running this code.

```r

download.file(url = "https://zenodo.org/record/8125141/files/E-MTAB-11349.sdrf.txt",
              destfile = "data/E-MTAB-11349.sdrf.txt")
              
download.file(url = "https://zenodo.org/record/8125141/files/E-MTAB-11349.counts.matrix.csv",
              destfile = "data/E-MTAB-11349.counts.matrix.csv")

```


We'll start by reading in the dataset from our `data` subfolder.



```r
samp.info.ibd <- read.table(file="./data/E-MTAB-11349.sdrf.txt", sep="\t", header=T, fill=T, check.names=F)

raw.counts.ibd <- read.table(file="./data/E-MTAB-11349.counts.matrix.csv", sep="," , header=T, fill=T, check.names=F)
```


## IDB dataset: Sample Information Clean Up

In real RNA-Seq experiments, sample information (metadata) has likely been collected from a range of sources, for example medical records, clinical and laboratory readings and then manually compiled in a matrix format, often a spreadsheet, before being submitted to a public respository. Each of these steps has the potential to introduce errors and inconsistency in the data that must be investigated prior to any downstream machine learning analysis. Whether you are working with sample information received directly from co-workers or downloaded from a public repsository, the following guidelines will help identify and resolve issues that may impact downstream analysis.


:::::::::::::::::::::::::::::::::::::: checklist

### Sample Information Readiness Checklist

Item  | Check For... | Rationale
- | ------ | -------
1 |  **Unique identifier** matching between counts matrix and sample information | Target variables must be matched to predictor variables for each sample. Target variables are the thing we are trying to predict such as disease state. Predictor variables are the things we are going to predict based on, such as the expression levels of particular genes.
2 | **Appropriate target variable** present in the data | The sample information needs to contain the target variable you are trying to predict with the machine learning model
3 |  **Unique variable names**, without spaces or special characters | Predictor variables may be drawn from the sample information. ML algorithms require unique variables. Spaces etc. may cause errors with bioinformatics and machine learning libraries used downstream
4 |  **Machine readable and consistent encoding**  of categorical variables | Inconsistent spellings between samples will be treated as separate values. Application specific formatting such as colour coding and notes in spreadsheets will be lost when converting to plain text formats.
5 |  **Numerical class values** formatted as a factor (for classification) | Some machine learning algorithms/libraries (e.g. SVM) require classes defined by a number (e.g. -1 and +1 )
6 |  **Balance of classes** in the target variable (for classification) | Imbalanced datasets may negatively impact results of machine learning algorithms


:::::::::::::::::::::::::::::::::::::: 

<br>

We'll now apply these steps sequentially to the sample information for the IBD dataset, contained in our variable `samp.info.ibd`. First let's take a look at the data:


```r
dplyr::glimpse(samp.info.ibd)
```

```{.output}
Rows: 590
Columns: 32
$ `Source Name`                             <chr> "Sample 1", "Sample 2", "Sam…
$ `Characteristics[organism]`               <chr> "Homo sapiens", "Homo sapien…
$ `Characteristics[age]`                    <int> 34, 49, 27, 9, 34, 60, 54, 2…
$ `Unit[time unit]`                         <chr> "year", "year", "year", "yea…
$ `Term Source REF`                         <chr> "EFO", "EFO", "EFO", "EFO", …
$ `Term Accession Number`                   <chr> "UO_0000036", "UO_0000036", …
$ `Characteristics[developmental stage]`    <chr> "adult", "adult", "adult", "…
$ `Characteristics[sex]`                    <chr> "male", "female", "male", "m…
$ `Characteristics[individual]`             <int> 1, 2, 3, 4, 5, 6, 7, 8, 9, 1…
$ `Characteristics[disease]`                <chr> "Crohns disease\tblood\torga…
$ `Characteristics[organism part]`          <chr> "", "blood", "blood", "blood…
$ `Material Type`                           <chr> "", "organism part", "organi…
$ `Protocol REF`                            <chr> "", "P-MTAB-117946", "P-MTAB…
$ `Protocol REF`                            <chr> "", "P-MTAB-117947", "P-MTAB…
$ `Protocol REF`                            <chr> "", "P-MTAB-117948", "P-MTAB…
$ `Extract Name`                            <chr> "", "Sample 2", "Sample 3", …
$ `Material Type`                           <chr> "", "RNA", "RNA", "RNA", "RN…
$ `Comment[LIBRARY_LAYOUT]`                 <chr> "", "PAIRED", "PAIRED", "PAI…
$ `Comment[LIBRARY_SELECTION]`              <chr> "", "other", "other", "other…
$ `Comment[LIBRARY_SOURCE]`                 <chr> "", "TRANSCRIPTOMIC", "TRANS…
$ `Comment[LIBRARY_STRATEGY]`               <chr> "", "Targeted-Capture", "Tar…
$ `Protocol REF`                            <chr> "", "P-MTAB-117950", "P-MTAB…
$ Performer                                 <chr> "", "Wellcome Trust Clinical…
$ `Assay Name`                              <chr> "", "Sample 2", "Sample 3", …
$ `Technology Type`                         <chr> "", "sequencing assay", "seq…
$ `Protocol REF`                            <chr> "", "P-MTAB-117949", "P-MTAB…
$ `Derived Array Data File`                 <chr> "", "ArrayExpress-normalized…
$ `Comment [Derived ArrayExpress FTP file]` <chr> "", "ftp://ftp.ebi.ac.uk/pub…
$ `Protocol REF`                            <chr> "", "P-MTAB-117949", "P-MTAB…
$ `Derived Array Data File`                 <chr> "", "ArrayExpress-raw.csv", …
$ `Comment [Derived ArrayExpress FTP file]` <chr> "", "ftp://ftp.ebi.ac.uk/pub…
$ `Factor Value[disease]`                   <chr> "", "normal", "normal", "nor…
```

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1:

On first glance, which of our checklist items do you think could be an issue in the sample information for this dataset?

:::::::::::::::::::::::: solution 

All of the checklist items apply in this case. Let's go through them...

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: 

<br>

1. Let's verify that there is a unique identifier that matches between the sample information and counts matrix, and identify the name of the column in the sample information. We manually find the names of the samples in `raw.counts.ibd` (the counts matrix), and set to a new variable `ibd.samp.names`. We then write a `for loop` that evaluates which columns in the sample information match the sample names. There is one match, the column named `Source Name`. Note that there are other columns such as `Assay Name` in this dataset that contain identifiers for some but not all samples. This is a good illustration of why it is important to check carefully to ensure you have a complete set of unique identifiers.


```r
ibd.samp.names <- colnames(raw.counts.ibd)[3:ncol(raw.counts.ibd)]

lst.colnames <- c()
for(i in seq_along(1:ncol(samp.info.ibd))){
  lst.colnames[i] <- all(samp.info.ibd[,i] == ibd.samp.names)
}

sprintf("The unique IDs that match the counts matrix are in column: %s", colnames(samp.info.ibd)[which(lst.colnames)])
```

```{.output}
[1] "The unique IDs that match the counts matrix are in column: Source Name"
```


<br>

2. We next extract this unique identifier along with the condition of interest, (in this case the `disease`), and any other variables that may be potential predictor variables or used to assess confounding factors in the subsequent analysis of the data (in this case patient `sex` and `age`). We will use functions from the `dplyr` package, part of the `tidyverse` collection of packages throughout this episode.


```r
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
```



```r
samp.info.ibd.sel <- dplyr::select(samp.info.ibd,
                              'Source Name',
                              'Characteristics[age]',
                              'Characteristics[sex]',
                              'Characteristics[disease]'
                              )
```

<br>

3. Let's then rename the variables to something more easy to interpret for humans, avoiding spaces and special characters. We'll save these selected columns to a new variable name `samp.info.ibd.sel`.


```r
samp.info.ibd.sel <- dplyr::rename(samp.info.ibd.sel,
                          'sampleID' = 'Source Name',
                          'age' = 'Characteristics[age]',
                          'sex' = 'Characteristics[sex]',
                          'condition' = 'Characteristics[disease]'
                          )
```

<br>

4. Now check how each of the variables are encoded in our reduced sample information data to identify any errors, data gaps and inconsistent coding of categorical variables. Firstly let's take a look at the data.


```r
dplyr::glimpse(samp.info.ibd.sel)
```

```{.output}
Rows: 590
Columns: 4
$ sampleID  <chr> "Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", …
$ age       <int> 34, 49, 27, 9, 34, 60, 54, 24, 39, 41, 38, 46, 7, 46, 21, 35…
$ sex       <chr> "male", "female", "male", "male", "female", "female", "femal…
$ condition <chr> "Crohns disease\tblood\torganism part\tP-MTAB-117946\tP-MTAB…
```

<br>


::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 2:

What steps do we need to take to ensure that these four data columns are consistently formatted and machine readable? Are you able to write the code to address these issues?

:::::::::::::::::::::::: solution 

* There is clearly an issue with the coding of the `condition` Crohn's disease. We'll fix that first, and then check the consistency of the coding of categorical variables `sex` and `condition`.


```r
samp.info.ibd.sel$condition[agrep("Crohns", samp.info.ibd.sel$condition)] <- "crohns_disease"

unique(c(samp.info.ibd.sel$sex, samp.info.ibd.sel$condition))
```

```{.output}
[1] "male"               "female"             "crohns_disease"    
[4] "normal"             "ulcerative colitis"
```


* Both `sex` and `condition` are consistently encoded over all samples. We should however remove the spaces in the values for `sampleID` and for the `condition` value ulcerative colitis. Note: since we are changing our unique identifier, we'll need to make the same update later in the counts matrix.



```r
samp.info.ibd.sel <- samp.info.ibd.sel %>%
                        dplyr::mutate(sampleID = gsub(" ", "_", sampleID),
                                      condition = gsub(" ", "_", condition))
```

:::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::::::::::::::

<br>

5. As some machine learning classification algorithms such as support vector machines require a numerical input for each class, we'll add in a new numerical column to represent our target variable called `class` where all patients denoted as `normal` are given the value -1, and all patients denoted as either `crohns_disease` or `ulcerative_colitis` are given the vale +1. We'll also now convert the type of categorical columns to factors as this is preferred by downstream libraries used for data normalisation such as DESeq2.



```r
samp.info.ibd.sel <- samp.info.ibd.sel %>%
                        dplyr::mutate(class = if_else(condition == 'normal', -1, 1))
```



```r
samp.info.ibd.sel[c('sex', 'condition', 'class')] <- lapply(samp.info.ibd.sel[c('sex', 'condition', 'class')], factor)
```

<br>

6. Finally, let's check the distribution of the classes by creating a `table` from the `class` column.


```r
table(samp.info.ibd.sel$class)
```

```{.output}

 -1   1 
267 323 
```

<br>

The two classes are approximately equally represented, so let's check everything one last time.


```r
glimpse(samp.info.ibd.sel)
```

```{.output}
Rows: 590
Columns: 5
$ sampleID  <chr> "Sample_1", "Sample_2", "Sample_3", "Sample_4", "Sample_5", …
$ age       <int> 34, 49, 27, 9, 34, 60, 54, 24, 39, 41, 38, 46, 7, 46, 21, 35…
$ sex       <fct> male, female, male, male, female, female, female, male, male…
$ condition <fct> crohns_disease, normal, normal, normal, normal, ulcerative_c…
$ class     <fct> 1, -1, -1, -1, -1, 1, 1, -1, 1, 1, -1, 1, -1, -1, -1, 1, -1,…
```

<br>

And save the output as a file.

```r

write.table(samp.info.ibd.sel, file="./data/ibd.sample.info.txt", sep = '\t', quote = FALSE, row.names = TRUE)

```

## IDB dataset: Counts Matrix Clean Up

:::::::::::::::::::::::::::::::::::::: checklist

### Count Matrix Readiness Checklist

Item  | Check For... | Rationale
- | ------ | -------
1 |  **Extraneous information** that you don't want to include as a predictor variable | Blank rows, alternative gene/ transcript IDs and metadata about the genes/transcripts will be treated as a predictor variable unless you explicitly remove them
2 |  **Unique variable names**, without spaces or special characters | ML algorithms require unique variables. Spaces etc. may cause errors with bioinformatics and ml libraries used downstream (and they must match the sample information)
3 | **Duplicate data** | Duplicate variables may bias results of analysis and should be removed unless intentional
4 | **Missing values** | Missing values may cause errors and bias results. Missing values must be identified and appropriate treatment determined (e.g. drop vs impute)

:::::::::::::::::::::::::::::::::::::: 

<br>

1. Take a look at a sample of columns and rows to see what the downloaded file looks like


```r
raw.counts.ibd[1:10,1:8]
```

```{.output}
            read Sample 1 Sample 2 Sample 3 Sample 4 Sample 5 Sample 6
1   1          *    13961    16595    20722    17696    25703    20848
2   2 ERCC-00002        0        0        0        0        0        0
3   3 ERCC-00003        0        0        0        0        0        0
4   4 ERCC-00004        0        0        0        0        3        0
5   5 ERCC-00009        0        0        0        0        0        0
6   6 ERCC-00012        0        0        0        0        0        0
7   7 ERCC-00013        0        0        0        0        0        0
8   8 ERCC-00014        0        0        0        0        0        0
9   9 ERCC-00016        0        0        0        0        0        0
10 10 ERCC-00017        2        0        0        0        1        0
```

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 3:

Can you see the irrelevant information that we need to remove from the counts matrix right away?

:::::::::::::::::::::::: solution 

* Remove the first row, which is in this case totals of counts for each sample
* Remove the first columns, which just numbers the transcript IDs
* Move the column named `read` that contains to the transcript IDs to the row names


```r
counts.mat.ibd <- raw.counts.ibd[-1,-1]

rownames(counts.mat.ibd) <- NULL

counts.mat.ibd <-  counts.mat.ibd %>% column_to_rownames('read')

counts.mat.ibd[1:10,1:6]
```

```{.output}
           Sample 1 Sample 2 Sample 3 Sample 4 Sample 5 Sample 6
ERCC-00002        0        0        0        0        0        0
ERCC-00003        0        0        0        0        0        0
ERCC-00004        0        0        0        0        3        0
ERCC-00009        0        0        0        0        0        0
ERCC-00012        0        0        0        0        0        0
ERCC-00013        0        0        0        0        0        0
ERCC-00014        0        0        0        0        0        0
ERCC-00016        0        0        0        0        0        0
ERCC-00017        2        0        0        0        1        0
ERCC-00019        0        0        0        0        0        0
```

::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::

<br>

2. Let's check for duplicate sampleIDs and transcript IDs. Provided there are no duplicate row or column names, the following should return `interger(0)`.


```r
which(duplicated(rownames(counts.mat.ibd)))
```

```{.output}
integer(0)
```

```r
which(duplicated(colnames(counts.mat.ibd)))
```

```{.output}
integer(0)
```

<br>

3. We'll double check the sampleIDs match the sample information file.


```r
if(!identical(colnames(counts.mat.ibd), samp.info.ibd.sel$sampleID)){stop()}
```

```{.error}
Error in eval(expr, envir, enclos): 
```

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 4:

Why did we get an error here?

:::::::::::::::::::::::: solution 

We renamed the samples in the sample information to remove spaces, so we need to do the same here.


```r
colnames(counts.mat.ibd) <- gsub(x = colnames(counts.mat.ibd), pattern = "\ ", replacement = "_")
```

:::::::::::::::::::::::: 
::::::::::::::::::::::::::::::::::::: 

<br>

4. Finally, we'll check that there are no missing values in the counts matrix. Note, there will be many zeros, but we are ensuring that there are no `NA` values or black records. The following should return `FALSE`.



```r
allMissValues <- function(x){all(is.na(x) | x == "")}

allMissValues(counts.mat.ibd)
```

```{.output}
[1] FALSE
```

<br>

Take a final look at the cleaned up matrix.


```r
counts.mat.ibd[1:10,1:6]
```

```{.output}
           Sample_1 Sample_2 Sample_3 Sample_4 Sample_5 Sample_6
ERCC-00002        0        0        0        0        0        0
ERCC-00003        0        0        0        0        0        0
ERCC-00004        0        0        0        0        3        0
ERCC-00009        0        0        0        0        0        0
ERCC-00012        0        0        0        0        0        0
ERCC-00013        0        0        0        0        0        0
ERCC-00014        0        0        0        0        0        0
ERCC-00016        0        0        0        0        0        0
ERCC-00017        2        0        0        0        1        0
ERCC-00019        0        0        0        0        0        0
```

```r
sprintf("There are %i rows, corresponding to the transcript IDs", dim(counts.mat.ibd)[1])
```

```{.output}
[1] "There are 22750 rows, corresponding to the transcript IDs"
```

```r
sprintf("There are %i columns, corresponding to the samples", dim(counts.mat.ibd)[2])
```

```{.output}
[1] "There are 590 columns, corresponding to the samples"
```

<br>

And save the output as a file.

```r

write.table(counts.mat.ibd, file="./data/counts.mat.ibd.txt", sep = '\t', quote = FALSE, row.names = TRUE)

```

:::::::::::::::::::::::::::::::::::::: callout

## Document your reformatting steps

It is crucial that you document your data reformatting so that it is reproducible by other researchers. Create a project folder for your data processing and include:

* Raw input data files collected from ArrayExpress, GEO or similar
* The output counts matrix and target variable vector, stored in a common format such at .txt format
* The full R scripts used to perform all the data clean up steps. Store the script with the input and output data.
* A readme file that explains the steps taken, and how to run the script on the inputs to generate the outputs. This should also include details of all software versions used (e.g. R version, RStudio or Jupyter Notebooks version, and the versions of any code libraries used)

:::::::::::::::::::::::::::::::::::::: 

## The TB Dataset: Another Example

:::::::::::::::::::::::::::::::::::::: discussion

Take a look at [The Tuburculosis (TB) Dataset - E-MTAB-6845](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6845). You can download the sdrf file from zenodo and save it to your `data` directory by running the following code.


```r
download.file(url = "https://zenodo.org/record/8125141/files/E-MTAB-6845.sdrf.txt",
              destfile = "data/E-MTAB-6845.sdrf.txt")
```

Read the sdrf file into R and take a look at the data. Make a list of the potential issues with this data based on the checklist above? What steps would need to be taken to get them ready to train a machine learning classifier?

As a help, the code to read the file in from your data directory is:


```r
samp.info.tb <- read.table(file="./data/E-MTAB-6845.sdrf.txt", sep="\t", header=T, fill=T, check.names=F)
```

:::::::::::::::::::::::::::::::::::::: 

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: instructor

Here are a number of issues that need addressing in the TB dataset:

* One important reformatting point in this dataset is that there are two lines in the sdrf file per sample, as there is a line for each of the read in the paired end reads. The only difference is the file name given for the fastq file. To solve this we need to select the distinct rows with the variables of interest. Check the issue this line of code:

`which(duplicated(samp.info.tb$Comment[sample id])))`

* There are missing values for some variables that may be useful as predictor variables, (e.g. `Characteristics[sex]` missing values are represented by double spaces). It would be better to recode them as `NA` as this is universally understood to represent a missing value. This is easily seem by running this:

`unique(samp.info.tb$`Characteristics[sex]`)`

* A number of variables need renaming and special characters and spaces removed. For example, the variable name `progressor status` and the value `TB progressor` of this variable use a space that could replaced with an underscore

* And perhaps most importantly of all, the dataset is highly unbalanced, with 9 TB progressors and 351 non-progressors when you account for duplicate data. This dataset is very unlikely to perform well as a training dataset for a machine learning classifier!

`table(samp.info.tb$`Factor Value[progressor status (median follow-up 1.9 years)]`)`


::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


:::::::::::::::::::::::::::::::::: challenge

### Challenge Extension

Can you correctly reorder the following code snippets to create a pipeline to reformat the data in line with the sample information checklist? Use the `Characteristics[progressor status (median follow-up 1.9 years)]` as the target variable.

```r

dplyr::rename(
  'sampleID' = 'Comment[sample id]',
  'prog_status' = 'Factor Value[progressor status (median follow-up 1.9 years)]') %>%  

dplyr::select(
            'Comment[sample id]',
            'Factor Value[progressor status (median follow-up 1.9 years)]')  %>%


dplyr::mutate(prog_status = factor(prog_status, levels = c("non_progressor","tb_progressor")), 
             class = factor(class, levels = c(-1, 1))) %>%   

dplyr::mutate(prog_status = dplyr::recode(prog_status,                                
              "non-progressor" = "non_progressor",
              "TB progressor" = "tb_progressor"), 
              class = if_else(prog_status == 'non_progressor', -1, 1)                 
              ) %>% 

dplyr::distinct() %>%

samp.info.tb %>%

dplyr::glimpse()     

```

:::::::::::::::::::::::::::::::::: solution

### Example code to clean up the sample information



```r
samp.info.tb %>% 
  
      dplyr::select(
            'Comment[sample id]',
            'Factor Value[progressor status (median follow-up 1.9 years)]')  %>%             # select columns

      dplyr::rename(
        'sampleID' = 'Comment[sample id]',
        'prog_status' = 'Factor Value[progressor status (median follow-up 1.9 years)]') %>%  # rename them

      dplyr::distinct() %>%                                                                  # remove duplicate rows

      dplyr::mutate(prog_status = dplyr::recode(prog_status,                                 # recode target variable
                    "non-progressor" = "non_progressor",
                    "TB progressor" = "tb_progressor"), 
                    class = if_else(prog_status == 'non_progressor', -1, 1)                 # create numerical class
                    ) %>% 

       dplyr::mutate(prog_status = factor(prog_status, levels = c("non_progressor","tb_progressor")), 
                     class = factor(class, levels = c(-1, 1))) %>%                          # set variables to factors

       dplyr::glimpse()                                                                     # view output
```

```{.output}
Rows: 360
Columns: 3
$ sampleID    <chr> "PR123_S19", "PR096_S13", "PR146_S14", "PR158_S12", "PR095…
$ prog_status <fct> non_progressor, non_progressor, non_progressor, non_progre…
$ class       <fct> -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1…
```
::::::::::::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: keypoints 

* Machine learning algorithms require specific data input formats, and for data to be consistently formatted across variables and samples. A classification task for example requires a matrix of the value of each input variable for each sample, and the class label for each sample.
* Clinical sample information is often messy, with inconsistent formatting as a result of how it is collected. This applies to data downloaded from public repositories.
* You must carefully check all data for formatting and data integrity issues that may negatively impact your downstream ml analysis.
* Document your data reformatting and store the code pipeline used along with the raw and reformatted data to ensure your procedure is reproducible

::::::::::::::::::::::::::::::::::::::::::::::::
