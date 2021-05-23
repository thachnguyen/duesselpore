\pagebreak

# Customise the tutorial template

The **`R`** code to prepare your own report is included in the distributed **`Rmarkdown`** file. The **`Rmarkdown`** file can be loaded, viewed, and edited in the **`RStudio`** software. Within your **`conda`** environment (and within your tutorial folder), simply type

\fontsize{8}{12}
```
rstudio Nanopore_cDNA_Tutorial.Rmd
```
\fontsize{10}{14}




**Final thoughts.** Behind this **`Rmarkdown`** file is a modest amount of **`R code`** - please explore the **`Rmarkdown`** template; modify it and run with your own samples.

To extract the whole set of **`R code`** from the **`Rmarkdown`**, use the **`purl`** command - this will extract the R code into its own file.

```
knitr::purl("Nanopore_cDNA_Tutorial.Rmd", quiet=TRUE)
```


# Further explore the expression data

Running the tutorial Rmarkdown script saves the R data objects, presentation methods, and results in a sessionFile. This saved session can be opened directly and it is possible to further explore the data in a dynamic fashion.

To load the data

```
rstudio
```

Once the R console has loaded re-load the saved session with

```
library(session)
restore.session("Analysis/Results/NanoporeDESeq2.Rdata")
```

There are a number of R objects that could be of immediate interest for further exploration of the data

* **`geneCounts`** - a **`data.frame`** containing the number of sequence reads mapped to each of the annotated genes. Biological samples correspond to column names; the annotated gene ids correspond to the row names.
* **`deSeqObj`** - a **`DESeqDataSet`** containing the experimental design and the quantitative data for the genes filtered by the **`readCountMinThreshold`** - count data can be extracted using **`counts(deSeqObj)`**
* **`deSeqRes`** - a **`DeSeqResults`** containing the results from the statistical analysis of differential expression. This can be summarised with **`summary(deSeqRes)`**
* **`deData`** - **`data.frame`** describing the genes with mapped reads, mean expression, log2FoldChange and statistical validation of the observation both with and without false discovery correction

Analysis hints include

1. **List gene ids** for genes with a gene expression profile
```
rownames(deSeqObj)
```
2. **Get scaled expression** for a gene of interest 
```
plotCounts(deSeqObj, gene="ENSG00000108821", intgroup="group", returnData=TRUE)
```
3. **Get expression details** for a gene of interest
```
deData["ENSG00000108821",]
```
4. **Show plot of gene expression** for a gene of interest
```
plotExpressionForGene("ENSG00000108821")
```


\pagebreak

# Glossary of terms

* __knit__ is the command to render an Rmarkdown file. The knitr package is used to embed code, the results of R analyses and their figures within the typeset text from the document. 

* __L50__  describes the number of sequences (or contigs) that are longer than, or equal to, the N50 length and therefore include half the bases of the assembly

* __N50__  describes the length (read length, contig length etc) where half the bases of the sequence collection are contained within reads/contigs of this length or longer

* __Rmarkdown__ is an extension to markdown. Functional R code can be embedded in a plain-text document and subsequently rendered to other formats including the PDF format of this report.

* __QV__  the quality value, -log10(p) that any given base is incorrect. QV may be either at the individual base level, or may be averaged across whole sequences


\pagebreak



