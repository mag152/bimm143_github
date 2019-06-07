1.  Import countData and colData

``` r
#Begin a new Rmarkdown document and use the read.csv() function to read these count data and metadata files.
counts <- read.csv("data/airway_scaledcounts.csv", stringsAsFactors = FALSE)
metadata <-  read.csv("data/airway_metadata.csv", stringsAsFactors = FALSE)
```

Lets perform some exploratory differential gene expression analysis.

``` r
#Now, take a look at each.
View(counts)
View(metadata)
```

1.  Toy differential gene expression Lets perform some exploratory
    differential gene expression analysis.

``` r
#Look at the metadata object again to see which samples are control and which are drug treated
control <- metadata[metadata[,"dex"]=="control",]
control.mean <- rowSums( counts[ ,control$id] )/4 
names(control.mean) <- counts$ensgene
```

``` r
treated <- metadata[metadata[,"dex"]=="treated",]
treated.mean <- rowSums( counts[ ,treated$id] )/4 
names(treated.mean) <- counts$ensgene
```

``` r
#We will combine our meancount data for bookkeeping purposes.
meancounts <- data.frame(control.mean, treated.mean)
```

``` r
#Directly comparing the raw counts is going to be problematic if we just happened to sequence one group at a higher depth than another. Later on we’ll do this analysis properly, normalizing by sequencing depth per sample using a better approach. But for now, colSums() the data to show the sum of the mean counts across all genes for each group. Your answer should look like this:
colSums(meancounts)
```

    ## control.mean treated.mean 
    ##     23005324     22196524

``` r
#A scatter plot showing the mean of the treated samples against the mean of the control samples (plotting both axes on a log scale).
plot.default(log(control.mean),log(treated.mean))
```

![](class14_files/figure-markdown_github/unnamed-chunk-7-1.png)

We can find candidate differentially expressed genes by looking for
genes with a large change between control and dex-treated samples. We
usually look at the log2 of the fold change, because this has better
mathematical properties.

``` r
#Here we calculate log2foldchange, add it to our meancounts data.frame and inspect the results either with the head() or the View() function for example.

meancounts$log2fc <- log2(meancounts[,"treated.mean"]/meancounts[,"control.mean"])
```

There are a couple of “weird” results. Namely, the NaN (“not a number”“)
and -Inf (negative infinity) results.

The NaN is returned when you divide by zero and try to take the log. The
-Inf is returned when you try to take the log of zero. It turns out that
there are a lot of genes with zero expression. Let’s filter our data to
remove these genes. Again inspect your result (and the intermediate
steps) to see if things make sense to you

``` r
zero.vals <- which(meancounts[,1:2]==0, arr.ind=TRUE)

to.rm <- unique(zero.vals[,1])
mycounts <- meancounts[-to.rm,]
```

A common threshold used for calling something differentially expressed
is a log2(FoldChange) of greater than 2 or less than -2. Let’s filter
the dataset both ways to see how many genes are up or down-regulated.

``` r
up.ind <- mycounts$log2fc > 2
down.ind <- mycounts$log2fc < (-2)
```

``` r
#Determine how many up and down regulated genes we have at the greater than 2 fc level.
sum(as.numeric(up.ind))
```

    ## [1] 250

``` r
sum(as.numeric(down.ind))
```

    ## [1] 367

``` r
#250 upregulated
#367 down regulated
```

1.  Adding annotation data Our mycounts result table so far only
    contains the Ensembl gene IDs. However, alternative gene names and
    extra annotation are usually required for informative for
    interpretation.

We can add annotation from a supplied CSV file, such as those available
from ENSEMBLE or UCSC. The annotables\_grch38.csv annotation table links
the unambiguous Ensembl gene ID to other useful annotation like the gene
symbol, full gene name, location, Entrez gene ID, etc.

``` r
#read.csv("annotables_grch38.csv")
#anno <- read.csv("annotables_grch38.csv")
```

Ideally we want this annotation data mapped (or merged) with our
mycounts data. In a previous class on writing R functions we introduced
the merge() function, which is one common way to do this.

In cases where you don’t have a preferred annotation file at hand you
can use other Bioconductor packages for annotation.

Bioconductor’s annotation packages help with mapping various ID schemes
to each other. Here we load the AnnotationDbi package and the annotation
package org.Hs.eg.db.

``` r
library("AnnotationDbi")
```

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
    ##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    ##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
    ##     setdiff, sort, table, tapply, union, unique, unsplit, which,
    ##     which.max, which.min

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: IRanges

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

``` r
library("org.Hs.eg.db")
```

    ## 

This is the organism annotation package (“org”) for Homo sapiens (“Hs”),
organized as an AnnotationDbi database package (“db”), using Entrez Gene
IDs (“eg”) as primary key. To get a list of all available key types,
use:

``` r
columns(org.Hs.eg.db)
```

    ##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
    ##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
    ##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
    ## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
    ## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
    ## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    ## [25] "UNIGENE"      "UNIPROT"

We can use the mapIds() function to add individual columns to our
results table. We provide the row names of our results table as a key,
and specify that keytype=ENSEMBL. The column argument tells the mapIds()
function which information we want, and the multiVals argument tells the
function what to do if there are multiple possible values for a single
input value. Here we ask to just give us back the first one that occurs
in the database

``` r
mycounts$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(mycounts),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
mycounts$entrez<- mapIds(org.Hs.eg.db,
                     keys=row.names(mycounts),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
mycounts$uniprot <- mapIds(org.Hs.eg.db,
                     keys=row.names(mycounts),
                     column="UNIPROT",
                     keytype="ENSEMBL",
                     multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

1.  DESeq2 analysis Let’s do this the right way. DESeq2 is an R package
    for analyzing count-based NGS data like RNA-seq. It is available
    from Bioconductor. Bioconductor is a project to provide tools for
    analyzing high-throughput genomic data including RNA-seq, ChIP-seq
    and arrays. You can explore Bioconductor packages here.

Bioconductor packages usually have great documentation in the form of
vignettes. For a great example, take a look at the DESeq2 vignette for
analyzing count data. This 40+ page manual is packed full of examples on
using DESeq2, importing data, fitting models, creating visualizations,
references, etc.

Just like R packages from CRAN, you only need to install Bioconductor
packages once (instructions here), then load them every time you start a
new R session.

``` r
library(DESeq2)
```

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## Loading required package: BiocParallel

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

``` r
citation("DESeq2")
```

    ## 
    ##   Love, M.I., Huber, W., Anders, S. Moderated estimation of fold
    ##   change and dispersion for RNA-seq data with DESeq2 Genome
    ##   Biology 15(12):550 (2014)
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
    ##     author = {Michael I. Love and Wolfgang Huber and Simon Anders},
    ##     year = {2014},
    ##     journal = {Genome Biology},
    ##     doi = {10.1186/s13059-014-0550-8},
    ##     volume = {15},
    ##     issue = {12},
    ##     pages = {550},
    ##   }

Importing data Bioconductor software packages often define and use
custom class objects for storing data. This helps to ensure that all the
needed data for analysis (and the results) are available. DESeq works on
a particular type of object called a DESeqDataSet. The DESeqDataSet is a
single object that contains input values, intermediate calculations like
how things are normalized, and all results of a differential expression
analysis.

You can construct a DESeqDataSet from (1) a count matrix, (2) a metadata
file, and (3) a formula indicating the design of the experiment.

We have talked about (1) and (2) previously. The third needed item that
has to be specified at the beginning of the analysis is a design
formula. This tells DESeq2 which columns in the sample information table
(colData) specify the experimental design (i.e. which groups the samples
belong to) and how these factors should be used in the analysis.
Essentially, this formula expresses how the counts for each gene depend
on the variables in colData.

Take a look at metadata again. The thing we’re interested in is the dex
column, which tells us which samples are treated with dexamethasone
versus which samples are untreated controls. We’ll specify the design
with a tilde, like this: design=\~dex. (The tilde is the shifted key to
the left of the number 1 key on my keyboard. It looks like a little
squiggly line).

We will use the DESeqDataSetFromMatrix() function to build the required
DESeqDataSet object and call it dds, short for our DESeqDataSet. If you
get a warning about “some variables in design formula are characters,
converting to factors” don’t worry about it. Take a look at the dds
object once you create it.

``` r
dds <- DESeqDataSetFromMatrix(countData=counts, 
                              colData=metadata, 
                              design=~dex, 
                              tidy=TRUE)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

``` r
dds
```

    ## class: DESeqDataSet 
    ## dim: 38694 8 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(38694): ENSG00000000003 ENSG00000000005 ...
    ##   ENSG00000283120 ENSG00000283123
    ## rowData names(0):
    ## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
    ## colData names(4): id dex celltype geo_id

DESeq pipeline Next, let’s run the DESeq pipeline on the dataset, and
reassign the resulting object back to the same variable. Before we
start, dds is a bare-bones DESeqDataSet. The DESeq() function takes a
DESeqDataSet and returns a DESeqDataSet, but with lots of other
information filled in (normalization, dispersion estimates, differential
expression results, etc). Notice how if we try to access these objects
before running the analysis, nothing exists.

``` r
sizeFactors(dds)
```

    ## NULL

``` r
dispersions(dds)
```

    ## NULL

``` r
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

Getting results Since we’ve got a fairly simple design (single factor,
two groups, treated versus control), we can get results out of the
object simply by calling the results() function on the DESeqDataSet that
has been run through the pipeline. The help page for ?results and the
vignette both have extensive documentation about how to pull out the
results for more complicated models (multi-factor experiments, specific
contrasts, interaction terms, time courses, etc.).

``` r
res <- results(dds)
res
```

    ## log2 fold change (MLE): dex treated vs control 
    ## Wald test p-value: dex treated vs control 
    ## DataFrame with 38694 rows and 6 columns
    ##                          baseMean     log2FoldChange             lfcSE
    ##                         <numeric>          <numeric>         <numeric>
    ## ENSG00000000003  747.194195359907  -0.35070302068658 0.168245681332529
    ## ENSG00000000005                 0                 NA                NA
    ## ENSG00000000419  520.134160051965  0.206107766417862 0.101059218008052
    ## ENSG00000000457  322.664843927049 0.0245269479387466 0.145145067649248
    ## ENSG00000000460   87.682625164828  -0.14714204922212 0.257007253994673
    ## ...                           ...                ...               ...
    ## ENSG00000283115                 0                 NA                NA
    ## ENSG00000283116                 0                 NA                NA
    ## ENSG00000283119                 0                 NA                NA
    ## ENSG00000283120 0.974916032393564  -0.66825846051647  1.69456285241871
    ## ENSG00000283123                 0                 NA                NA
    ##                               stat             pvalue              padj
    ##                          <numeric>          <numeric>         <numeric>
    ## ENSG00000000003  -2.08446967499531 0.0371174658432818 0.163034808641677
    ## ENSG00000000005                 NA                 NA                NA
    ## ENSG00000000419   2.03947517584631 0.0414026263001157 0.176031664879167
    ## ENSG00000000457  0.168982303952742  0.865810560623564 0.961694238404392
    ## ENSG00000000460  -0.57252099672319  0.566969065257939 0.815848587637731
    ## ...                            ...                ...               ...
    ## ENSG00000283115                 NA                 NA                NA
    ## ENSG00000283116                 NA                 NA                NA
    ## ENSG00000283119                 NA                 NA                NA
    ## ENSG00000283120 -0.394354484734893  0.693319342566817                NA
    ## ENSG00000283123                 NA                 NA                NA

``` r
#We can summarize some basic tallies using the summary function.
summary(res)
```

    ## 
    ## out of 25258 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 1563, 6.2%
    ## LFC < 0 (down)     : 1188, 4.7%
    ## outliers [1]       : 142, 0.56%
    ## low counts [2]     : 9971, 39%
    ## (mean count < 10)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
#We can order our results table by the smallest p value:
resOrdered <- res[order(res$pvalue),]
```

The results function contains a number of arguments to customize the
results table. By default the argument alpha is set to 0.1. If the
adjusted p value cutoff will be a value other than 0.1, alpha should be
set to that value:

``` r
res05 <- results(dds, alpha=0.05)
summary(res05)
```

    ## 
    ## out of 25258 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 1236, 4.9%
    ## LFC < 0 (down)     : 933, 3.7%
    ## outliers [1]       : 142, 0.56%
    ## low counts [2]     : 9033, 36%
    ## (mean count < 6)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

The more generic way to access the actual subset of the data.frame
passing a threshold like this is with the subset() function, e.g.:

``` r
resSig05 <- subset(as.data.frame(res), padj < 0.05)
nrow(resSig05)
```

    ## [1] 2181

``` r
#Q9. How many are significant with an adjusted p-value < 0.05? How about 0.01? Save this last set of results as resSig01.
resSig01 <- subset(as.data.frame(res), padj < 0.01)
nrow(resSig01)
```

    ## [1] 1437

``` r
#Q10. Using either the previously generated anno object (annotations from the file  annotables_grch38.csv file) or the mapIds() function (from the AnnotationDbi package) add annotation to your res01 results data.frame.
resSig01$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(resSig01),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
ord <- order( resSig01$padj )
#View(res01[ord,])
```

Finally, let’s write out the ordered significant results with
annotations.

``` r
write.csv(resSig01[ord,], "signif01_results.csv")
```
