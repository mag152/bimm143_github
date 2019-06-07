Class 13
--------

Focus on te exican Ancestry in LA for Q5.

``` r
#Read CSV file from ENSEMBLE
#mxl <- read.csv("___put csv file here____")
#head(mxl)
```

How many of each genotype are there?

``` r
#table(mxl$Genotype..forward.strand.)
```

Proportion or percent of total for each genotype

``` r
#(table(mxl$Genotype..forward.strand.) / nrow(mxl)) * 100
```

Quality Scores in FASTQ files
-----------------------------

The fourth line of a FASTQ sequence format file encodes the quality
score that tells us how good the sequence at a given position is
(i.e. how likely it is to be correct based on the instrument)

``` r
#library(seqinr)
#library(gtools)

#s2c("DDDDCDEDCDDDDBBDDDCC@")
#acs(s2c("DDDDCDEDCDDDDBBDDDCC@")) - 33
#phred
```

\#Population Scale analysis Readx RNA-Seq count data with genotype
information results table

``` r
expr <- read.table("rs8067378_ENSG00000172057.6.txt")
head(expr)
```

    ##    sample geno      exp
    ## 1 HG00367  A/G 28.96038
    ## 2 NA20768  A/G 20.24449
    ## 3 HG00361  A/A 31.32628
    ## 4 HG00135  A/A 34.11169
    ## 5 NA18870  G/G 18.25141
    ## 6 NA11993  A/A 32.89721

``` r
summary(expr)
```

    ##      sample     geno          exp        
    ##  HG00096:  1   A/A:108   Min.   : 6.675  
    ##  HG00097:  1   A/G:233   1st Qu.:20.004  
    ##  HG00099:  1   G/G:121   Median :25.116  
    ##  HG00100:  1             Mean   :25.640  
    ##  HG00101:  1             3rd Qu.:30.779  
    ##  HG00102:  1             Max.   :51.518  
    ##  (Other):456

``` r
inds <- expr$geno == "G/G"
summary(expr[inds,"exp"])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   6.675  16.903  20.074  20.594  24.457  33.956

``` r
inds <- expr$geno == "A/G"
summary(expr[inds,"exp"])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   7.075  20.626  25.065  25.397  30.552  48.034

``` r
inds <- expr$geno == "A/A"
summary(expr[inds,"exp"])
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   11.40   27.02   31.25   31.82   35.92   51.52

``` r
boxplot(exp ~ geno , data=expr)
```

![](class13_files/figure-markdown_github/unnamed-chunk-10-1.png)

////////////////////////////////////////////////

Overview: The purpose of this lab session is to introduce a set of tools
used in highthroughput sequencing and the process of investigating
interesting gene variance in Genomics. High-throughput sequencing is now
routinely applied to gain insight into a wide range of important topics
in biology and medicine

Section 1: Identify genetic variants of interest There are a number of
gene variants associated with childhood asthma. A study from Verlaan et
al. (2009) shows that 4 candidate SNPs demonstrate significant evidence
for association. You want to find what they are by visiting OMIM
(<a href="http://www.omim.org" class="uri">http://www.omim.org</a>) and
locating the Verlaan et al. paper description.

``` r
#Q1: What are those 4 candidate SNPs?
#[HINT, you will may want to check the first few links of search result] 
#"In association studies of 4 candidate SNPs (rs12936231, rs8067378, rs9303277, and rs7216389)"
```

``` r
#Q2: What three genes do these variants overlap or effect?
#[HINT, you can find the information from the ENSEMBLE page as shown in the image below with red rectangles]
#ZPBP2, GSDMB, and ORMDL3 (???????????????????????????????????)
```

Now, you want to know the location of SNPs and genes in the genome. You
can find the coordinates for the SNP itself on the Ensemble page along
with overlapping genes or whether it is intergenic (i.e. between genes).
However, to explore the surrounding regions and neighboring SNPs you
will need to visit the linked Ensemble genome browser by clicking on the
Location tab (highlighted with a yellow rectangle above).

``` r
#Q3: What is the location of rs8067378 and what are the different alleles for rs8067378?
#[HINT, alleles and location are listed at the top of the the Ensemble page. You may search in a genome browser to find this information]
#
```

Section 2
