Class 5 R graphics
==================

2A. Line plot
=============

``` r
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
```

above to use console: command+return
====================================

default is false for header -&gt; need to change to TRUE if 1st line is header
==============================================================================

``` r
weight$Age
```

    ##  [1] 0 1 2 3 4 5 6 7 8 9

``` r
weight$Weight
```

    ##  [1] 3.6 4.4 5.2 6.0 6.6 7.2 7.8 8.4 8.8 9.2

``` r
plot(weight, pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight (kg)", main="Baby Weight by Month")
```

![](class05_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
plot(weight$Age, weight$Weight, type = "b", pch=15, cex=1.5, lwd=2, ylim=c(2,10), xlab="Age (months)", ylab="Weight (kg)", main="Baby Weight by Month")
```

![](class05_files/figure-markdown_github/unnamed-chunk-2-2.png)

``` r
# Above is more concise -> favorable if there are more variables
```

\#2B. Line plot

``` r
feat <- read.table("bimm143_05_rstats/feature_counts.txt", 
                   header= TRUE, sep = "\t")
```

sep=" separates info by tab
===========================

``` r
barplot(feat$Count)
```

![](class05_files/figure-markdown_github/unnamed-chunk-4-1.png)

I need to argue w this plot to make it nicer
============================================

``` r
old.par <- par(mar = c(0, 0, 0, 0))
par(mar=c(4,11,1,1))
barplot(feat$Count, horiz = TRUE, xlab="", names.arg = feat$Feature,
        main="Number of features in the mouse GRCm38 genome", las=1, xlim = c(0,80000))
```

![](class05_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
par(mar=old.par)
```

\#?? use par (the parameter) in the mouse plot

Section 3
=========

``` r
counts <- read.table("bimm143_05_rstats/male_female_counts.txt",
                     sep="\t", header=TRUE)
```

use read.delim() to get same result as above
============================================

``` r
counts <-read.delim("bimm143_05_rstats/male_female_counts.txt")

barplot(counts$Count, names.arg = counts$Sample, las=2,
        col=rainbow(10))
```

![](class05_files/figure-markdown_github/unnamed-chunk-7-1.png)
