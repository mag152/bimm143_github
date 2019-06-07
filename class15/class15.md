Section 1. Differential Expression Analysis You can download the count
data and associated metadata from here: GSE37704\_featurecounts.csv and
GSE37704\_metadata.csv. This is similar to our starting point for the
last class where we used DESeq2 for the first time. We will use it again
today!

\#\#RNA Seq Analysis

The data for for hands-on session comes from GEO entry: GSE37704, which
is associated with the following publication:

> Trapnell C, Hendrickson DG, Sauvageau M, Goff L et al. “Differential
> analysis of gene regulation at transcript resolution with RNA-seq”.
> Nat Biotechnol 2013 Jan;31(1):46-53. PMID: 23222703

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

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

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

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
metaFile <- "GSE37704_metadata.csv"
countFile <- "GSE37704_featurecounts.csv"
#I removed /data/ from the .csv file for the data to run

# Import metadata and take a peak
colData = read.csv(metaFile, row.names=1)
head(colData)
```

    ##               condition
    ## SRR493366 control_sirna
    ## SRR493367 control_sirna
    ## SRR493368 control_sirna
    ## SRR493369      hoxa1_kd
    ## SRR493370      hoxa1_kd
    ## SRR493371      hoxa1_kd

``` r
# Import countdata
countData = read.csv(countFile, row.names=1)
head(countData)
```

    ##                 length SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ## ENSG00000186092    918         0         0         0         0         0
    ## ENSG00000279928    718         0         0         0         0         0
    ## ENSG00000279457   1982        23        28        29        29        28
    ## ENSG00000278566    939         0         0         0         0         0
    ## ENSG00000273547    939         0         0         0         0         0
    ## ENSG00000187634   3214       124       123       205       207       212
    ##                 SRR493371
    ## ENSG00000186092         0
    ## ENSG00000279928         0
    ## ENSG00000279457        46
    ## ENSG00000278566         0
    ## ENSG00000273547         0
    ## ENSG00000187634       258

Hmm… remember that we need the countData and colData files to match up
so we will need to remove that odd first column in countData namely
contData$length.

``` r
#remove first col from countData
countData <- as.matrix(countData[,-1])
head(countData)
```

    ##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ## ENSG00000186092         0         0         0         0         0
    ## ENSG00000279928         0         0         0         0         0
    ## ENSG00000279457        23        28        29        29        28
    ## ENSG00000278566         0         0         0         0         0
    ## ENSG00000273547         0         0         0         0         0
    ## ENSG00000187634       124       123       205       207       212
    ##                 SRR493371
    ## ENSG00000186092         0
    ## ENSG00000279928         0
    ## ENSG00000279457        46
    ## ENSG00000278566         0
    ## ENSG00000273547         0
    ## ENSG00000187634       258

``` r
#Don't run above repeatedly bc this will continue to delete the first column every run.
```

Lets remove the rows with zero counts in all experiments (i.e. columns)

``` r
nonzero.rows <- rowSums(countData) != 0
countData <- (countData[nonzero.rows,])
#This is a much better alternative
```

\#Show how many rows we have after removing columns w/ zeros

``` r
nrow(countData)
```

    ## [1] 15975

Running DESeq2 Nice now lets setup the DESeqDataSet object required for
the DESeq() function and then run the DESeq pipeline. This is again
similar to our last days hands-on session.

``` r
dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)
dds = DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
dds
```

    ## class: DESeqDataSet 
    ## dim: 15975 6 
    ## metadata(1): version
    ## assays(4): counts mu H cooks
    ## rownames(15975): ENSG00000279457 ENSG00000187634 ...
    ##   ENSG00000276345 ENSG00000271254
    ## rowData names(22): baseMean baseVar ... deviance maxCooks
    ## colnames(6): SRR493366 SRR493367 ... SRR493370 SRR493371
    ## colData names(2): condition sizeFactor

Next, get results for the HoxA1 knockdown versus control siRNA (remember
that these were labeled as “hoxa1\_kd” and “control\_sirna” in our
original colData metaFile input to DESeq, you can check this above and
by running resultsNames(dds) command).

``` r
res = results(dds, contrast=c("condition", "hoxa1_kd", "control_sirna"))
summary(res)
```

    ## 
    ## out of 15975 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 4349, 27%
    ## LFC < 0 (down)     : 4396, 28%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 1237, 7.7%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

Volcono plot Now we will make a volcano plot, a commonly produced
visualization from this type of data that we introduced last day.
Basically it’s a plot of log2 fold change vs -log adjusted p-value.

``` r
#Volcano plot is a plot of log2 fold change vs -log adjusted p-value.
plot( res$log2FoldChange, -log(res$padj) )


abline(v=c(-2,2), col="green", lty=2, lwd=2)
abline(h=-log(0.01),col="green", lty=2, lwd=2)
```

![](class15_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
# Make a color vector for all genes
mycols <- rep("gray", nrow(res) )

# Color red the genes with absolute fold change above 2
mycols[ abs(res$log2FoldChange) > 2 ] <- "red"

# Color blue those with adjusted p-value less than 0.01
#  and absolute fold change more than 2
inds <- (res$padj < 0.01) & (abs(res$log2FoldChange) > 2 )
mycols[ inds ] <- "blue"

#use "col=mycols" to put the color above into the same plot
plot( res$log2FoldChange, -log(res$padj), col=mycols, xlab="Log2(FoldChange)", ylab="-Log(P-value)" )
```

![](class15_files/figure-markdown_github/unnamed-chunk-10-2.png)

Adding gene annotation Since we mapped and counted against the Ensembl
annotation, our results only have information about Ensembl gene IDs.
However, our pathway analysis downstream will use KEGG pathways, and
genes in KEGG pathways are annotated with Entrez gene IDs. So lets add
them as we did the last day.

``` r
#Q. Use the mapIDs() function multiple times to add SYMBOL, ENTREZID and GENENAME annotation to our results.
library("AnnotationDbi")
library("org.Hs.eg.db")
```

    ## 

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

``` r
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(res), 
                    keytype="ENSEMBL",
                    column="SYMBOL",
                    multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    keytype="ENSEMBL",
                    column="ENTREZID",
                    multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    keytype="ENSEMBL",
                    column="GENENAME",
                    multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
head(res, 10)
```

    ## log2 fold change (MLE): condition hoxa1_kd vs control_sirna 
    ## Wald test p-value: condition hoxa1 kd vs control sirna 
    ## DataFrame with 10 rows and 9 columns
    ##                          baseMean     log2FoldChange              lfcSE
    ##                         <numeric>          <numeric>          <numeric>
    ## ENSG00000279457  29.9135794276176  0.179257083672691  0.324821565250145
    ## ENSG00000187634  183.229649921658   0.42645711840331  0.140265820376891
    ## ENSG00000188976  1651.18807619944 -0.692720464846371 0.0548465415913881
    ## ENSG00000187961  209.637938486147  0.729755610585229  0.131859899969346
    ## ENSG00000187583  47.2551232589398 0.0405765278756319  0.271892808601774
    ## ENSG00000187642  11.9797501642461  0.542810491577362  0.521559849534146
    ## ENSG00000188290  108.922127976716    2.0570638345631  0.196905312993835
    ## ENSG00000187608   350.71686801731  0.257383686481775  0.102726560033547
    ## ENSG00000188157    9128.439421961  0.389908792022771 0.0467163395511497
    ## ENSG00000237330 0.158192358990472  0.785955208142751    4.0804728567969
    ##                              stat               pvalue
    ##                         <numeric>            <numeric>
    ## ENSG00000279457 0.551863246932652    0.581042050747029
    ## ENSG00000187634  3.04034951107426  0.00236303749730955
    ## ENSG00000188976 -12.6301576133497  1.4398954015367e-36
    ## ENSG00000187961  5.53432552849562 3.12428248077692e-08
    ## ENSG00000187583  0.14923722361139    0.881366448669145
    ## ENSG00000187642  1.04074439790984    0.297994191720984
    ## ENSG00000188290   10.446969679419 1.51281875407292e-25
    ## ENSG00000187608  2.50552229528295   0.0122270689409891
    ## ENSG00000188157  8.34630443585717 7.04321148771393e-17
    ## ENSG00000237330 0.192613757210411    0.847261469988086
    ##                                 padj      symbol      entrez
    ##                            <numeric> <character> <character>
    ## ENSG00000279457    0.686554777832896          NA          NA
    ## ENSG00000187634  0.00515718149494272      SAMD11      148398
    ## ENSG00000188976 1.76548905389749e-35       NOC2L       26155
    ## ENSG00000187961 1.13412993107612e-07      KLHL17      339451
    ## ENSG00000187583    0.919030615571379     PLEKHN1       84069
    ## ENSG00000187642     0.40337930975409       PERM1       84808
    ## ENSG00000188290 1.30538189681069e-24        HES4       57801
    ## ENSG00000187608   0.0237452288908021       ISG15        9636
    ## ENSG00000188157 4.21962808560682e-16        AGRN      375790
    ## ENSG00000237330                   NA      RNF223      401934
    ##                                                                     name
    ##                                                              <character>
    ## ENSG00000279457                                                       NA
    ## ENSG00000187634                 sterile alpha motif domain containing 11
    ## ENSG00000188976 NOC2 like nucleolar associated transcriptional repressor
    ## ENSG00000187961                              kelch like family member 17
    ## ENSG00000187583                 pleckstrin homology domain containing N1
    ## ENSG00000187642             PPARGC1 and ESRR induced regulator, muscle 1
    ## ENSG00000188290                   hes family bHLH transcription factor 4
    ## ENSG00000187608                            ISG15 ubiquitin like modifier
    ## ENSG00000188157                                                    agrin
    ## ENSG00000237330                                  ring finger protein 223

``` r
head(as.data.frame(res))
```

    ##                   baseMean log2FoldChange      lfcSE        stat
    ## ENSG00000279457   29.91358     0.17925708 0.32482157   0.5518632
    ## ENSG00000187634  183.22965     0.42645712 0.14026582   3.0403495
    ## ENSG00000188976 1651.18808    -0.69272046 0.05484654 -12.6301576
    ## ENSG00000187961  209.63794     0.72975561 0.13185990   5.5343255
    ## ENSG00000187583   47.25512     0.04057653 0.27189281   0.1492372
    ## ENSG00000187642   11.97975     0.54281049 0.52155985   1.0407444
    ##                       pvalue         padj  symbol entrez
    ## ENSG00000279457 5.810421e-01 6.865548e-01    <NA>   <NA>
    ## ENSG00000187634 2.363037e-03 5.157181e-03  SAMD11 148398
    ## ENSG00000188976 1.439895e-36 1.765489e-35   NOC2L  26155
    ## ENSG00000187961 3.124282e-08 1.134130e-07  KLHL17 339451
    ## ENSG00000187583 8.813664e-01 9.190306e-01 PLEKHN1  84069
    ## ENSG00000187642 2.979942e-01 4.033793e-01   PERM1  84808
    ##                                                                     name
    ## ENSG00000279457                                                     <NA>
    ## ENSG00000187634                 sterile alpha motif domain containing 11
    ## ENSG00000188976 NOC2 like nucleolar associated transcriptional repressor
    ## ENSG00000187961                              kelch like family member 17
    ## ENSG00000187583                 pleckstrin homology domain containing N1
    ## ENSG00000187642             PPARGC1 and ESRR induced regulator, muscle 1

Writeout our ordered and annotated results object

``` r
#Q. Finally for this section let's reorder these results by adjusted p-value and save them to a CSV file in your current project directory.
res = res[order(res$pvalue),]
write.csv(res, file="deseq_results.csv")
```

Great, this is looking good so far. Now lets see how pathway analysis
can help us make further sense out of this ranked list of differentially
expressed genes.

\#\#Section 2. Pathway analysis Here we are going to use the gage
package for pathway analysis. Once we have a list of enriched pathways,
we’re going to use the pathview package to draw pathway diagrams,
shading the molecules in the pathway by their degree of
up/down-regulation.

``` r
#Install packages in your CONSOLE
#BiocManager::install( c("pathview", "gage", "gageData") )
```

``` r
library(pathview)
```

    ## ##############################################################################
    ## Pathview is an open source software package distributed under GNU General
    ## Public License version 3 (GPLv3). Details of GPLv3 is available at
    ## http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    ## formally cite the original Pathview paper (not just mention it) in publications
    ## or products. For details, do citation("pathview") within R.
    ## 
    ## The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    ## license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ## ##############################################################################

``` r
library(gage)
library(gageData)
```

``` r
data(kegg.sets.hs)
data(sigmet.idx.hs)

# Focus on signaling and metabolic pathways only
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
```

``` r
# Examine the first 3 pathways
head(kegg.sets.hs)
```

    ## $`hsa00232 Caffeine metabolism`
    ## [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   
    ## 
    ## $`hsa00983 Drug metabolism - other enzymes`
    ##  [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"  
    ##  [8] "1551"   "1553"   "1576"   "1577"   "1806"   "1807"   "1890"  
    ## [15] "221223" "2990"   "3251"   "3614"   "3615"   "3704"   "51733" 
    ## [22] "54490"  "54575"  "54576"  "54577"  "54578"  "54579"  "54600" 
    ## [29] "54657"  "54658"  "54659"  "54963"  "574537" "64816"  "7083"  
    ## [36] "7084"   "7172"   "7363"   "7364"   "7365"   "7366"   "7367"  
    ## [43] "7371"   "7372"   "7378"   "7498"   "79799"  "83549"  "8824"  
    ## [50] "8833"   "9"      "978"   
    ## 
    ## $`hsa00230 Purine metabolism`
    ##   [1] "100"    "10201"  "10606"  "10621"  "10622"  "10623"  "107"   
    ##   [8] "10714"  "108"    "10846"  "109"    "111"    "11128"  "11164" 
    ##  [15] "112"    "113"    "114"    "115"    "122481" "122622" "124583"
    ##  [22] "132"    "158"    "159"    "1633"   "171568" "1716"   "196883"
    ##  [29] "203"    "204"    "205"    "221823" "2272"   "22978"  "23649" 
    ##  [36] "246721" "25885"  "2618"   "26289"  "270"    "271"    "27115" 
    ##  [43] "272"    "2766"   "2977"   "2982"   "2983"   "2984"   "2986"  
    ##  [50] "2987"   "29922"  "3000"   "30833"  "30834"  "318"    "3251"  
    ##  [57] "353"    "3614"   "3615"   "3704"   "377841" "471"    "4830"  
    ##  [64] "4831"   "4832"   "4833"   "4860"   "4881"   "4882"   "4907"  
    ##  [71] "50484"  "50940"  "51082"  "51251"  "51292"  "5136"   "5137"  
    ##  [78] "5138"   "5139"   "5140"   "5141"   "5142"   "5143"   "5144"  
    ##  [85] "5145"   "5146"   "5147"   "5148"   "5149"   "5150"   "5151"  
    ##  [92] "5152"   "5153"   "5158"   "5167"   "5169"   "51728"  "5198"  
    ##  [99] "5236"   "5313"   "5315"   "53343"  "54107"  "5422"   "5424"  
    ## [106] "5425"   "5426"   "5427"   "5430"   "5431"   "5432"   "5433"  
    ## [113] "5434"   "5435"   "5436"   "5437"   "5438"   "5439"   "5440"  
    ## [120] "5441"   "5471"   "548644" "55276"  "5557"   "5558"   "55703" 
    ## [127] "55811"  "55821"  "5631"   "5634"   "56655"  "56953"  "56985" 
    ## [134] "57804"  "58497"  "6240"   "6241"   "64425"  "646625" "654364"
    ## [141] "661"    "7498"   "8382"   "84172"  "84265"  "84284"  "84618" 
    ## [148] "8622"   "8654"   "87178"  "8833"   "9060"   "9061"   "93034" 
    ## [155] "953"    "9533"   "954"    "955"    "956"    "957"    "9583"  
    ## [162] "9615"  
    ## 
    ## $`hsa04514 Cell adhesion molecules (CAMs)`
    ##   [1] "1000"      "1001"      "100133583" "1002"      "1003"     
    ##   [6] "100506658" "1013"      "10666"     "10686"     "1272"     
    ##  [11] "1364"      "1365"      "1366"      "137075"    "1462"     
    ##  [16] "1493"      "149461"    "214"       "22871"     "23114"    
    ##  [21] "23308"     "23562"     "23705"     "24146"     "257194"   
    ##  [26] "25945"     "26047"     "26285"     "2734"      "29126"    
    ##  [31] "29851"     "3105"      "3106"      "3107"      "3108"     
    ##  [36] "3109"      "3111"      "3112"      "3113"      "3115"     
    ##  [41] "3117"      "3118"      "3119"      "3122"      "3123"     
    ##  [46] "3125"      "3126"      "3127"      "3133"      "3134"     
    ##  [51] "3135"      "3383"      "3384"      "3385"      "3655"     
    ##  [56] "3676"      "3680"      "3683"      "3684"      "3685"     
    ##  [61] "3688"      "3689"      "3695"      "3696"      "3897"     
    ##  [66] "4099"      "4267"      "4359"      "4684"      "4685"     
    ##  [71] "4756"      "4897"      "4950"      "49861"     "5010"     
    ##  [76] "50848"     "51208"     "5133"      "5175"      "53842"    
    ##  [81] "54413"     "57502"     "57555"     "57863"     "5788"     
    ##  [86] "5792"      "5797"      "5817"      "5818"      "5819"     
    ##  [91] "58494"     "6382"      "6383"      "6385"      "6401"     
    ##  [96] "6402"      "6403"      "6404"      "652614"    "6614"     
    ## [101] "6693"      "6900"      "7122"      "7412"      "80380"    
    ## [106] "80381"     "8174"      "83700"     "8506"      "8516"     
    ## [111] "9019"      "9071"      "9073"      "9074"      "9075"     
    ## [116] "9076"      "9080"      "90952"     "914"       "920"      
    ## [121] "923"       "925"       "926"       "933"       "9369"     
    ## [126] "9378"      "9379"      "940"       "941"       "942"      
    ## [131] "947"       "958"       "959"       "965"       "9672"     
    ## [136] "999"      
    ## 
    ## $`hsa04010 MAPK signaling pathway`
    ##   [1] "10000"     "100137049" "10125"     "10235"     "10368"    
    ##   [6] "10369"     "10454"     "10746"     "10912"     "11072"    
    ##  [11] "11184"     "11221"     "11261"     "1147"      "115727"   
    ##  [16] "123745"    "1326"      "1386"      "1398"      "1399"     
    ##  [21] "1432"      "1616"      "1647"      "1649"      "1843"     
    ##  [26] "1844"      "1845"      "1846"      "1847"      "1848"     
    ##  [31] "1849"      "1850"      "1852"      "1950"      "1956"     
    ##  [36] "2002"      "2005"      "207"       "208"       "2122"     
    ##  [41] "2246"      "2247"      "2248"      "2249"      "2250"     
    ##  [46] "2251"      "2252"      "2253"      "2254"      "2255"     
    ##  [51] "2256"      "2257"      "2258"      "2259"      "2260"     
    ##  [56] "2261"      "2263"      "2264"      "22800"     "22808"    
    ##  [61] "23118"     "2316"      "23162"     "2317"      "2318"     
    ##  [66] "2353"      "23542"     "25780"     "26279"     "26281"    
    ##  [71] "26291"     "27006"     "27091"     "27092"     "27330"    
    ##  [76] "2768"      "2872"      "2885"      "30814"     "3164"     
    ##  [81] "3265"      "3303"      "3304"      "3305"      "3306"     
    ##  [86] "3310"      "3312"      "3315"      "355"       "3551"     
    ##  [91] "3552"      "3553"      "3554"      "356"       "3725"     
    ##  [96] "3727"      "3845"      "391013"    "3925"      "408"      
    ## [101] "409"       "4137"      "4149"      "4208"      "4214"     
    ## [106] "4215"      "4216"      "4217"      "4296"      "4342"     
    ## [111] "4609"      "4616"      "468"       "4763"      "4773"     
    ## [116] "4776"      "4790"      "4791"      "4803"      "4893"     
    ## [121] "4908"      "4909"      "4914"      "4915"      "50487"    
    ## [126] "5058"      "5062"      "51295"     "51347"     "5154"     
    ## [131] "5155"      "5156"      "5159"      "51701"     "51776"    
    ## [136] "5319"      "5320"      "5321"      "5322"      "5494"     
    ## [141] "5495"      "5530"      "5532"      "5533"      "5534"     
    ## [146] "5535"      "5536"      "5566"      "5567"      "5568"     
    ## [151] "5578"      "5579"      "55799"     "5582"      "5594"     
    ## [156] "5595"      "55970"     "5598"      "5599"      "5600"     
    ## [161] "5601"      "5602"      "5603"      "5604"      "5605"     
    ## [166] "5606"      "5607"      "5608"      "5609"      "5613"     
    ## [171] "56940"     "57551"     "5778"      "5801"      "5871"     
    ## [176] "5879"      "5880"      "5881"      "5894"      "5906"     
    ## [181] "5908"      "5921"      "5922"      "5923"      "5924"     
    ## [186] "59283"     "59284"     "59285"     "5970"      "5971"     
    ## [191] "6195"      "6196"      "6197"      "6237"      "627"      
    ## [196] "6300"      "63928"     "6416"      "64600"     "6654"     
    ## [201] "6655"      "6722"      "673"       "6788"      "6789"     
    ## [206] "6885"      "7040"      "7042"      "7043"      "7046"     
    ## [211] "7048"      "7124"      "7132"      "7157"      "7186"     
    ## [216] "7189"      "773"       "774"       "775"       "776"      
    ## [221] "777"       "778"       "7786"      "779"       "781"      
    ## [226] "782"       "783"       "784"       "785"       "7850"     
    ## [231] "786"       "7867"      "8074"      "80824"     "81579"    
    ## [236] "836"       "8398"      "8399"      "84647"     "84867"    
    ## [241] "8491"      "8517"      "8550"      "8569"      "8649"     
    ## [246] "8681"      "8817"      "8822"      "8823"      "8911"     
    ## [251] "8912"      "8913"      "8986"      "9020"      "9064"     
    ## [256] "9175"      "9252"      "9254"      "9261"      "929"      
    ## [261] "9344"      "93589"     "9448"      "9479"      "9693"     
    ## [266] "994"       "9965"      "998"      
    ## 
    ## $`hsa04012 ErbB signaling pathway`
    ##  [1] "10000"  "1026"   "1027"   "10298"  "10718"  "1398"   "1399"  
    ##  [8] "145957" "1839"   "1950"   "1956"   "1978"   "2002"   "2064"  
    ## [15] "2065"   "2066"   "2069"   "207"    "208"    "23533"  "23624" 
    ## [22] "2475"   "25"     "2549"   "25759"  "27"     "2885"   "2932"  
    ## [29] "3084"   "3265"   "369"    "3725"   "374"    "3845"   "399694"
    ## [36] "4609"   "4690"   "4893"   "5058"   "5062"   "5063"   "5290"  
    ## [43] "5291"   "5293"   "5294"   "5295"   "5296"   "5335"   "53358" 
    ## [50] "5336"   "5578"   "5579"   "5582"   "5594"   "5595"   "5599"  
    ## [57] "5601"   "5602"   "5604"   "5605"   "5609"   "56924"  "57144" 
    ## [64] "572"    "5747"   "5894"   "6198"   "6199"   "6416"   "6464"  
    ## [71] "6654"   "6655"   "6714"   "673"    "6776"   "6777"   "685"   
    ## [78] "7039"   "815"    "816"    "817"    "818"    "8440"   "8503"  
    ## [85] "867"    "868"    "9542"

``` r
#Shows first pathway involves caffeine metabolism. Next is drug metabolism. Next purine metabolis, etc... -> This is before overlap
#We need to see if pathways overlap to see any "enriched pathways"
```

\#The main gage() function requires a named vector of fold changes,
where the names of the values are the Entrez gene IDs.

\#Note that we used the mapIDs() function above to obtain Entrez gene
IDs (stored in
res*e**n**t**r**e**z*)*a**n**d**w**e**h**a**v**e**t**h**e**f**o**l**d**c**h**a**n**g**e**r**e**s**u**l**t**s**f**r**o**m**D**E**S**e**q*2*a**n**a**l**y**s**i**s*(*s**t**o**r**e**d**i**n**r**e**s*log2FoldChange).

Create a vector of FoldChange values that has ENTREZ identifiers as the
names of the vector. This is the format that the **gage()** function
wants.

``` r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
```

    ##      1266     54855      1465     51232      2034      2317 
    ## -2.422719  3.201955 -2.313738 -2.059631 -1.888019 -1.649792

Now, let’s run the gage pathway analysis.

``` r
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs)

attributes(keggres)
```

    ## $names
    ## [1] "greater" "less"    "stats"

See help on the gage function with ?gage. Specifically, you might want
to try changing the value of same.dir. This value determines whether to
test for changes in a gene set toward a single direction (all genes up
or down regulated) or changes towards both directions simultaneously
(i.e. any genes in the pathway dysregulated). Here, we’re using the
default same.dir=TRUE, which will give us separate lists for pathways
that are upregulated versus pathways that are down-regulated.

Now lets look at the object returned from gage().

``` r
attributes(keggres)
```

    ## $names
    ## [1] "greater" "less"    "stats"

\#It is a list with three elements, “greater”, “less” and “stats”.

\#You can also see this in your Environmnet panel/tab window of RStudio
or use the R command str(keggres).

\#Like any list we can use the dollar syntax to access a named element,
e.g. head(keggres*g**r**e**a**t**e**r*)*a**n**d**h**e**a**d*(*k**e**g**g**r**e**s*less).

\#Lets look at the first few down (less) pathway results:

``` r
# Look at the first few down (less) pathways
head(keggres$less)
```

    ##                                          p.geomean stat.mean        p.val
    ## hsa04110 Cell cycle                   8.995727e-06 -4.378644 8.995727e-06
    ## hsa03030 DNA replication              9.424076e-05 -3.951803 9.424076e-05
    ## hsa03013 RNA transport                1.375901e-03 -3.028500 1.375901e-03
    ## hsa03440 Homologous recombination     3.066756e-03 -2.852899 3.066756e-03
    ## hsa04114 Oocyte meiosis               3.784520e-03 -2.698128 3.784520e-03
    ## hsa00010 Glycolysis / Gluconeogenesis 8.961413e-03 -2.405398 8.961413e-03
    ##                                             q.val set.size         exp1
    ## hsa04110 Cell cycle                   0.001448312      121 8.995727e-06
    ## hsa03030 DNA replication              0.007586381       36 9.424076e-05
    ## hsa03013 RNA transport                0.073840037      144 1.375901e-03
    ## hsa03440 Homologous recombination     0.121861535       28 3.066756e-03
    ## hsa04114 Oocyte meiosis               0.121861535      102 3.784520e-03
    ## hsa00010 Glycolysis / Gluconeogenesis 0.212222694       53 8.961413e-03

\#Each keggres*l**e**s**s**a**n**d**k**e**g**g**r**e**s*greater object
is data matrix with gene sets as rows sorted by p-value.

\#The top “less/down” pathways is “Cell cycle” with the KEGG pathway
identifier hsa04110.

\#Now, let’s try out the pathview() function from the pathview package
to make a pathway plot with our RNA-Seq expression results shown in
color. To begin with lets manually supply a pathway.id (namely the first
part of the “hsa04110 Cell cycle”) that we could see from the print out
above.

``` r
pathview(gene.data=foldchanges, pathway.id="hsa04110")
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ## Info: Working in directory /Users/marygarcia/Desktop/bimm143_github/class15

    ## Info: Writing image file hsa04110.pathview.png

![](hsa04110.pathview.png)

Section 4. Reactome Analysis Reactome is database consisting of
biological molecules and their relation to pathways and processes.
Reactome, such as many other tools, has an online software available
(<a href="https://reactome.org/" class="uri">https://reactome.org/</a>)
and R package available
(<a href="https://bioconductor.org/packages/release/bioc/html/ReactomePA.html" class="uri">https://bioconductor.org/packages/release/bioc/html/ReactomePA.html</a>).

If you would like more information, the documentation is available here:
<a href="https://reactome.org/user/guide" class="uri">https://reactome.org/user/guide</a>

Let’s now conduct over-representation enrichment analysis and
pathway-topology analysis with Reactome using the previous list of
significant genes generated from our differential expression results
above.

First, Using R, output the list of significant genes at the 0.05 level
as a plain text file:

``` r
sig_genes <- res[res$padj <= 0.05 & !is.na(res$padj), "symbol"]
print(paste("Total number of significant genes:", length(sig_genes)))
```

    ## [1] "Total number of significant genes: 8147"

``` r
write.table(sig_genes, file="significant_genes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)
```

Then, to perform pathway analysis online go to the Reactome website
(<a href="https://reactome.org/PathwayBrowser/#TOOL=AT" class="uri">https://reactome.org/PathwayBrowser/#TOOL=AT</a>).
Select “choose file” to upload your significant gene list. Then, select
the parameters “Project to Humans”, then click “Analyze”.
