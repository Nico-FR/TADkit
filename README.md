
<!-- README.md is generated from README.Rmd. Please edit that file -->

# 1 TADkit

The main objective of the TADkit R package is to make the visualization
of data associated with HiC analyses accessible to biologists without
the need for bioinformatics skills. To this end, two pairs of functions
(based on gviz and ggplot packages) have been created to visualize :

-   domains such as TADs (topological Associated Domain) or compartments
    (compartment A and B),
-   interaction matrices.

Other tools has been added in order to analysed genomic annotations in
the light of the 3D organisation.

In this tutorial, we are going to start with the visualization of
domains and matrices with the datas included in the package. Note that
most of the functions refer to the so called TAD but it also work for
other kind of domains (e.g compartments).

## 1.1 Installation

You can install the development version of TADkit from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Nico-FR/TADkit")
```

# 2 Data format

First, we will see the different types of data and their formats that
can be visualized with TADkit. Datas of 2 chromosomes and 2 cows are
available within the package. Let’s start by loading the packages needed
for that tutorial.

``` r
library(TADkit)
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
```

## 2.1 matrix

Interaction count (HiC matrix) for a genome is a very large size and can
be store in a multiple formats. HiC matrix is symmetrical with respect
to the diagonal, so only the upper part of the diagonal is loaded in
memory for only one chromosome (intra chromosomal interactions).

### 2.1.1 data frame

A basic storage is a compressed data frame for each chromosome in which
each row and columns represent a bin. Bovine matrices for chromosome 25
and 26 are available in the TADkit package as a sparse `dgCMatrix`
format (i.e upper part of the matrix without the zeros).

``` r
class(TADkit::matrix_1_chr25_50kb)
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"

#raw count between bins 60 to 70
TADkit::matrix_1_chr25_50kb[60:70,60:70]
#> 11 x 11 sparse Matrix of class "dgCMatrix"
#>                                                             
#>  [1,] 3959 1024  312  256  325  269  181  133  165  269   83
#>  [2,]    . 2528 1046  480  596  435  302  273  258  541  208
#>  [3,]    .    . 4581 1967  982  526  325  345  316  529  186
#>  [4,]    .    .    . 4136 1832  577  359  291  302  467  165
#>  [5,]    .    .    .    . 3721 1597  631  423  463  643  166
#>  [6,]    .    .    .    .    . 4037 1557  723  650  982  156
#>  [7,]    .    .    .    .    .    . 4171 1825 1091 1241  211
#>  [8,]    .    .    .    .    .    .    . 4616 1687 1019  251
#>  [9,]    .    .    .    .    .    .    .    . 4397 1941  255
#> [10,]    .    .    .    .    .    .    .    .    . 4281  661
#> [11,]    .    .    .    .    .    .    .    .    .    . 4690
```

Let’s write the matrix as data frame on the current directory:

``` r
write.table(as.data.frame(as.matrix(TADkit::matrix_1_chr25_50kb)), #dgCMatrix 2 data frame
            "./matrix_1_chr25_50kb.df", 
            row.names = TRUE, #add row names
            col.names = TRUE, #add column names
            quote = FALSE, sep = "\t") 
```

Those data frames can be load in R as a matrix:

``` r
mat1.df = read.table("./matrix_1_chr25_50kb.df",
                 skip = 1, #skip col.names
                 sep = "\t",
                 row.names = 1 #specify the columns number with row names
                 ) 
mat1.mat = as.matrix(mat1.df) # translate to matrix
mat1.dgcmat = as(mat1.mat, "dgCMatrix") #translate to dgCMatrix
```

### 2.1.2 .cool / .mcool

One of the most common HiC array storage formats is the .cool or .mcool
format which store matrices for 1 resolution and multiple resolutions
respectively. To load HiC matrix for one chromosome from cool/.mcool as
a matrix (`dgCMatrix`) use `coolFetch()`:

``` r
#read .cool file
my_matrix.10kb = coolFetch("my_matrix.10kb.cool", chr = "chr_name")

#read .mcool file
my_matrix.10kb = coolFetch("my_matrix.mcool", chr = "chr_name", bin.width = 10e3)
```

## 2.2 .bed

Domains are usely stored in bed file format which contains at least 3
columns (chromosome name, the start and the end of the TAD).

``` r
head(TADkit::tad_1_10kb.bed) #TAD for bovine 1 estimated from 10kb matrix
#>   chr   start     end
#> 1  25  360000  550000
#> 2  25  550000  680000
#> 3  25  680000  870000
#> 4  25  870000 1050000
#> 5  25 1050000 1360000
#> 6  25 1360000 1520000
```

Other columns can be added like the stand, any score or character. In
most of the TADkit functions bed files must be stored in `GRanges` in
which we can store the size of the chromosomes we will need.  
To create a `GRanges` with the size of the chromosomes use
`dataframe2grange()`:

``` r
#chromosomes sizes
TADkit::chromsize
#>   chr     size
#> 1  25 42350435
#> 2  26 51992305

#create GRanges with chr sizes:
tad_1_10kb.gr = dataframes2grange(TADkit::tad_1_10kb.bed, TADkit::chromsize)
tad_1_10kb.gr
#> GRanges object with 262 ranges and 0 metadata columns:
#>               seqnames            ranges strand
#>                  <Rle>         <IRanges>  <Rle>
#>     25_360000       25     360000-550000      *
#>     25_550000       25     550000-680000      *
#>     25_680000       25     680000-870000      *
#>     25_870000       25    870000-1050000      *
#>    25_1050000       25   1050000-1360000      *
#>           ...      ...               ...    ...
#>   26_50505000       26 50505000-50615000      *
#>   26_50615000       26 50615000-50890000      *
#>   26_50890000       26 50890000-51140000      *
#>   26_51140000       26 51140000-51250000      *
#>   26_51250000       26 51250000-51670000      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome

#chromsize within GRanges object
seqlengths(tad_1_10kb.gr)
#>       25       26 
#> 42350435 51992305
```

TADs are often stored in the form of domains and each domain gives the
beginning and the end of a TAD. Sometime TADs are stored in the form of
boundaries which gives bins in which borders are located. The width of
the boundaries corresponds to the size of the bins of the matrix (i.e
10kb in our datas).

``` r
#TADs as a form of boundaries:
boundaries = data.frame(chr = seqnames(tad_1_10kb.gr),
                      start = start(tad_1_10kb.gr) - 5e3,
                      end = start(tad_1_10kb.gr) + 5e3)
head(boundaries)
#>   chr   start     end
#> 1  25  355000  365000
#> 2  25  545000  555000
#> 3  25  675000  685000
#> 4  25  865000  875000
#> 5  25 1045000 1055000
#> 6  25 1355000 1365000
```

Then boundaries can be translated as domains with `boundary2domain()`:

``` r
boundary2domain(boundaries)
#> GRanges object with 260 ranges and 0 metadata columns:
#>               seqnames            ranges strand
#>                  <Rle>         <IRanges>  <Rle>
#>     25_360000       25     360000-550000      *
#>     25_550000       25     550000-680000      *
#>     25_680000       25     680000-870000      *
#>     25_870000       25    870000-1050000      *
#>    25_1050000       25   1050000-1360000      *
#>           ...      ...               ...    ...
#>   26_50125000       26 50125000-50505000      *
#>   26_50505000       26 50505000-50615000      *
#>   26_50615000       26 50615000-50890000      *
#>   26_50890000       26 50890000-51140000      *
#>   26_51140000       26 51140000-51250000      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```

In addition to TADs, we can also plot any type of annotations. For this
tutorial we will use genes annotations available on Ensembl:

``` r
txdb <- makeTxDbFromBiomart(biomart="ensembl", dataset="btaurus_gene_ensembl")
#> Download and preprocess the 'transcripts' data frame ... OK
#> Download and preprocess the 'chrominfo' data frame ... OK
#> Download and preprocess the 'splicings' data frame ... OK
#> Download and preprocess the 'genes' data frame ... OK
#> Prepare the 'metadata' data frame ... OK
#> Make the TxDb object ... OK
genomic.gr = genes(txdb, columns=c("TXTYPE"))
genomic.gr
#> GRanges object with 27607 ranges and 1 metadata column:
#>                      seqnames              ranges strand |          TXTYPE
#>                         <Rle>           <IRanges>  <Rle> | <CharacterList>
#>   ENSBTAG00000000005       17   65389743-65505336      + |  protein_coding
#>   ENSBTAG00000000008       29   32214439-32244810      - |  protein_coding
#>   ENSBTAG00000000009       18   12338037-12342272      + |  protein_coding
#>   ENSBTAG00000000010       21   34209956-34223394      + |  protein_coding
#>   ENSBTAG00000000011        8     7950815-7971600      - |  protein_coding
#>                  ...      ...                 ...    ... .             ...
#>   ENSBTAG00000055312        2   96551552-96557130      + |  protein_coding
#>   ENSBTAG00000055313       21   65169462-65169520      + |           miRNA
#>   ENSBTAG00000055314       19   24021135-24051219      - |  protein_coding
#>   ENSBTAG00000055315        8 102574538-102575815      + |  protein_coding
#>   ENSBTAG00000055316        2   94738416-94738887      + |  protein_coding
#>   -------
#>   seqinfo: 2211 sequences (1 circular) from an unspecified genome
```

## 2.3 .bedgraph

Bedgraph are used to store a score for each bin, like insulation score.
Insulation score (IS) is defined for a bin as an average number of
intyeractions that occur across this bin in some vicinity of the bin
(Crane et al., 2015). Insulation scores are used to call TAD boundaries.
In TADkit plot functions, bedgraph inputs can be in 3 formats:

-   `dataframe`,
-   `GRanges`,
-   or the path of the file.

To avoid having to load too much data in the R environment, it is
possible to specify the file paths instead of an R object:

``` r
# insulation score for indiv 1:
head(TADkit::IS_1_10kb.bedgraph)
#>   chr start   end        IS
#> 1  25  5000 15000 -1.646012
#> 2  25 15000 25000 -1.621941
#> 3  25 25000 35000 -1.552033
#> 4  25 35000 45000 -1.444361
#> 5  25 45000 55000 -1.265504
#> 6  25 55000 65000 -1.203068

#write bedgraph in a file
write.table(TADkit::IS_1_10kb.bedgraph, "./IS_1_10kb.bedgraph", col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t")
```

## 2.4 .bigwig

Similar to bedgraph, bigwig files are used to store a score for each
interval in an indexed and compressed form. Insulation scores stored as
a bigwig format can be translated to a `GRranges` with `bw2grange()`.

For this tutorial we are going to use coverage data from RNA sequencing
experiment available on Ensembl:

``` r
#download expression data (RNAseq) of bovine:
download.file("http://ftp.ensembl.org/pub/release-104/bamcov/bos_taurus/genebuild/ARS-UCD1.2.ENA.heart.1.bam.bw", destfile = "./rna_seq_bov.bw", method = "curl") #if it fails: remove the method parameter or download it manually in the current directory
```

# 3 Domains plot

Two functions have been built to plot TADs: for one indiviual
`TADplot()` and multiple individual `mTADplot()`.

In addition to domains, few others tracks can be plotted:

-   bedgraph to plot bin scores (e.g insulation score…),
-   bigwig to plot coverage datas (e.g RNAseq…),
-   bed to plot annotations (e.g genes…).

## 3.1 TADplot

Let start with the most basic usage:

``` r
TADplot(tad.gr = tad_1_10kb.gr, chr = 25, start = 13e6, stop = 15e6)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

Note that the area is extended to the first and last TAD of the window.

Add insulation score (bedgraph), RNAseq (bigwig) and genomic annotations
(genes grouped according to the “TXTYPE” column, i.e first metadata
column) to the plot:

``` r
TADplot(tad.gr = tad_1_10kb.gr, chr = 25, start = 13e6, stop = 15e6,
        bedgraph = IS_1_10kb.bedgraph,
        bigwig.path = "./rna_seq_bov.bw", 
        annot.gr = genomic.gr, 
        annot.col = 1, #columns to group annotations
        bigwig.yaxis = "log2" #log2 of RNAseq values
        )
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="100%" />

## 3.2 mTADplot

`mTADplot()` allow to plot 1 or more individual/sample in a same plot.
All tracks inputs must be stored in a list, and all object in the list
must have names:

To avoid having to load too much data in the R environment, we will see
how to use the paths of bedgraph files. Datas for the 2 bovines:

``` r
#create GRanges with TADs
tad_2_10kb.gr = dataframes2grange(TADkit::tad_2_10kb.bed, TADkit::chromsize) 

#write bedgraph as file
write.table(TADkit::IS_2_10kb.bedgraph, "./IS_2_10kb.bedgraph", col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t") 
write.table(TADkit::PC1_1_50kb.bedgraph, "./PC1_1_50kb.bedgraph", col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t") 
write.table(TADkit::PC1_2_50kb.bedgraph, "./PC1_2_50kb.bedgraph", col.names = FALSE, quote = FALSE, row.names = FALSE, sep = "\t") 
```

### 3.2.1 Create list

Let’s create lists for the 2 individuals (then named as “bov1” and
“bov2”):

-   TADs:

``` r
#list of TADs
tad.lst = list(tad_1_10kb.gr, tad_2_10kb.gr)
names(tad.lst) = c("bov1", "bov2")

#insulation score (path) list
IS.lst = list("./IS_1_10kb.bedgraph", "./IS_2_10kb.bedgraph")
names(IS.lst) = c("bov1", "bov2")
```

-   insulation scores:

``` r
#insulation score (path) list
IS.lst = list("./IS_1_10kb.bedgraph", "./IS_2_10kb.bedgraph")
names(IS.lst) = c("bov1", "bov2")
```

Do the plot:

``` r
mTADplot2(tad.lst = tad.lst, chr = 25, start = 13e6, stop = 15e6,
        bedgraph.lst = IS.lst, bedgraph.name = "IS")
```

<img src="man/figures/README-unnamed-chunk-19-1.png" width="100%" />

### 3.2.2 Create multiple lists

In case of several tracks of bedgraph, we can create a list containing
few other lists. Let’s create 2 lists with the insulation scores and the
PC1 scores:

``` r
#insulation scores (path) as list
IS.lst = list("./IS_1_10kb.bedgraph", "./IS_2_10kb.bedgraph")
names(IS.lst) = c("bov1", "bov2")

#PC1 scores (path) as list
PC1.lst = list("./PC1_1_50kb.bedgraph", "./PC1_2_50kb.bedgraph")
names(PC1.lst) = c("bov1", "bov2")
```

Now let’s create a list with those 2 lists:

``` r
#multiple lists fo bedgraph
bedgraphs.lst = list(IS.lst, PC1.lst)
names(bedgraphs.lst) = c("IS", "PC1")

bedgraphs.lst
#> $IS
#> $IS$bov1
#> [1] "./IS_1_10kb.bedgraph"
#> 
#> $IS$bov2
#> [1] "./IS_2_10kb.bedgraph"
#> 
#> 
#> $PC1
#> $PC1$bov1
#> [1] "./PC1_1_50kb.bedgraph"
#> 
#> $PC1$bov2
#> [1] "./PC1_2_50kb.bedgraph"
```

Do the plot:

``` r
mTADplot2(tad.lst = tad.lst, chr = 25, start = 13e6, stop = 15e6,
        bedgraph.lst = bedgraphs.lst)
```

<img src="man/figures/README-unnamed-chunk-22-1.png" width="100%" />

## 3.3 Options

When you master TADplot(), mTADplot() is much more powerful and flexible
to make the desired graphics. Moreover there are many options to change
some of the graphic parameters:

### 3.3.1 bigwigPath.lst

-   bigwig.binsize: Bin sizes for the histogram. Default = 1e3.

-   bigwig.xaxis: Function used to transform the x-axis among each
    bigwig.binsize. Defaults = “median”. Alternatively, other predefined
    functions can be supplied as character (“mean”, “median”, “sum”,
    “min”, “max” or “extreme”).

-   bigwig.chr: Chromosome name used to filter chromosome names that can
    be different from chr (e.g chr = “1” and bigwig.chr = “chr1”).
    Default = NULL to used the same name as chr.

-   bigwig.yaxis: Function used to transforming the y-axis values.
    Default = NULL. Use “log2” to use the function log2(x + 1) to
    transform the y-axis (for RNA seq) or provide any other function.

### 3.3.2 annotation.lst

-   annot.col: Column number of the metadata from annot.gr file(s) used
    to group the annotation tracks. Default = NULL and the name of each
    annotation is added.

### 3.3.3 bedgraphPath.lst

-   bedgraph.name: Name of the bedgraph track when there is only one
    track (default = “bedgraph”). Otherwise it takes the names of each
    list.

-   bedgraph_outliers: Ratio to remove outliers of all each bedgraph
    values. Default is 0 (ie no filter). To remove the first and last 2
    percentiles (i.e extreme values / outliers) use 0.02.

# 4 Matrix plot

Two functions have been built to plot interaction count (HiC matrix) for
one matrix `MATplot()` and 2 matrices `mMATplot()` (on the upper or
lower part) as a ggplot graph. In addition to matrices, 3 others tracks
can be plotted:

-   bed file type to visualized domains as lines or a triangles,
-   bedpe file type to hightlight interactions between 2 areas (loops).

## 4.1 MATplot

### 4.1.1 matrix + triangles

Plot of the matrix at 50kb resolution with TADs as triangles on the
upper part of the matrix.

``` r
MATplot(matrix = matrix_1_chr25_50kb, start = 10e6, stop = 30e6,
        bin.width = 50e3, log2 = T,
        tad.upper.tri = tad_1_10kb.gr,
        tad.chr = 25, #filter TADs for chr 25
        scale.colors = "H", #color of matrix, try "D" or "H"
        annotations.color = "red")+
  ggtitle("log2(raw count): chr25 bov1")
#> Warning: Removed 1 rows containing missing values (geom_segment).
#> Removed 1 rows containing missing values (geom_segment).
```

<img src="man/figures/README-unnamed-chunk-23-1.png" width="100%" />

### 4.1.2 matrix + triangles + loops + lines

bedpe file are used to highlight 2 areas. For example let’s create a
bedpe file to highlight a region on chromosome 25 between :

-   11Mb to 13Mb,
-   26Mb to 27Mb:

``` r
bedpe = data.frame(chr1 = "25", start1 = 11e6, end1 = 13e6,
                   chr2 = "25", start2 = 26e6, end2 = 27e6)
bedpe
#>   chr1  start1    end1 chr2  start2    end2
#> 1   25 1.1e+07 1.3e+07   25 2.6e+07 2.7e+07
```

In addition to bedpe file, let’s see 2 others functions:

-   `matObsExp()` to produce the observed / expected ratio of
    interaction count,
-   `PC1calling()` to call the compartments A and B from PC1 values.

``` r
MATplot(matrix = matObsExp(matrix_1_chr25_50kb), 
        start = 10e6, stop = 30e6,
        bin.width = 50e3, log2 = T,
        tad.upper.tri = tad_1_10kb.gr,
        tad.chr = 25, #filter TADs for chr 25
        scale.colors = "OE", 
        annotations.color = "red",
        loop.bedpe = bedpe,
        tad.upper.line = PC1calling(PC1_1_50kb.bedgraph), #compartments
        tad.line.col = 1 #use the fist metadata columns with vectors A and B
        )+
  ggtitle("log2(obs/exp): chr25 bov1")
#> Warning: Removed 1 rows containing missing values (geom_segment).
#> Removed 1 rows containing missing values (geom_segment).
```

<img src="man/figures/README-unnamed-chunk-25-1.png" width="100%" />

## 4.2 mMATplot

### 4.2.1 2 matrices + triangles

Like `MATplot()`, `mMATplot()` allow to plot 2 matrices: one on the
upper part of the plot and the other one on the lower part. Let’s plot
the 2 matrices and TADs for bov1 and bov2:

``` r
mMATplot(matrix.upper = matrix_1_chr25_50kb,
         matrix.lower = matrix_2_chr25_50kb,
         matrix.upper.txt = "bov1",
         matrix.lower.txt = "bov2",
         start = 10e6, stop = 30e6,
         bin.width = 50e3, log2 = T,
         tad.upper.tri = tad_1_10kb.gr,
         tad.lower.tri = tad_2_10kb.gr,
         tad.chr = 25)+
  ggtitle("log2(raw count): chr25 bov1 vs bov2")
#> Warning: Removed 1 rows containing missing values (geom_segment).
#> Removed 1 rows containing missing values (geom_segment).
#> Removed 1 rows containing missing values (geom_segment).
#> Removed 1 rows containing missing values (geom_segment).
```

<img src="man/figures/README-unnamed-chunk-26-1.png" width="100%" />

# 5 Clear files

``` r
file.remove(list.files(full.names = TRUE, pattern = ".bw"))
#> [1] TRUE
file.remove(list.files(full.names = TRUE, pattern = ".bedgraph"))
#> [1] TRUE TRUE TRUE TRUE
file.remove(list.files(full.names = TRUE, pattern = ".df"))
#> [1] TRUE
```
