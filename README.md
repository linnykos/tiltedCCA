# Purpose

This repository contains all the functions to perform Tilted-CCA and the downstream analysis, for the paper "Quantifying common and distinct information in single-cell multimodal data with Tilted-CCA". See the companion GitHub package https://github.com/linnykos/tiltedCCA_analysis for all the analyses performed in the paper.

This code was developed and tested primarily on R 4.1.2. on a Macbook (macOS 11.6.8 Big Sur) equipped with an i7 processor.

# Installation

This package can be installed through `devtools` in R.

```{r}
library("devtools")
devtools::install_github("linnykos/tiltedCCA")
```
The package itself depends on several packages. These include `MASS`, `RSpectra`, `Matrix`, `Seurat`, `SeuratObject`, `irlba`, `RANN`, `dbscan`, `scales`, `sparseMatrixStats`, `RColorBrewer`, `ggplot2`, `ggrepel`, and `quadprog`. See the last section of this README to see where (i.e., CRAN, Bioconductor, or GitHub) to download all such packages.

# Setup

The following shows the suggested package versions that the developer (GitHub username: linnykos) used when developing the Tilted-CCA package.

```
> devtools::session_info()
─ Session info ─────────────────────────────────────────────────────
 setting  value
 version  R version 4.1.2 (2021-11-01)
 os       Red Hat Enterprise Linux
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/New_York
 date     2022-10-04
 pandoc   1.12.3.1 @ /usr/bin/pandoc

─ Packages ─────────────────────────────────────────────────────────
 package           * version    date (UTC) lib source
 abind               1.4-5      2016-07-21 [1] CRAN (R 4.1.2)
 assertthat          0.2.1      2019-03-21 [1] CRAN (R 4.1.2)
 brio                1.1.3      2021-11-30 [1] CRAN (R 4.1.2)
 cachem              1.0.6      2021-08-19 [1] CRAN (R 4.1.2)
 callr               3.7.1      2022-07-13 [1] CRAN (R 4.1.2)
 cli                 3.3.0      2022-04-25 [1] CRAN (R 4.1.2)
 cluster             2.1.2      2021-04-17 [2] CRAN (R 4.1.2)
 codetools           0.2-18     2020-11-04 [2] CRAN (R 4.1.2)
 colorspace          2.0-3      2022-02-21 [1] CRAN (R 4.1.2)
 cowplot             1.1.1      2020-12-30 [1] CRAN (R 4.1.2)
 crayon              1.5.1      2022-03-26 [1] CRAN (R 4.1.2)
 data.table          1.14.2     2021-09-27 [1] CRAN (R 4.1.2)
 DBI                 1.1.3      2022-06-18 [1] CRAN (R 4.1.2)
 dbscan            * 1.1-10     2022-01-15 [1] CRAN (R 4.1.2)
 deldir              1.0-6      2021-10-23 [1] CRAN (R 4.1.2)
 devtools            2.4.4      2022-07-20 [1] CRAN (R 4.1.2)
 digest              0.6.29     2021-12-01 [1] CRAN (R 4.1.2)
 dplyr               1.0.9      2022-04-28 [1] CRAN (R 4.1.2)
 ellipsis            0.3.2      2021-04-29 [1] CRAN (R 4.1.2)
 fansi               1.0.3      2022-03-24 [1] CRAN (R 4.1.2)
 fastmap             1.1.0      2021-01-25 [1] CRAN (R 4.1.2)
 fitdistrplus        1.1-8      2022-03-10 [1] CRAN (R 4.1.2)
 fs                  1.5.2      2021-12-08 [1] CRAN (R 4.1.2)
 future              1.27.0     2022-07-22 [1] CRAN (R 4.1.2)
 future.apply        1.9.0      2022-04-25 [1] CRAN (R 4.1.2)
 generics            0.1.3      2022-07-05 [1] CRAN (R 4.1.2)
 ggplot2           * 3.3.6      2022-05-03 [1] CRAN (R 4.1.2)
 ggrepel           * 0.9.1      2021-01-15 [1] CRAN (R 4.1.2)
 ggridges            0.5.3      2021-01-08 [1] CRAN (R 4.1.2)
 globals             0.16.0     2022-08-05 [1] CRAN (R 4.1.2)
 glue                1.6.2      2022-02-24 [1] CRAN (R 4.1.2)
 goftest             1.2-3      2021-10-07 [1] CRAN (R 4.1.2)
 gridExtra           2.3        2017-09-09 [1] CRAN (R 4.1.2)
 gtable              0.3.0      2019-03-25 [1] CRAN (R 4.1.2)
 htmltools           0.5.3      2022-07-18 [1] CRAN (R 4.1.2)
 htmlwidgets         1.5.4      2021-09-08 [1] CRAN (R 4.1.2)
 httpuv              1.6.5      2022-01-05 [1] CRAN (R 4.1.2)
 httr                1.4.3      2022-05-04 [1] CRAN (R 4.1.2)
 ica                 1.0-3      2022-07-08 [1] CRAN (R 4.1.2)
 igraph              1.3.4      2022-07-19 [1] CRAN (R 4.1.2)
 irlba             * 2.3.5      2021-12-06 [1] CRAN (R 4.1.2)
 jsonlite            1.8.0      2022-02-22 [1] CRAN (R 4.1.2)
 KernSmooth          2.23-20    2021-05-03 [2] CRAN (R 4.1.2)
 later               1.3.0      2021-08-18 [1] CRAN (R 4.1.2)
 lattice             0.20-45    2021-09-22 [2] CRAN (R 4.1.2)
 lazyeval            0.2.2      2019-03-15 [1] CRAN (R 4.1.2)
 leiden              0.4.2      2022-05-09 [1] CRAN (R 4.1.2)
 lifecycle           1.0.1      2021-09-24 [1] CRAN (R 4.1.2)
 listenv             0.8.0      2019-12-05 [1] CRAN (R 4.1.2)
 lmtest              0.9-40     2022-03-21 [1] CRAN (R 4.1.2)
 magrittr            2.0.3      2022-03-30 [1] CRAN (R 4.1.2)
 MASS              * 7.3-54     2021-05-03 [2] CRAN (R 4.1.2)
 Matrix            * 1.3-4      2021-06-01 [2] CRAN (R 4.1.2)
 MatrixGenerics    * 1.6.0      2021-10-26 [1] Bioconductor
 matrixStats       * 0.62.0     2022-04-19 [1] CRAN (R 4.1.2)
 memoise             2.0.1      2021-11-26 [1] CRAN (R 4.1.2)
 mgcv                1.8-38     2021-10-06 [2] CRAN (R 4.1.2)
 mime                0.12       2021-09-28 [1] CRAN (R 4.1.2)
 miniUI              0.1.1.1    2018-05-18 [1] CRAN (R 4.1.2)
 munsell             0.5.0      2018-06-12 [1] CRAN (R 4.1.2)
 nlme                3.1-153    2021-09-07 [2] CRAN (R 4.1.2)
 parallelly          1.32.1     2022-07-21 [1] CRAN (R 4.1.2)
 patchwork           1.1.1      2020-12-17 [1] CRAN (R 4.1.2)
 pbapply             1.5-0      2021-09-16 [1] CRAN (R 4.1.2)
 pillar              1.8.0      2022-07-18 [1] CRAN (R 4.1.2)
 pkgbuild            1.3.1      2021-12-20 [1] CRAN (R 4.1.2)
 pkgconfig           2.0.3      2019-09-22 [1] CRAN (R 4.1.2)
 pkgload             1.3.0      2022-06-27 [1] CRAN (R 4.1.2)
 plotly              4.10.0     2021-10-09 [1] CRAN (R 4.1.2)
 plyr                1.8.7      2022-03-24 [1] CRAN (R 4.1.2)
 png                 0.1-7      2013-12-03 [1] CRAN (R 4.1.2)
 polyclip            1.10-0     2019-03-14 [1] CRAN (R 4.1.2)
 prettyunits         1.1.1      2020-01-24 [1] CRAN (R 4.1.2)
 processx            3.7.0      2022-07-07 [1] CRAN (R 4.1.2)
 profvis             0.3.7      2020-11-02 [1] CRAN (R 4.1.2)
 progressr           0.10.1     2022-06-03 [1] CRAN (R 4.1.2)
 promises            1.2.0.1    2021-02-11 [1] CRAN (R 4.1.2)
 ps                  1.7.0      2022-04-23 [1] CRAN (R 4.1.2)
 purrr               0.3.4      2020-04-17 [1] CRAN (R 4.1.2)
 quadprog          * 1.5-8      2019-11-20 [1] CRAN (R 4.1.2)
 R6                  2.5.1      2021-08-19 [1] CRAN (R 4.1.2)
 RANN              * 2.6.1      2019-01-08 [1] CRAN (R 4.1.2)
 RColorBrewer      * 1.1-3      2022-04-03 [1] CRAN (R 4.1.2)
 Rcpp                1.0.9      2022-07-08 [1] CRAN (R 4.1.2)
 RcppAnnoy           0.0.19     2021-07-30 [1] CRAN (R 4.1.2)
 remotes             2.4.2      2021-11-30 [1] CRAN (R 4.1.2)
 reshape2            1.4.4      2020-04-09 [1] CRAN (R 4.1.2)
 reticulate          1.25       2022-05-11 [1] CRAN (R 4.1.2)
 rgeos               0.5-9      2021-12-15 [1] CRAN (R 4.1.2)
 rlang               1.0.4      2022-07-12 [1] CRAN (R 4.1.2)
 ROCR                1.0-11     2020-05-02 [1] CRAN (R 4.1.2)
 rpart               4.1-15     2019-04-12 [2] CRAN (R 4.1.2)
 RSpectra          * 0.16-1     2022-04-24 [1] CRAN (R 4.1.2)
 Rtsne               0.16       2022-04-17 [1] CRAN (R 4.1.2)
 scales            * 1.2.0      2022-04-13 [1] CRAN (R 4.1.2)
 scattermore         0.8        2022-02-14 [1] CRAN (R 4.1.2)
 sctransform         0.3.3      2022-01-13 [1] CRAN (R 4.1.2)
 sessioninfo         1.2.2      2021-12-06 [1] CRAN (R 4.1.2)
 Seurat            * 4.1.1      2022-05-02 [1] CRAN (R 4.1.2)
 SeuratObject      * 4.1.0      2022-05-01 [1] CRAN (R 4.1.2)
 shiny               1.7.2      2022-07-19 [1] CRAN (R 4.1.2)
 sp                * 1.5-0      2022-06-05 [1] CRAN (R 4.1.2)
 sparseMatrixStats * 1.6.0      2021-10-26 [1] Bioconductor
 spatstat.core       2.4-2      2022-04-01 [1] CRAN (R 4.1.2)
 spatstat.data       2.2-0      2022-04-18 [1] CRAN (R 4.1.2)
 spatstat.geom       2.4-0      2022-03-29 [1] CRAN (R 4.1.2)
 spatstat.random     2.2-0      2022-03-30 [1] CRAN (R 4.1.2)
 spatstat.sparse     2.1-1      2022-04-18 [1] CRAN (R 4.1.2)
 spatstat.utils      2.3-1      2022-05-06 [1] CRAN (R 4.1.2)
 stringi             1.7.8      2022-07-11 [1] CRAN (R 4.1.2)
 stringr             1.4.0      2019-02-10 [1] CRAN (R 4.1.2)
 survival            3.2-13     2021-08-24 [2] CRAN (R 4.1.2)
 tensor              1.5        2012-05-05 [1] CRAN (R 4.1.2)
 testthat          * 3.1.4      2022-04-26 [1] CRAN (R 4.1.2)
 tibble              3.1.8      2022-07-22 [1] CRAN (R 4.1.2)
 tidyr               1.2.0      2022-02-01 [1] CRAN (R 4.1.2)
 tidyselect          1.1.2      2022-02-21 [1] CRAN (R 4.1.2)
 tiltedCCA         * 0.0.0.0114 2022-09-06 [1] local
 urlchecker          1.0.1      2021-11-30 [1] CRAN (R 4.1.2)
 usethis             2.1.6      2022-05-25 [1] CRAN (R 4.1.2)
 utf8                1.2.2      2021-07-24 [1] CRAN (R 4.1.2)
 uwot                0.1.13     2022-08-16 [1] CRAN (R 4.1.2)
 vctrs               0.4.1      2022-04-13 [1] CRAN (R 4.1.2)
 viridisLite         0.4.0      2021-04-13 [1] CRAN (R 4.1.2)
 withr               2.5.0      2022-03-03 [1] CRAN (R 4.1.2)
 xtable              1.8-4      2019-04-21 [1] CRAN (R 4.1.2)
 zoo                 1.8-10     2022-04-15 [1] CRAN (R 4.1.2)
```