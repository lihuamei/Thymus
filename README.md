# Thymus
Unraveling the spatial organization and development of human thymocytes through integration of spatial transcriptomics and single-cell multi-omics profiling

<p align="center">
	<img src="vignette_files/Thymus.jpg" alt="Resized Image" width="800">
</p>

<b> All the analysis codes used in our manuscript are provided in the `source.code` directory. Raw and preprocessed data can be obtained from the provided URL: https://ngdc.cncb.ac.cn/bioproject/ </b>

## TSO-His

TSO-his is a specialized tool designed for the identification of the cortex, medulla, and corticomedullary junction regions in thymic ST sections. This tool leverages statistical testing and network search algorithms to achieve its functionality.

<p align="center">
	<img src="vignette_files/TSO.His.jpg" alt="Resized Image" width="800">
</p>

## 1. How to install

``` r
library(devtools)
install_github("lihuamei/Thymus/thymusTSO")

``` 

## 2. Loading `thymusTSO` package and testing
``` r
library(thymusTSO)
sp.obj <- system.file('data/thymus_T2.RDS', package = 'thymusTSO') %>% readRDS
sp.obj <- tsoHis(sp.obj)
SpatialPlot(sp.obj[[1]], group.by = 'HE.Labels', cols = c('grey', 'red', 'green', 'pink', 'yellow') %>% `names<-`(unique(sp.obj$HE.Labels))

fitDistLinesByGlm(sp.obj, plot.tar = 'CCL25', degree = 10)
fitDistLinesByWindows(sp.obj, plot.tar = c('CCL25', 'CCL19', 'CD19', 'RAG1'), win = 20)

``` 
<p align="center">
	<img src="vignette_files/exam.1.jpg" alt="Resized Image" width="800">
</p>

If calling XGBoost model for predicting cortical and medullary spots, please set `call.xgb = TRUE`.

``` r
sp.obj <- tsoHis(sp.obj, call.xgb = TRUE)

``` 