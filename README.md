**rewind**: **R**econstructing **E**tiology **w**ith B**in**ary **D**ecomposition

An R package for fitting binary factor analysis models. 

**Maintainer**: Zhenke Wu, zhenkewu@umich.edu

## Table of content
- [1. Installation](#id-section1)
- [2. Example](#id-section2)

<div id='id-section1'/>

Installation
--------------
```r
# install.packages("devtools",repos="https://cloud.r-project.org")
devtools::install_github("zhenkewu/rewind")
```


<div id='id-section2'/>

Examples
---------

- 1. binary factor analysis with **pre-specified # of factors** but unknown # of clusters
```r
library(rewind)
example(sampler)
```

- 2. binary factor analysis with **unknown # of factors** but unknown # of clusters
```r
library(rewind)
example(slice_sampler)
```
