**rewind**: **R**econstructing **E**tiology **w**ith B**in**ary **D**ecomposition

An R package for fitting Bayesian restricted latent class models. 

zhenkewu badges:
[![Travis CI Build Status](https://travis-ci.org/zhenkewu/rewind.svg?branch=master)](https://travis-ci.org/zhenkewu/rewind)

**Maintainer**: Zhenke Wu, zhenkewu@umich.edu

**References**: If you are using **rewind** for clustering multivariate binary
observations with restricted latent class models that accounts for differential
measurement errors, please cite the following paper:

|       | Citation     | Paper Link
| -------------  | -------------  | -------------  |
| Bayesian restricted LCM    | Wu, Z. and Zeger, S. L. (2018+), Clustering Multivariate Binary Outcomes with Restricted Latent Class Models: A Bayesian Approach. Working Paper   |[Link]()| 


## Table of content
- [1. Installation](#id-section1)
- [2. Overview](#id-section2)
- [2. Example](#id-section3)

<div id='id-section1'/>

Installation
--------------
```r
install.packages("devtools",repos="https://cloud.r-project.org")
devtools::install_github("zhenkewu/rewind")
```
<div id='id-section2'/>

Overview
----------
It is common to observe multiple binary data where the goal is to filter observations into homogeneous groups, or clustering. The measurements typically exhibit differential errors in measuring traits such as the true presence or absence of a list of pathogens in disease epidemiology or of antibodies in medicine, or the true capability or incapability of answering test questions in educational assessment. This paper studies restricted latent class models designed to account for such differential errors. In three parts, the models are specified by 1) unobserved low-dimensional multivariate binary states for each class, 2) a binary design matrix specifying the most accurately measured class(es) for each trait, and 3) a product-Bernoulli likelihood. Identifiablilty conditions based on likelihood are discussed which inform the posterior algorithm design. In finite samples and dimensions, we encourage *in a priori* low model complexity via few classes with distinct sparse latent states. The posterior distribution of the number of clusters and the patterns of latent states are estimated from data and by design tends to concentrate on smaller values and sparser patterns, respectively. We use Markov chain Monte Carlo to approximate the posterior distribution of unknown parameters to estimate clusters, individual latent states and the design matrix with uncertainty. The model readily extends to studies where some subjects' true latent classes are known or important prior knowledge about measurement accuracy is available from external sources. The methods are illustrated with analyses of two data sets from medicine and educational assessment.


<div id='id-section3'/>

Examples
---------

- 1. binary factor analysis with **pre-specified # of factors** but unknown # of clusters
```r
library(rewind)
example(sampler)
```

