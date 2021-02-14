### PipeOne integrates multi-modal data from RNA-seq to identify disease features and clinically-relevant subtypes

![](https://img.shields.io/badge/license-MIT-brightgreen)
[![Documentation Status](https://readthedocs.org/projects/pipeone/badge/?version=latest)](https://pipeone.readthedocs.io/en/latest/?badge=latest)
[![](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A520.07.1.5413-brightgreen)](https://www.nextflow.io/)


[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
#### Overview

RNA sequencing (RNA-seq) represents one of the most widely-used technologies to investigate the transcriptome, and various algorithms have been developed to analyze RNA-seq data. However, a workflow integrating multi-modal information from these algorithms to study sequenced samples comprehensively is still lacking. Here, we present PipeOne, a cross-platform one-stop analysis workflow for a large number of RNA-seq samples. It includes three modules, data preprocessing and feature matrices construction by combining eight RNA analysis tools, disease feature prioritization by machine learning, and disease subtyping by clustering and survival analysis. PipeOne can be easily applied to other cancer types and complex diseases and extended to analyzing other types of high-throughput data. PipeOne is freely available at [https://github.com/nongbaoting/PipeOne](https://github.com/nongbaoting/PipeOne).

----

<img src = "./Figs/s0.png", width = "800">

#### Module 1: RNA-seq Processing

-----

<img src = "./Figs/s1.png", width = "800">


#### Module 2: Feature Prioritization

-----

<img src = "./Figs/s2.png">


#### Module 3: Subtype Analysis

----------

<img src = "./Figs/s3.png", width = "800">