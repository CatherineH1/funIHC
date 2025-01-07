# Addressing class imbalance in functional data clustering
The goal of functional clustering is twofold: first, to categorize curves with similar temporal behaviors into separate clusters, and second, to obtain a representative curve that summarizes the typical temporal behavior within each cluster. An important challenge in current functional clustering techniques is class imbalance, where some clusters contain a significantly greater number of curves than others. While class imbalance is extensively addressed in supervised classification, it remains relatively unexplored in unsupervised contexts. To address this gap, we propose adapting the iterative hierarchical clustering approach, originally designed for multivariate data, to the context of functional data. Thus introducing a novel method called functional iterative hierarchical clustering (funIHC) to effectively handle the clustering of imbalanced functional data. Through comprehensive simulation studies and benchmarking datasets, we demonstrate the effectiveness of the funIHC approach. Utilizing funIHC on gene expression data related to human influenza infection induced by the H3N2 virus, we identify five distinct and biologically meaningful patterns of gene expression. 

Higgins, C., Carey, M. Addressing class imbalance in functional data clustering. Adv Data Anal Classif (2024). https://doi.org/10.1007/s11634-024-00611-8


## Overview
This repository contains R and MATLAB implementations of funIHC.

An example demonstrating how to use both the R and MATLAB versions, along with the data required to reproduce the simulations described in our paper, is provided.

## Data
The Data.zip file contains the datasets required to reproduce one iteration of the simulations descirbed in Section 3.1 and 3.2.

## Code
funIHC.R contains the R implementation, and example.R provides an example of running funIHC

funIHC.m contains the MATLAB implementation, and example.m provides and example of running funIHC

## Code options
As described in our paper, we propose three distinct distance metrics for use with
funIHC:
1. funIHC (curves): where the functional form of the curves is directly utilized with
a measure of functional distance to quantify dissimilarity between them. Set type=0
2. funIHC (1st derivative): the first-order derivative of the curves is utilized with
a measure of functional distance to quantify dissimilarity between them. Set type=1
3. funIHC (coefficients): which involves a dimension reduction via a basis function
expansion of the original curves followed by computation of distances based on
the resulting vectors. Set type=2
