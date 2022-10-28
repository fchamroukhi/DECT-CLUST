# DECT-CLUST: 


DECT-CLUST contains an ensemble of original unsupervised learning algorithms dedicated to **D**ual-**E**nergy **C**omputed **T**omography image clustering and application to **H**ead and **N**eck **S**quamous **C**ell **C**arcinoma segmentation.

The archive contains the **source codes** used in the **DECT-CLUST paper**. In principle the data used in this context are not sharable but you can still reach out if you think we may have a common interest and we preserve data constraints...


The core codes are written by Faïcel Chamroukhi (faicel DOT chamroukhi AT math DOT cnrs DOT fr) and the routines written to perform the image experiments (including GPU components) are written by Ségolène Brivet


### Environment
This project is developed in Matlab. A remote GPU machine was used for heavy computations.


### Data
DECT are 4D data: a 3D body volume over a range of X-ray energy levels.  
*the data folder in this repo is empty and need to be filled on local machines.*  


### Models and Algorithms

DECT-CLUST is based on a family of Spatial Functional Regression Mixture Models for t+3d/t+2D image Segmentation.

It implements unsupervised learning techniques based on mixture models and functional data analysis (FDA) models to the clustering of DECT images. It designs functional mixture models that integrate spatial image context in mixture weights, with mixture component densities being constructed upon the DECT energy decay curves as functional observations. The developed algorithms are dedicated expectation-maximization (EM) algorithms for the maximum likelihood estimation (MLE) of the parameters of the constructed models. To our knowledge, this is the first article to develop statistical FDA and model-based clustering techniques to take advantage of the full spectral information provided by DECT. 




### 1. Image preparation
(Optional)
Specific to our data, we pre-process DICOM files in this step.  

In `dataset_builder` folder, a large section around tumors is cropped from DECT scans, as well as its ground truth segmentation, and the two are saved as .mat files in `data_tumor` folder.  


### 2. Clustering
In `clustering` folder, `main.m` script:
- loads a 4D-image file (.mat) and its associated 3D ground truth segmentation (.mat), 
- learn image clustering models with different tunable options,
- post-process clustering results and plot visualizations,

`build_cluster_separ_indx` script:
- computes metrics to assess the clustering quality,

`build_results_tables` script:
- compiles previously computed metrics from different methods,

`build_results_analysis` script:
- outputs and plots statistical comparative analysis (e.g. t-test, boxplots)

