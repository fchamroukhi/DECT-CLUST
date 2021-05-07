SFRM

Spatial Functional Regression Mixture Models for t+2D/t+3d image Segmentation

Includes at least the following strategies:

Strategy 1 (FunSRM): A mixture of functional i-polynomial, ii-Spline, or ii-B-Spline regressions, with spatial-depending mixing proportions.
Strategy 2 (ExtraFeat-GMM): The 'temporal' curves of the image are first represented by A mixture of functional i-polynomial, ii-Spline, or ii-B-Spline regressions. A mixture of multivariate Gaussian mixture, with spatial-depending mixing proportions, is then applied to the regression parameters.

- The spatial dependency can be modelled via a softmax function, or a normalised weighted Gaussian functions, given the spatial voxel coordinates

- The statistical inference is performed via dedicated expectation-maximisation (EM) algorithms for conditional maximum likelihood or joint maximum likelihood

- An additionnel stage of Markov Random Field (MRF) smoothing can be added as an afterwards step to each of these strategies.

