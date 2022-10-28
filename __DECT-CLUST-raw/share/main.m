%
% Curve clustering with the MixFRHLP model and the EM (or a CEM-like)
% algorithm
%
%%%% Faicel Chamroukhi (2011)%%%%%%%
%
%   When using this code please cite the following papers : The two first
%   ones concern the model and its use in clusterng and the two last ones
%   concern the model and its use in discrimination, and a review.
%
% @article{Chamroukhi-RHLP-2009,
% 	Author = {Chamroukhi, F. and Sam\'{e}, A. and Govaert, G. and Aknin, P.},
% 	Journal = {Neural Networks},
% 	Number = {5-6},
% 	Pages = {593--602},
% 	Publisher = {Elsevier Science Ltd.},
% 	Title = {Time series modeling by a regression approach based on a latent process},
% 	Volume = {22},
% 	Year = {2009}
%     }
%
% @article{Chamroukhi-MixRHLP-2011,
% 	Author = {Sam{\'e}, A. and Chamroukhi, F. and Govaert, G{\'e}rard and Aknin, P.},
% 	Issue = 4,
% 	Journal = {Advances in Data Analysis and Classification},
% 	Pages = {301--321},
% 	Publisher = {Springer Berlin / Heidelberg},
% 	Title = {Model-based clustering and segmentation of time series with changes in regime},
% 	Volume = 5,
% 	Year = {2011}
%     }
%
%
%
% @article{Chamroukhi-FDA-2018,
% 	Journal = {Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery},
% 	Author = {Faicel Chamroukhi and Hien D. Nguyen},
% 	Note = {DOI: 10.1002/widm.1298.},
% 	Volume = {},
% 	Title = {Model-Based Clustering and Classification of Functional Data},
% 	Year = {2018},
% 	Month = {Dec}}
%
%
% @article{Chamroukhi-RHLP-FLDA,
% 	Author = {Chamroukhi, F. and Sam\'{e}, A. and Govaert, G. and Aknin, P.},
% 	Journal = {Neurocomputing},
% 	Number = {7-9},
% 	Pages = {1210--1221},
% 	Title = {A hidden process regression model for functional data description. Application to curve discrimination},
% 	Volume = {73},
% 	Year = {2010}
%     }
%
% @article{Chamroukhi-FMDA-2013,
% 	Author = {Chamroukhi, F. and Glotin, H. and Sam{\'e}, A.},
% 	Journal = {Neurocomputing},
% 	Pages = {153-163},
% 	Title = {Model-based functional mixture discriminant analysis with hidden process regression for curve classification},
% 	Volume = {112},
% 	Year = {2013}
%     }


clear;
close all;
clc;


G = 10;% number of clusters
K = 3;% number of regimes (polynomial regression components)
p = 2;% degree of the polynomials
q = 1;% order of the logistic regression (by default 1 for contiguous segmentation)


load('data/subject8_tumor_GT.mat');
load('data/subject8_tumor.mat');

% segm_vol_full: ....?
% subject: ....?   

% Constructing the data matrix structures
[nrows, ncols, t] = size(segm_vol_full);
nlevels = length(subject);

nviews = size(subject{1}, 3);

view = 1;

%decay_curves
for i=1:nrows
    for j=1:ncols
        for kev=1:nlevels
            decay_curves(i,j,kev) = subject{kev}(i,j,view);
        end
        plot(squeeze(decay_curves(i,j,:)))
        pause
        clf
    end
end
% voxel coordinates.

% Energy curves

% Select a ROI on one 2D slice of the 3D volume
slic_min = 115; % for subject8 (lower slice containing a tumor)
slic_nb = 116; % select a number between 116 and 123 (slices with tumor)
[col_obj,row_obj] = find(imdilate(segm_vol_full(:,:,slic_nb),strel('disk',30)));  % ROI around tumor
col_obj = col_obj + 10; % small shift for subject8

% Store at line 'i' the corresponding coordinates and curves
coord = [col_obj,row_obj];
decay_curves = zeros(length(col_obj),21);
for r=1:length(col_obj)
    for kev=1:21
        decay_curves(r,kev) = subject{kev}(row_obj(r),col_obj(r),slic_nb-slic_min);
    end
end

Y = decay_curves;
T = 1:21;

%type_variance = 'common';
type_variance = 'free';
n_tries = 2;
max_iter = 1000;
init_kmeans = 1;
threshold = 1e-5;
verbose = 1;
verbose_IRLS = 0;

[n, m]=size(Y);


t0 = tic;
solution =  MixFRHLP_EM(T, Y, G , K, p, q, type_variance, init_kmeans, n_tries, max_iter, threshold, verbose, verbose_IRLS);
%solution =  MixFRHLP_CEM(Y, G , K, p, q, type_variance, init_kmeans, n_tries, max_iter, threshold, verbose, verbose_IRLS);
disp(['Elapsed time: ',num2str(toc(t0))])

% bic(K)=solution.BIC;
% aic(K)=solution.AIC;
% icl(K)=solution.ICL1;

show_MixRHLP_results(Y, solution)


% show results on slice
cl_num = G;
clr = jet(cl_num);
slic = map_to_0_1(subject{1}(:,:,1))';
im_tum = cat(3, slic, slic, slic);

figure;
imshow(im_tum)
title("slice " + num2str(slic_nb-1)) % index showed on 3DSlicer
hold on
for cl_id=1:cl_num
    
    ind = solution.klas==cl_id;
    plot(row_obj(ind),col_obj(ind),'.', 'color',clr(cl_id,:),'MarkerSize',5);
    
end

