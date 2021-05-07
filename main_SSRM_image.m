%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Learning spatial mixture of functional regression models for spectral image clustering%% FC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear all; close all; clc;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          %% choose a regression type %%                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%model = 'PRM'; % Polynomial regression mixturemodel = 'SRM'; % Spline regression mixturemodel = 'bSRM';% B-Spline regression mixturemodel_selection = 0;%% data (spectral image)load('data/subject8_tumor_GT.mat');load('data/subject8_tumor.mat');% Select a ROI on one 2D slice of the 3D volumeslic_min = 115; % for subject8 (lower slice containing a tumor)slic_nb = 116; % select a number between 116 and 123 (slices with tumor)%[row_obj, col_obj] = find(imdilate(segm_vol_full(100:200,200:350,slic_nb)',strel('disk',30)));  % ROI around tumor[col_obj, row_obj] = find(imdilate(segm_vol_full(:,:,slic_nb),strel('disk',30)));  % ROI around tumor%col_obj = col_obj + 10; % small shift for subject8% Store at line 'i' the corresponding coordinates and curvescoord = [col_obj, row_obj];decay_curves = zeros(length(row_obj),21);for r=1:length(col_obj)    for kev=1:21        decay_curves(r,kev) = subject{kev}(row_obj(r),col_obj(r),slic_nb-slic_min);    endend% Decay_curves = zeros(length(col_obj), length(row_obj), 21);% for i=1:length(col_obj)% for j=1:length(row_obj)%     for kev=1:21%         Decay_curves(i,j,kev) = subject{kev}(row_obj(i),col_obj(j),slic_nb-slic_min);%     end% end% end% coord = grid2d(512, 512);% Y = reshape(Decay_curves,[],21);%Y = decay_curves;V = coord/max(max(coord)); % scale coordinates to keep them in [0,1]T = 1:21;T = T/max(max(T)); % scale x sampling values to keep them in [0,1] %linspace(0, 1, m);% % % %% Uncomment to apply to Other (non-spatial) curve data sets (for algo testing)%   dataname = 'waveform'; load(dataname); Y = waveform;% %  dataname = 'satellite';load(dataname); Y = satellite;% % % dataname = 'yeast_cellcycle';load(dataname); Y = yeast_cellcycle;% % %  dataname = 'phonemes'; load(dataname); Y = phonemes;%   [n, m] = size(Y); T = linspace(0,1, m); V = ones(length(Y), 2);%not spatial data%%[n, m] = size(Y);%Y = Y - ones(n,1)*mean(Y,1);%% data matrices% Spatial coordinatesCurves.spatialcoord = V;% data.VoxelCoordinates = V;% CurvesCurves.abscissas = T;% data.WavelengthLevels = T;Curves.ordinates =  Y;% data.ReflectanceValues = Y; %% SRM model specificationK            = 20; % number of clusters in the data% algo setting nbr_EM_runs = 1;switch(model)    case('PRM')        p = 4; % polynomial regression degree                regressionOptions.basis = 'polynomial';        regressionOptions.p = p;    case 'SRM'        spline_order = 4; % 2: linear, 3: qudratic, 4: cubic, etc, for (b-)spline spatial regression mixture        nknots       = 10; % fixed number of internal knots, for (b-)spline spatial regression mixture                regressionOptions.basis = 'spline';        regressionOptions.spline_order = spline_order;        regressionOptions.nknots = nknots;    case('bSRM')        Bspline_order = 4; % 2: linear, 3: qudratic, 4: cubic, etc, for (b-)spline spatial regression mixture        nknots       = 10; % fixed number of internal knots, for (b-)spline spatial regression mixture                regressionOptions.basis = 'B-spline';        regressionOptions.Bspline_order = Bspline_order;        regressionOptions.nknots = nknots;    otherwise        error('unknown model type');end%% SRM Model fittingtic;if (~model_selection) % no model selection (fixed K)    [mixModel, mixStats] = learn_SRM_EM(Curves, K, mixingOption, regressionOptions, nbr_EM_runs);else % Select K from the data    current_BIC = -inf;    Kmin = 5; Kmax = 30;    Krange = Kmin:Kmax;        BIC = zeros(1, length(Krange));    for K = Krange(1):Krange(end)        [mixModel_K, mixStats_K] = learn_SRM_EM(Curves, K, mixingOption, regressionOptions, nbr_EM_runs);        fprintf(1,'Number of segments K: %d | BIC %f \n', K, mixStats_K.BIC);        if mixStats_K.BIC>current_BIC            mixModel = mixModel_K;            mixStats = mixStats_K;            %            current_BIC = mixStats_K.BIC;        end        BIC(K - Krange(1)+1) = mixStats_K.BIC;    endendfprintf('Elapsed time %f sec \n', toc);    %% plot clustering resultsshow_SRM_results(Curves, mixModel, mixStats);