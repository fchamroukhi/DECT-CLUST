
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                             Compile results                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segolene Brivet (segolene.brivet@mail.mcgill.ca)


% The general workflow is:
% 1. main.m: learn model and get clustering results in 'mixstats_red' variable 
% 2. build_cluster_separ_idx: compute clustering index scores and Dice scores for several images and update 'mixstats_red' variable with these scores
% 3.Here: build_results_table: compile results for several methods and save each of them as a .mat file
% 4. resutls_analysis: boxplot and statistical values to analyse methods


results_folder_name = 'results-SgMFR-Bspl';
% results_folder_name = 'results-SgMFR-poly';
% results_folder_name = 'results-SgMVFR-Bspl';
% results_folder_name = 'results-SsMFR';
% results_folder_name = 'results-GMM';
% results_folder_name = 'results-kmeans';
% results_folder_name = 'results-SelSearch';

pat_results = {};
nb_K = [];
time_results = [];

dice_results = [];
DB_spat_results = [];
DB_spec_results = [];
DBt_spat_results = [];
DBt_spec_results = [];




% patient_names = ["HNSCC21"];
patient_names = ["HNSCC2","HNSCC3","HNSCC5","HNSCC8","HNSCC9","HNSCC10",...
    "HNSCC11","HNSCC12","HNSCC13","HNSCC15","HNSCC15A","HNSCC17","HNSCC17A","HNSCC18","HNSCC20",...
    "HNSCC21","HNSCC22","HNSCC22A","HNSCC25","HNSCC26","HNSCC27","HNSCC29","HNSCC30",...
    "HNSCC31A","HNSCC32","HNSCC33","HNSCC34","HNSCC35","HNSCC36","HNSCC37A","HNSCC38","HNSCC39",...
    "HNSCC41","HNSCC42","HNSCC44","HNSCC44AM","HNSCC45","HNSCC46","HNSCC47","HNSCC48","HNSCC49",...
    "HNSCC51","HNSCC52","HNSCC52AM","HNSCC53","HNSCC55","HNSCC56","HNSCC57",...
    "HNSCC61A","HNSCC62","HNSCC63","HNSCC63A","HNSCC64A","HNSCC65A","HNSCC66","HNSCC67","HNSCC68","HNSCC69",...
    "HNSCC70A","HNSCC71","HNSCC72A","HNSCC73","HNSCC74","HNSCC75","HNSCC76","HNSCC77","HNSCC78","HNSCC79","HNSCC80",...
    "HNSCC81","HNSCC82","HNSCC83","HNSCC84","HNSCC85","HNSCC87","HNSCC88","HNSCC89","HNSCC90",...
    "HNSCC91","HNSCC92","HNSCC95","HNSCC96","HNSCC97","HNSCC98",...
    "HNSCC100","HNSCC101","HNSCC103","HNSCC105","HNSCC106","HNSCC108","HNSCC109"];

for patient_name = patient_names
    
    lambda = 0.075;
    roi_size = 150;
    
    fls = dir(fullfile(results_folder_name,[char(patient_name),'_*ROI',num2str(roi_size),'_*',num2str(lambda),'_mixstats_red.mat']));
    if length(fls) > 1
        disp("Several results are in folder for this case");
    elseif length(fls) == 0
        nb_K = [nb_K, NaN];
        dice_results = [dice_results, NaN];
        time_results = [time_results, NaN];
        DB_spat_results = [DB_spat_results, NaN];
        DB_spec_results = [DB_spec_results, NaN];
        DBt_spat_results = [DBt_spat_results, NaN];
        DBt_spec_results = [DBt_spec_results, NaN];
    else
        mixstatspath = fullfile(results_folder_name,fls(1).name);
        load(mixstatspath);
        
        nb_K = [nb_K, length(unique(mixstats_red.klas))];  % for MFR, MVFR, GMM data format
        % nb_K = [nb_K, length(unique(mixstats_red.T))];   % for kmeans and SelSearch data format.
        
        dice_results = [dice_results, mixstats_red.dice];
        time_results = [time_results, mixstats_red.elaps_time];
        
        DB_spat_results = [DB_spat_results, mixstats_red.DB_spat];
        DB_spec_results = [DB_spec_results, mixstats_red.DB_spec];
        DBt_spat_results = [DBt_spat_results, mixstats_red.DBt_spat];
        DBt_spec_results = [DBt_spec_results, mixstats_red.DBt_spec];
    end
    
end

results.nb_K = nb_K;
results.dice = dice_results;
results.time = time_results;
results.DB_spat = DB_spat_results;
results.DB_spec = DB_spec_results;
results.DBt_spat = DBt_spat_results;
results.DBt_spec = DBt_spec_results;
results.patient_names = patient_names;

% save('results_SgMFR-Bspl','results')
% save('results_SgMFR-poly','results')
% save('results_SgMVFR-Bspl','results')
% save('results-SsMFR-Bspl','results')
% save('results_GMM','results')
% save('results_kmeans','results')
% save('results_SelSearch','results')

mean_nbK = nanmean(nb_K)
std_nbK = nanstd(nb_K)
mean_dice = nanmean(dice_results)
std_dice = nanstd(dice_results)
mean_time = nanmean(time_results)
std_time = nanstd(time_results)


