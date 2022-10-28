
function [segm_vol_full, medic_img_struct]  = load_hnscc_segm(segm_folder, patient_name, segm_type)
% find the desired segmentation file based on folder names

%%% ARGUMENTS
% data_folder           : string containing the data folder path
% patient_name          : string containing the patient identification name (ex: "HNSCC5")
% segm_type             : string containing the segmentation type

% find segmentation file
switch segm_type
    case 'tumor'
        pn = char(patient_name);
        tumor_folder = fullfile(segm_folder,patient_name,'ROI mask', ['Hnscc',pn(6:end)]);
    otherwise
        error("Only tumor segmentations exist for HNSCC.");
end

% segmentation volume
segm_vol_full = read_DICOM_from_dir(tumor_folder, []);
segm_vol_full = permute(segm_vol_full, [2 1 3]);

% segmentation header
dname=dir(fullfile(tumor_folder,'IM*'));
medic_img_struct = dicominfo(fullfile(tumor_folder,dname(1).name));

% Cut tumor if too big
% % 1. Select only the major tumor part if splitted along Z axis,
% % 2. Select only 20 middles slices if bigger than that
segm_vol_full = cut_if_big_object(segm_vol_full);

end
