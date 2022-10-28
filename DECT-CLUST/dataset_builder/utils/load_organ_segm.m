
function [nrrd_struct, segm_vol_full, segm_vol_full_list, study_type] = ...
    load_organ_segm(segm_folder, patient_name, segm_type, organs_id)
% find the desired segmentation file and its corresponding study based on folder names

%%% ARGUMENTS
% data_folder           : string containing the data folder path
% patient_name          : string containing the patient identification name (ex: "SRG8_MultiEnergy")
% segm_type             : string containing the segmentation type
% organ_id              : cell of cells containing the organ ids to focus on, {} for tumor or lymph node

% legend for organ_id:
% 
% 1 - LeftSubmandibularGland 
% 2 - RightSubmandibularGland 
% 3 - LeftParotidGland
% 4 - RightParotidGland
% 5 - ThyroidGland
% 6 - ThyroidCartilage
% 7 - LungApex 
% 8 - MuscleTissue 
% 9 - C3VertebralBody 


% find segmentation file
filelist = dir(fullfile(segm_folder,patient_name,'**/*.nrrd'));
switch segm_type
    case 'tumor'
        organ_pos = regexp({filelist.name},'Primary ?[Tt]umor-label.nrrd');
    case 'lymphnode'
        organ_pos = regexp({filelist.name},'Lymph ?node-label.nrrd');
    case 'organs'
        organ_pos = regexp({filelist.name},'[Oo]rgans-label.nrrd');
    otherwise
        error("Unrecognized segmentation type.");
end
organ_ind = find(~cellfun(@isempty,organ_pos));
if isempty(organ_ind)
    error("This patient has no file of " + segm_type + " segmentation in directory or subdir of " + fullfile(segm_folder,patient_name));
end
org_i = organ_ind(1);
organ_ind(1) = [];
organ_folder = filelist(org_i).folder; % take the first to appear by default if multiple files
organ_file_name = fullfile(organ_folder,filelist(org_i).name);

% match segmentation with DECT study
while 1
    if ~isempty(regexp(organ_folder, '[Oo]pened[ _]mouth'))
        study_type = 'OM';
        break;
    elseif ~isempty(regexp(organ_folder, '[Tt]ongue[ _]out'))
        study_type = 'OM';
        % try to avoid this case and find another segmentation at next iteration: (not sure this is OM study each time)
        if ~isempty(organ_ind)
            org_i = organ_ind(1);
            organ_ind(1) = [];
            organ_folder = filelist(org_i).folder; % take the first to appear by default if multiple files
            organ_file_name = fullfile(organ_folder,filelist(org_i).name);
        else
            break;
        end
    else % i.e. nothing on folder name (default) or "closed mouth"
        study_type = 'Mono';
        break;
    end
end

% segmentation volume
nrrd_struct = nhdr_nrrd_read(organ_file_name, 1);

% find object of interest if organs
if strcmp(segm_type,"organs")
    segm_vol_full = [];
    segm_vol_full_list = {};
    for c=1:length(organs_id)
        org_cell = organs_id{c};
        org_id = org_cell(1);
        segm_focus = nrrd_struct.data == org_id{1};
        for i=2:length(org_cell)
            org_id = org_cell(i);
            segm_focus = segm_focus + (nrrd_struct.data == org_id{1});
        end
        segm_vol_full_list{c} = permute(int8(segm_focus), [2 1 3]);
    end
else
    segm_vol_full = permute(nrrd_struct.data, [2 1 3]);
    segm_vol_full_list = {};
end

nrrd_struct.data = [];



% % % % % % % % % %%%%%  REFINEMENTS IF TOO BIG  %%%%% % % % % % % % % % %

% Cut volume section if too big
% % 1. Select only the major tumor part if splitted along Z axis,
% % 2. Select only 20 middles slices if bigger than that

if isempty(segm_vol_full_list) && ~isempty(segm_vol_full)
    segm_vol_full = cut_if_big_object(segm_vol_full);
else
    for c = 1:length(segm_vol_full_list)
        segm_vol_full_list{c} = cut_if_big_object(segm_vol_full_list{c});
    end
end




end




