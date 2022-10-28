
function [char_subj, subject_nb, subject_id, org_str_list] = get_subject_names(patient_name, segm_type, organs_id)

%%% patient_name: string for patient name (like "HNSCC5")

org_str_list = {};

char_subj = char(patient_name);
if strcmp(char_subj(1:5),'HNSCC')
    subject_nb = char_subj(6:end);
    subject_id = char(patient_name);
else
    subject_nb = char(regexp(patient_name,"\d*",'match'));
    subject_id = ['subject',char(subject_nb),'_',segm_type];
    if ~isempty(organs_id)
        for c=1:length(organs_id)
            org_cell = organs_id{c};
            org_id = org_cell(1);
            org_str = num2str(org_id{1});
            for i=2:length(org_cell)
                org_id = org_cell(i);
                org_str = [org_str,'-',num2str(org_id{1})];
            end
            org_str_list{c} = org_str;
        end
    end
end

end