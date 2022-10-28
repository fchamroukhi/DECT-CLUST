
function [kev_list_Mono, kev_list_OM] = get_kev_lists(data_info)
    % get the energy values list (kev)  ## 2 records of energy are sometimes in files: 'Mono' and maybe something else
    kev_list_Mono = [];
    kev_list_OM = [];

    for energy_idx=1:length(data_info{1}.study{1}.series)
        descr = split(data_info{1}.study{1}.series{energy_idx}.info.SeriesDescription);
        for i=1:length(descr)
            k = str2double(descr{i});
            if ~isnan(k)
                if strcmp(descr{1},'Mono')
                    kev_list_Mono = [kev_list_Mono, k];
                else
                    kev_list_OM = [kev_list_OM, k];
                end
            end
        end
    end

    kev_list_Mono = sort(kev_list_Mono);
    kev_list_OM = sort(kev_list_OM);

end