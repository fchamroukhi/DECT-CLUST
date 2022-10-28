
% Cut ground truth object if too big:
% % 1. Select only the major object part if splitted along Z axis,
% % 2. Select only 20 middles slices if bigger than that


function segm_vol_full = cut_if_big_object(segm_vol_full)

    [~,~,slic_z] = ind2sub(size(segm_vol_full),find(segm_vol_full));
    slic_z = unique(slic_z);
    
    
    %%% If object is in several parts with discontinuities on Z dim,
    %%% select only the biggest part of Z slices
    
    ds = diff(slic_z);
    ch = find(ds-1);  % indices of discontinuities
    if ~isempty(ch)
        ch = [0; ch; length(slic_z)];
        [~,ind] = max(ch(2:end)-ch(1:end-1));
        selected_slices = slic_z((ch(ind)+1):ch(ind+1));
        slic_z = selected_slices;
    end


    %%%% if object is too big, take only middle slices

    if length(slic_z) > 20
        mid = ceil(length(slic_z)/2);
        selected_slices = slic_z(mid-10 : mid+9); % 20 middles slices
    end


    %%% Replace by 0 all other parts of object if necessary

    if exist('selected_slices','var')
        segm_vol = segm_vol_full(:,:,selected_slices);
        segm_vol_full = zeros(size(segm_vol_full)); 
        segm_vol_full(:,:,selected_slices) = segm_vol;
    end

end
