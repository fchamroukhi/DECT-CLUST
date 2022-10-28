

function [img, focus, xy_min, xy_max] = mask_outOfBody(subj)

% INPUT: 
%           subj: 3D volume of an image (typically a 3D scan at one energy level)
% 
% OUTPUT:
%           img: 3D image volume of same size as 'subj' where background has been fixed to -1000
% 
% optional outputs:
%           foucs: 3D volume of same size as 'subj' where matrix entries are boolean: 1 if voxel is in foreground, 0 if voxel is in background
%           xy_min, xy_max: min and max coordinates that frame the foreground (right, top, left, bottom)


img = zeros(size(subj));

for z=1:size(subj,3)
    slic = subj(:,:,z);
    bigmask = bwareafilt( imfill( slic>-500 ,'holes') ,1);
    foc = -1000*ones(size(subj,1:2));   % fix background to -1000
    foc(bigmask) = slic(bigmask);   % populate matrix with relevant scan values 
    img(:,:,z) = foc;
end


if nargout > 1

    focus = zeros(size(subj));
    xcoord = []; ycoord = [];
    for z=1:size(subj,3)
        focus(:,:,z) = img(:,:,z) > -500;
        [x,y] = find(focus(:,:,z));
        xcoord = [xcoord;x]; ycoord = [ycoord;y];
    end
    xy_min = [min(xcoord), min(ycoord)]; xy_max = [max(xcoord), max(ycoord)];   % right, top, left, bottom of patient (when transposed img after)

end

end



