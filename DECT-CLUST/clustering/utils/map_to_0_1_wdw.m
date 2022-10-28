
function img = map_to_0_1_wdw(img,mnx)

img(img<mnx(1)) = mnx(1);
img(img>mnx(2)) = mnx(2);

img = img - mnx(1);  % subtract min value
img = double(img) / (mnx(2)-mnx(1));  % divide by max value

end