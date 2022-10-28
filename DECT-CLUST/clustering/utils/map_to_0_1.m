
function img = map_to_0_1(img)

img = img + 1024;  % subtract min value = -1024
img = double(img) / (3071+1024);  % divide by max value = 3071

end

