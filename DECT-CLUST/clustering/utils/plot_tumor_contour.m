
% plot tumor contour
function plot_tumor_contour(tumor_contour, xy_min, clr)
if ~isempty(tumor_contour)
    for c=1:length(tumor_contour)
        tum_cont = cell2mat(tumor_contour(c));
        col_tum = tum_cont(:,1)-xy_min(1)+1;
        row_tum = tum_cont(:,2)-xy_min(2)+1;
        for i=1:length(row_tum) % tumor in blue
            plot(row_tum(i),col_tum(i),'.','color',clr,'MarkerSize',6);  % light blue: [0,0.5,1], white: [0.99,0.99,0.99]
        end
    end
end
end