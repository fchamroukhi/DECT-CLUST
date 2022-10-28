function [select_ind, vrb] = adjust_GT_subj_slices(sz_gt, sz_subj)

if sz_subj < sz_gt
    ind = linspace(1,sz_gt,sz_subj);
    select_ind = round(ind);
	vrb = 'gt';
	
else  % sz_gt <= sz_subj
    ind = linspace(1,sz_subj,sz_gt);
    select_ind = round(ind);
	vrb = 'subj';
end