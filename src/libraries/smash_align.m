function [sma_slice_cln,sma_slice_ice] = smash_align(p0_total,n_total,slc_cln,slc_ice,display_fig)

    [sma_slice_cln] = smash_slice(p0_total,n_total,slc_cln,display_fig);
    [sma_slice_ice] = smash_slice(p0_total,n_total,slc_ice,display_fig);

    % Take longest segment (ignore small ice cones)
    for slc_idx = 1:size(sma_slice_ice,2)
        fprintf('%d \n',slc_idx)
        seg = segmentation(sma_slice_ice{slc_idx})';
        sma_slice_ice{slc_idx} = pts2Ixy(seg{1});
        
        seg = segmentation(sma_slice_cln{slc_idx})';
        sma_slice_cln{slc_idx} = pts2Ixy(seg{1});
    end

    % Align all pairs of ice and clean slice
    for i = 1:size(sma_slice_cln,2)
        fprintf('%d \n',i)
        sma_slice_cln{i} = shift_slice_2d(sma_slice_cln{i},sma_slice_ice{i});
    end

end

function [Ixy] = pts2Ixy(seg)
    
    Ixy = [];
    Ixy.Ix_s = [seg(:,1),[seg(end,1);seg(1:end-1,1)]];
    Ixy.Iy_s = [seg(:,2),[seg(end,2);seg(1:end-1,2)]];

end