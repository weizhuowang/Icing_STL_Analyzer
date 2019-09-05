function [segment] = segmentation(slice)
%Segmentation take in slice, put out segments
%   Input slice in Ixy format, but not in order
%   Slice contains numerous amount of 2 point segments, segment returns
%   connect lines in cells of pts format

% Currently 2D only, might need some mods if 3D

    total = zeros([size(slice.Ix_s),2]);
    total(:,:,1) = slice.Ix_s;
    total(:,:,2) = slice.Iy_s;
    
%     total = total(1:100,:,:); %
    for i = 1:size(total,1)
        segment{i} = squeeze(total(i,:,:));
    end
    
    prev_sz = 0;
    i=0;
    
    counter = 5;
    while counter > 0
        if abs(size(segment,2)-prev_sz ) > 0
            counter = 5;
        else
            counter = counter - 1;
        end
%     for j = 1:600
        i = i+1;
        prev_sz = size(segment,2);
        [segment] = iter_combine(segment);
        
        [~,I] = sort(cellfun(@length,segment));
        segment = segment(flip(I));
    
        if mod(i,2) == 1
            fprintf('Iteration %d, current_size: %d \n',i,prev_sz)
        end
    end
        
end