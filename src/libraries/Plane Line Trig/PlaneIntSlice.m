% slice is in combined format
function [IntP] = PlaneIntSlice(n,p0,slice_comb)

    if size(slice_comb,3) == 2
        slice_comb(:,:,3) = [0];
    end
    
    IntP = [];
    for i = 1:size(slice_comb,1)
%         pts = squeeze(slice_comb(i,:,:))
        pts = [slice_comb(i,1,1),slice_comb(i,1,2),slice_comb(i,1,3);slice_comb(i,2,1),slice_comb(i,2,2),slice_comb(i,2,3)];
        [Intersect,IP] = PlaneIntLine(n,p0,pts);
        if Intersect
            if IP(1) <= p0(1)
                IntP(end+1,:) = IP;
            end
        end
    end

end