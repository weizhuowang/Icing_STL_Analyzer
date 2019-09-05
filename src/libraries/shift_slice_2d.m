% shift the clean slice to match the icing one
function [slice_cln] = shift_slice_2d(slice_cln,slice_ice)

    slice_cln_comb = [];
    slice_ice_comb = [];
    
    slice_cln_comb(:,:,1) = slice_cln.Ix_s;
    slice_cln_comb(:,:,2) = slice_cln.Iy_s;
    slice_cln_comb(:,:,3) = zeros(size(slice_cln.Iy_s));
    
    slice_ice_comb(:,:,1) = slice_ice.Ix_s;
    slice_ice_comb(:,:,2) = slice_ice.Iy_s;
    slice_ice_comb(:,:,3) = zeros(size(slice_ice.Iy_s));

    align_pos = max(max(slice_ice.Iy_s))-0.02;
    corr = slice_diff(slice_cln_comb,slice_ice_comb,[0,1,0],[0,align_pos,0]);
    fprintf('| %2.2f |',corr);fprintf('\n');

    % Shift the slice to fit iceing slice
    slice_cln.Ix_s = slice_cln.Ix_s+corr(1);
    slice_cln.Iy_s = slice_cln.Iy_s+corr(2);
    
end

function correction = slice_diff(slice_cln,slice_ice,n,p0)

    slice_ana = slice_cln;
    base_pt = [];
    for i = 1:size(slice_ana,1)
        pts = squeeze(slice_ana(i,:,:));
        [Intersect,IntP] = PlaneIntLine(n,p0,pts);
        if Intersect
            base_pt = IntP;
        end
    end
    
    slice_ana = slice_ice;
    shifted_pts = [];
    for i = 1:size(slice_ana,1)
        pts = squeeze(slice_ana(i,:,:));
        [Intersect,IntP] = PlaneIntLine(n,p0,pts);
        if Intersect
            shifted_pts(end+1,:) = IntP;
        end
    end
    
    shifted_pts;
    
    correction = [0,0,0];
    % Closest
%     for i = 1:size(shifted_pts,1)
%         if (norm(correction,2) == 0)
%             correction = shifted_pts(i,:) - base_pt;
%         else
%             if norm(shifted_pts(i,:) - base_pt,2) < norm(correction,2)
%                 correction = shifted_pts(i,:) - base_pt;
%             end
%         end
%     end

    % Left most
    if length(shifted_pts) ~= 0
        min_idx = find(shifted_pts(:,1)==min(shifted_pts(:,1)));
        correction = shifted_pts(min_idx,:)-base_pt;
    else
        fprintf('No intersection, adjust line position\n')
    end
    
end

function [Intersect,IntP] = PlaneIntLine(n,p0,pts)
%PlaneIntLine Plane intersection with line
%   n=normal of plane, p0=point plane passing, p1,p2=two points defining
%   the line

    Intersect = false;
    IntP = [];

    p1 = pts(1,:);
    p2 = pts(2,:);

    d1 = n*(p1-p0)';
    d2 = n*(p2-p0)';

    % One point on the plane
    if d1*d2 == 0
        Intersect = true;
        if (d1 == 0) && (d2 == 0)
            IntP = [p1;p2];
        elseif d1 == 0
            IntP = p1;
        elseif d2 == 0 
            IntP = p2;
        end

    % Intersection in the middle
    elseif d1*d2 < 0
        Intersect = true;
        x = (p2-p1)/(norm(p2-p1));
        if n*x' == 0
            disp('Division by 0');
        else
            IntP = p2-x*(d2/(n*x'));
        end

    % No intersection Do nothing
    elseif d1*d2 > 0
    end

end