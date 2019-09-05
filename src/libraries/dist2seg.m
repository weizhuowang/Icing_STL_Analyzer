function [distance,shift] = dist2seg(p0,seg)

    % generate distance
    distance = [];
    if direct_above(p0,seg(1,:),seg(2,:))
        [Intersect,IntP] = PlaneIntLine(seg(2,:)-seg(1,:),p0,seg);
        if Intersect
            distance = norm(IntP-p0,2);
        end
    else
        dist1 = norm(p0-seg(1,:),2);
        dist2 = norm(p0-seg(2,:),2);
        distance = min(dist1,dist2);
    end
    
    % Generate shift
    seg_dir_unit = (seg(1,:)-seg(2,:))/norm(seg(2,:)-seg(1,:),2);
    shift = (p0-seg(2,:))*seg_dir_unit.';

end

% 1.1s for 2M calls
% return true if pt is directly above line seg
function abv_flg = direct_above(inq_pt,seg_pt1,seg_pt2)

    abv_flg = false;
    abv_flg = abv_flg || is_obtuse(seg_pt1,inq_pt,seg_pt2);
    abv_flg = abv_flg || is_obtuse(seg_pt2,inq_pt,seg_pt1);
    abv_flg = ~abv_flg;
    
end

% print out 1 when not (right/acute)
function obt_flg = is_obtuse(center_pt,pt1,pt2)

    leg1 = pt1 - center_pt;
    leg2 = pt2 - center_pt;
    obt_flg = leg1*leg2'<0;

end 