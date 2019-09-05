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