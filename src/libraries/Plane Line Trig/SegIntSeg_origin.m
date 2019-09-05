% %% test seg v seg
% 
% p1 = [0,0;
%       0.5,0.72];
% p2 = [0,1;
%       1,0];
% 
% figure(1);clf;
% plot(p1(:,1),p1(:,2),'r');grid on;hold on;
% plot(p2(:,1),p2(:,2)),'b';hold on;
% 
% tic
% for i = 1:1000000
%     [Int,Int_log] = SegIntSeg(p1, p2);
% end
% toc
% if Int_log
%     Int
%     scatter(Int(1),Int(2));hold on;
% end

% 100w 0.427s

function [Int,Int_log] = SegIntSeg(line1, line2)
    
    p1 = line1(1,:); % [x1,y1]
    q1 = line1(2,:); % [x2,y2]
    p2 = line2(1,:); % [x3,y3]
    q2 = line2(2,:); % [x4,y4]
    
    o1 = seg_orien(p1,q1,p2);
    o2 = seg_orien(p1,q1,q2);
    o3 = seg_orien(p2,q2,p1);
    o4 = seg_orien(p2,q2,q1);
    
    Int_log = false;
    if (o1 ~= o2 && o3 ~= o4) 
        Int_log = true; 
    end
    
    % p1, q1 and p2 are colinear and p2 lies on segment p1q1 
    if (o1 == 0 && onSegment(p1, p2, q1)) 
        Int_log = true; 
    end
  
    % p1, q1 and q2 are colinear and q2 lies on segment p1q1 
    if (o2 == 0 && onSegment(p1, q2, q1)) 
        Int_log = true; 
    end
    
    % p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) 
        Int_log = true; 
    end
    
    % p2, q2 and q1 are colinear and q1 lies on segment p2q2 
    if (o4 == 0 && onSegment(p2, q1, q2)) 
        Int_log = true; 
    end
    
    % Interception point
    Int = [];
    
% %     Could potentially fail if any of the seg is vertical
%     if Int_log
%         k1 = (q1(2)-p1(2))/(q1(1)-p1(1));
%         k2 = (q2(2)-p2(2))/(q2(1)-p2(1));
%         b1 = p1(2)-k1*p1(1);
%         b2 = p2(2)-k1*p2(1);
%         x = (b2-b1)/(k1-k2);
%         y = k1*x+b1;
%         Int = [x,y];
%     end

% Could potentially failed if two segs are colinear
    if Int_log
        u = (   (q2(2)-p2(2))*(p1(1)-p2(1)) - (q2(1)-p2(1))*(p1(2)-p2(2))   )/...
            (   (q2(1)-p2(1))*(q1(2)-p1(2)) - (q2(2)-p2(2))*(q1(1)-p1(1))   ) ;
        Int = [(q1(1)-p1(1))*u+p1(1) , (q1(2)-p1(2))*u+p1(2)];
    end
    
end


% Given three colinear points p, q, r, the function checks if 
% point q lies on line segment 'pr' 

function [onSeg] = onSegment(p,q,r) 

    onSeg = false;
    if (q(1) <= max(p(1), r(1)) && q(1) >= min(p(1), r(1)) && q(2) <= max(p(2), r(2)) && q(2) >= min(p(2), r(2))) 
        onSeg = true; 
    end

end
    
% To find orientation of ordered triplet (p, q, r). 
% The function returns following values 
% 0 --> p, q and r are colinear 
% 1 --> CW 
% 2 --> CCW
function [Orien] = seg_orien(p,q,r)

    val = (q(2) - p(2)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(2) - q(2));
    if val == 0
        Orien = 0;
    else
        if val > 0
            Orien = 1;
        else
            Orien = 2;
        end
    end

end
