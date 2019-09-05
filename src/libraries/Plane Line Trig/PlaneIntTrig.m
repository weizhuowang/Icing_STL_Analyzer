function [Intersect,IntP] = PlaneIntTrig(n,p0,pts)
%PlaneIntTrig Plane intersection with Triangle
%   n=normal of plane, p0=point plane passing, p1,p2,p3=three points defining
%   the triangle

p1 = pts(1,:);
p2 = pts(2,:);
p3 = pts(3,:);

l1 = [p1;p2];
l2 = [p2;p3];
l3 = [p3;p1];

[tf1,Ip1] = PlaneIntLine(n,p0,l1);
[tf2,Ip2] = PlaneIntLine(n,p0,l2);
[tf3,Ip3] = PlaneIntLine(n,p0,l3);

Intersect = tf1||tf2||tf3;
if Intersect
    Ips = [Ip1;Ip2;Ip3];
    Ips = DeleteRepeat(Ips);
    IntP = Ips;
else
    IntP = [];
end

% % Debug
% fprintf('Intersect: %d | ',[tf1,tf2,tf3])
% fprintf('\nTrig_intersect: %d\n',Intersect)
% 
% line = [p1;Ip1;p2;Ip2;p3;Ip3;p1];
% x = line(:,1);
% y = line(:,2);
% z = line(:,3);
% plot3(x,y,z);grid on;hold on;
% scatter3(x,y,z)
% xlabel('x');ylabel('y');zlabel('z');
% 
% fprintf('\nIp: %2.2f %2.2f %2.2f\n',Ip1)
% fprintf('Ip: %2.2f %2.2f %2.2f\n',Ip2)
% fprintf('Ip: %2.2f %2.2f %2.2f\n',Ip3)

end

function [NoRepeatIps] = DeleteRepeat(Ips)

    NoRepeatIps = Ips;
    for i = 1:size(Ips,1)
        pt1 = Ips(i,:);
        for j = i+1:size(Ips,1)
            pt2 = Ips(j,:);
            theta = norm(cross(pt1,pt2));
            if (theta < 1e-10) && (norm(pt1-pt2)< 1e-10)
                NoRepeatIps(j,:) = [];
            end
        end
    end
%     NoRepeatIps

end
