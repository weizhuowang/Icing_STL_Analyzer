% Algorithm:
% Concave signed area by cross product

clc,clear;
addpath(genpath('..'))

pts = [-1,-1;
       -1 ,0.5 ;
       1 ,-1;
       1,-2]

figure(1);clf;
plot([pts(:,1);pts(1,1)],[pts(:,2);pts(1,2)]);grid on;hold on;
scatter(0,0);hold on;
xlim([min(pts(:,1)),max(pts(:,1))]*1.5);ylim([min(pts(:,2)),max(pts(:,2))]*1.5);

area = 0;
for i = 1:size(pts,1)
    if i ~= size(pts,1)
        i1 = i+1;
    else
        i1 = 1;
    end
    xi = pts(i,1);
    yi = pts(i,2);
    xi1 = pts(i1,1);
    yi1 = pts(i1,2);
    area = area + 0.5*(xi*yi1 - xi1*yi);
end

area = abs(area)




%% test seg v seg

p1 = [0.4,0;
      0.4,0.7];
p2 = [0,1;
      1,0];

figure(1);clf;
plot(p1(:,1),p1(:,2),'r');grid on;hold on;
plot(p2(:,1),p2(:,2)),'b';hold on;

tic
for i = 1:1000000
    [Int,Int_log] = SegIntSeg(p1, p2);
end
toc
if Int_log
    Int
    scatter(Int(1),Int(2));hold on;
end