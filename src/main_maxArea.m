clc,clear;
% =========Slicer Set up==========
addpath(genpath('..'))
p0_total = []; n_total = [];
tot_idx = 1:1:3;
for i = tot_idx
% for i = 1
    p0_total(end+1,:) = [0,i,0];
    n_total(end+1, :) = [0,1,0];
end

display_fig = false;
save_video = true;

% stls = ["Scallop-Max_Redo Root to Mid 200A stretch.stl",...
%         "Scallop-Max_Redo Root to Mid 200B stretch.stl",...
%         "Scallop-Max_Redo Mid to Tip 200A stretch.stl",...
%         "Scallop-Max_Redo Mid to Tip 200B stretch.stl"];
    
stls = ["WB33_Root to Mid 200A stretch.stl",...
        "WB33_Root to Mid 200B stretch.stl",...
        "WB33_Mid to Tip 200A stretch.stl",...
        "WB33_Mid to Tip 200B stretch.stl"];

stls_ref = ["clean_leading_edge_surf.stl"];
% ================================



% ==========Slice=================
[tt,slices] = STL_Main(stls,p0_total,n_total,display_fig);
[tt,slc_ref] = STL_Main(stls_ref,p0_total(1,:),n_total(1,:),display_fig);

% ================================



%% =========Pre Process Slice=========

% Smash slice to 2D
fprintf('\n=======\nSmashing 3d slice to 2d\n=======\n')

[sma_slice] = smash_slice(p0_total,n_total,slices,display_fig)
[sma_slice_ref] = smash_slice(p0_total(1,:),n_total(1,:),slc_ref,display_fig)

tic
% Take longest segment (ignore small ice cones)

for slc_idx = 1:size(sma_slice,2)
    
    slc_idx
    seg = segmentation(sma_slice{slc_idx})'
    sma_slice{slc_idx} = seg{1}; % pts2Ixy(seg{1});

end

seg = segmentation(sma_slice_ref{1})'
sma_slice_ref = seg{1}; % pts2Ixy(seg{1});
% sma_slice_ref{1}.Ix_s = sma_slice_ref{1}.Ix_s(2:end,:);
% sma_slice_ref{1}.Iy_s = sma_slice_ref{1}.Iy_s(2:end,:);
toc

figure(5);clf;
for i = 1:size(sma_slice,2)
    ff = sma_slice{i};
    if size(ff,1)>0
        plot(ff(:,1),ff(:,2),'linewidth',2);hold on;grid on;
    end
end
plot(sma_slice_ref(:,1),sma_slice_ref(:,2),'k','linewidth',2);hold on;grid on;

% For video purpose, find frame range
if save_video
    frame_range = [min(sma_slice_ref);max(sma_slice_ref)];
    for i = 1:length(sma_slice)
        frame_range(1,:) = min( [frame_range(1,:);min(sma_slice{i})] );
        frame_range(2,:) = max( [frame_range(2,:);max(sma_slice{i})] );
    end
    frame_range(1,:) = frame_range(1,:)-0.05;
    frame_range(2,:) = frame_range(2,:)+0.05;
end

%% =========Find Max Area Slice=========
% There are three formats of points:
% Ix, Iy 2D*2
% Pts [pt1x,pt1y;
%      pt2x,pt2y...] 2D
% Combined deg1:segments deg2:x,y deg3:start,end 3D

if save_video
    pastnow = num2str(now,12);
%     datestr(pastnow)
    mov_name = ['../results/maxarea_',pastnow,'.mp4'];
    delete(mov_name)
    myV = VideoWriter(mov_name,'MPEG-4');
    myV.Quality = 100;
    myV.FrameRate = 2;
    open(myV);
end


% Init
max_slice = sma_slice_ref;
for slc_idx = 1:length(sma_slice)
    temp_slice = sma_slice{slc_idx};

    % Split from intersections
    [slice_points,cutted_max,cutted_temp] = split_slices(max_slice,temp_slice);

    figure(1);clf
    for i = 1:min(length(cutted_max),length(cutted_temp))
        plot(cutted_max{i}(:,1),cutted_max{i}(:,2),'linewidth',2);hold on;grid on;
        plot(cutted_temp{i}(:,1),cutted_temp{i}(:,2),'linewidth',2);hold on;grid on;
    end

    % Find the curve that encloses larger area and stitch that in
    baseline_A = S_area(max_slice,false);
    vpa(baseline_A)
    found_flg = 0;

    % replace segments in max_slice to see if area is larger
    figure(2);clf
    for i_tmp = 1:length(cutted_temp)

        cur_seg = cutted_temp{i_tmp};
        cur_ht = [cur_seg(1,:);cur_seg(end,:)];
        cur_ht = sort_by_norm(cur_ht);
        for i_max = 1:length(cutted_max)

            max_ht = [cutted_max{i_max}(1,:);cutted_max{i_max}(end,:)];
            max_ht = sort_by_norm(max_ht);
            if isequal(cur_ht,max_ht)
                replace_idx = i_max;
                found_flg = 1;
            end

        end

        if found_flg
            figure(2);clf;
            tmp_segs = cutted_max;
            tmp_segs{replace_idx} = cur_seg;
            tmp_segs_pts = combine_segments(tmp_segs);
            compete_A = S_area(tmp_segs_pts,false);
            cutted_max_pts = combine_segments(cutted_max);
            
            if compete_A > baseline_A
                
                cutted_max{replace_idx} = cur_seg;
                baseline_A = compete_A;
%                 plot(tmp_segs_pts(:,1),tmp_segs_pts(:,2),'g');hold on;grid on;
                plot(cur_seg(:,1),cur_seg(:,2),'g');hold on;grid on;
            else
%                 plot(tmp_segs_pts(:,1),tmp_segs_pts(:,2),'r');hold on;grid on;
                plot(cur_seg(:,1),cur_seg(:,2),'r');hold on;grid on;
            end
            found_flg = 0;
            
            % Debug plot
            plot(cutted_max_pts(:,1),cutted_max_pts(:,2),'b');hold on;grid on;
            
            if save_video
                xlim(frame_range(:,1)');ylim(frame_range(:,2)');
                frame = getframe(gcf);
                writeVideo(myV,frame);
            end

        end
        max_slice = combine_segments(cutted_max);
    end
end

if save_video
    close(myV);
end

%%
figure(3);clf;
plot(max_slice(:,1),max_slice(:,2),'--b','linewidth',5);hold on;grid on;
for i = 1:size(sma_slice,2)
    ff = sma_slice{i};
    if size(ff,1)>0
        plot(ff(:,1),ff(:,2),'linewidth',2);hold on;grid on;
    end
end
plot(sma_slice_ref(:,1),sma_slice_ref(:,2),'k','linewidth',2);hold on;grid on;

%% test
% % Signed area
% pts = [-1,-1;
%        -1 ,0.5 ;
%        1 ,-1;
%        1,-2];
% 
% area = S_area(sma_slice_ref,display_fig)
% area = S_area(sma_slice{2},display_fig)

% % pt2line distance
% distance = pt2line(pt,v1,v2)

%% 
%  =======================================================
%  =================Helper Functions======================
%  =======================================================

function comb_seg = combine_segments(segs)

    sz_tmp = length(segs);
    comb_seg = segs;
    diff = sz_tmp;
    while diff>0
        comb_seg = iter_combine(comb_seg);
        diff = abs(length(comb_seg)-sz_tmp);
        sz_tmp = length(comb_seg);
    end
    if length(comb_seg)>1
        fprintf('potential data loss in combine_segments \n');
    end
    comb_seg = comb_seg{1};
end

function [sorted] = sort_by_norm(before)
    
    % sort by row norm, small comes first
    norm_before = vecnorm(before,2,2);
    [~,order] = sort(norm_before);
    sorted = before(order,:);

end

% function [slice_points,cutted_max,cutted_temp] = split_slices(max_slice,temp_slice)
% 
%     slice_points = [];
% 
% %     Find Intersection points
%     for idx_max_1 = 1:size(max_slice,1)
%         idx_max_2 = mod(idx_max_1,size(max_slice,1))+1;
%         for idx_temp_1 = 1:size(temp_slice,1)
%             idx_temp_2 = mod(idx_temp_1,size(temp_slice,1))+1;
%             
%             line1 = [max_slice(idx_max_1,:);max_slice(idx_max_2,:)];
%             line2 = [temp_slice(idx_temp_1,:);temp_slice(idx_temp_2,:)];
%             [Int,Int_log] = SegIntSeg(line1, line2);
%             if Int_log
%                 slice_points(end+1,:) = [Int,idx_max_1,idx_temp_1];
%             end
%             
%         end
%     end
%     
% %     Now split the slices into slice_pairs
%     %     Add in the intersection points
%     tic
%     for i = 1:size(slice_points,1)
%         temp_pt = slice_points(i,1:2);
%         max_slice  = insert_point(temp_pt,max_slice);
%         temp_slice = insert_point(temp_pt,temp_slice);
%     end
%     
%     
%     %     debug plot
%     figure(10);clf;
%     plot3(max_slice(:,1),max_slice(:,2),1:size(max_slice,1),'k','linewidth',2);hold on;grid on;
%     plot(temp_slice(:,1),temp_slice(:,2),'linewidth',2);hold on;grid on;
%     scatter(slice_points(:,1),slice_points(:,2))
%     
%     %     Split open slices at intersection points
%         %     Find cut points
%     cutpt_max = [];
%     cutpt_temp = [];
%     for i = 1:size(max_slice,1)
%         distan = vecnorm(slice_points(:,1:2)-max_slice(i,:),2,2);
%         if any(distan<1e-10)
%             cutpt_max(end+1) = i;
%         end
%     end
%     cutpt_max'
%     
%     for i = 1:size(temp_slice,1)
%         distan = vecnorm(slice_points(:,1:2)-temp_slice(i,:),2,2);
%         if any(distan<1e-10)
%             cutpt_temp(end+1) = i;
%         end
%     end
%     cutpt_temp'
%     
%         %     Cut at cut points
%     cutted_max = {};
%     for i = 1:length(cutpt_max)
%         i1 = cutpt_max(i);
%         i2 = cutpt_max(mod(i,length(cutpt_max))+1);
%         fprintf('i1: %d i2: %d\n',[i1,i2])
%         if i2<i1
%             cutted_max{end+1} = ([max_slice(i1:end,:);max_slice(1:i2,:)]);
%         else
%             cutted_max{end+1} = (max_slice(i1:i2,:));
%         end
%     end
%     
%     cutted_temp = {};
%     for i = 1:length(cutpt_temp)
%         i1 = cutpt_temp(i);
%         i2 = cutpt_temp(mod(i,length(cutpt_temp))+1);
%         if i2<i1
%             cutted_temp{end+1} = ([temp_slice(i1:end,:);temp_slice(1:i2,:)]);
%         else
%             cutted_temp{end+1} = (temp_slice(i1:i2,:));
%         end
%     end
%     
%     
%     toc
% end
% 
% function [slice] = insert_point(pt,slice)
% 
%     for i = 1:size(slice,1)
%         p1 = slice(i,:);
%         p2 = slice(mod(i,size(slice,1))+1,:);
%         if onSegment(p1,pt,p2)
%             slice = [slice(1:i,:);pt;slice(i+1:end,:)];
%             break
%         end
%     end
% 
% end
% 
% 
% % Given three colinear points p, q, r, the function checks if 
% % point q lies on line segment 'pr' 
% 
% function [onSeg] = onSegment(p,q,r) 
% 
%     onSeg = false;
%     if (q(1) <= max(p(1), r(1)) && q(1) >= min(p(1), r(1)) && q(2) <= max(p(2), r(2)) && q(2) >= min(p(2), r(2))) 
%         onSeg = true; 
%     end
% 
% end

function [dist] = pt2line(pt,v1,v2)

    a = v1 - v2;
    b = pt - v2;
    dist = norm(cross(a,b)) / norm(a);
    
end

% Note:Third col of pts will be ignored
% Use Concave signed area by cross product algorithm to calculte the signed
% area of any polygon.
function [area] = S_area(pts,disp_fig)
       
    if disp_fig
        figure(7);clf;
        plot([pts(:,1);pts(1,1)],[pts(:,2);pts(1,2)]);grid on;hold on;
        scatter(0,0);hold on;
        xlim([min(pts(:,1)),max(pts(:,1))]*1.5);ylim([min(pts(:,2)),max(pts(:,2))]*1.5);
    end
    
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

    area = abs(area);

end

function [pts] = Ixy2pts(slice)

    pts = [slice.Ix_s(:,1),slice.Iy_s(:,1)];

end

function [Ixy] = pts2Ixy(seg)
    
    Ixy = [];
    Ixy.Ix_s = [seg(:,1),[seg(end,1);seg(1:end-1,1)]];
    Ixy.Iy_s = [seg(:,2),[seg(end,2);seg(1:end-1,2)]];

end

function slice_comb = pts2combine(slice)

    slice_comb = [];
    slice_Ixy = pts2Ixy(slice);
    slice_comb(:,:,1) = slice_Ixy.Ix_s;
    slice_comb(:,:,2) = slice_Ixy.Iy_s;

end

function slice_pts = combine2pts(slice_comb)

    slice_Ixy = [];
    slice_Ixy.Ix_s = slice_comb(:,:,1);
    slice_Ixy.Iy_s = slice_comb(:,:,2);
    slice_pts = Ixy2pts(slice_Ixy);
    
end

function [slice_comb] = xy2combine(slice)

    slice_comb = [];
    slice_comb(:,:,1) = slice.Ix_s;
    slice_comb(:,:,2) = slice.Iy_s;

end

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

function [total] = proj2plane(p0,n,total)

    n = n/(norm(n));
    for i = 1:2
        dist = total(:,:,i)-p0;
        total(:,:,i) = total(:,:,i)-dist*n.'*n;
    end
    
end

function [Ix_s,Iy_s] = smash2D(total,e1,e2)
    
    Ix_s = zeros(size(total,1),size(total,3));
    Iy_s = zeros(size(total,1),size(total,3));
    
    for i = 1:2
        Ix_s(:,i) = total(:,:,i)*e1';
        Iy_s(:,i) = total(:,:,i)*e2';
    end

end