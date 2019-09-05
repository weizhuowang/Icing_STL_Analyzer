clc,clear;
% =========Slicer Set up==========
addpath(genpath('..'))
p0_total = []; n_total = [];
tot_idx = 10;
for i = tot_idx
% for i = 1
    p0_total(end+1,:) = [0,0,i];
    n_total(end+1, :) = [0,0,1];
end
display_fig = false;
save_video = false;
        
stls = ["acc_south.stl"];
%% ==========Slice=================
[tt,slc_trac,fv_trac] = STL_Main(stls,p0_total,n_total,display_fig);

%% smash into 2D

[sma_slice_trac] = smash_slice(p0_total,n_total,slc_trac,display_fig);
sma_slice_trac = segmentation(sma_slice_trac{1})';
sma_slice_trac_i = sma_slice_trac{1};
sma_slice_trac_o = sma_slice_trac{2};

scale_fac = 11.93/4.6;

sma_slice_trac_i(:,2) = 300-sma_slice_trac_i(:,2);
sma_slice_trac_o(:,2) = 300-sma_slice_trac_o(:,2);

sma_slice_trac_i = sma_slice_trac_i*scale_fac;
sma_slice_trac_o = sma_slice_trac_o*scale_fac;

%% plot slices
figure(2);clf
plot(sma_slice_trac_i(:,1),sma_slice_trac_i(:,2),'b');hold on;grid on;axis equal
plot(sma_slice_trac_o(:,1),sma_slice_trac_o(:,2),'r');hold on;
scatter(sma_slice_trac_i(:,1),sma_slice_trac_i(:,2),'b');hold on;
scatter(sma_slice_trac_o(:,1),sma_slice_trac_o(:,2),'r');hold on;

%% Definitions
max_acc = 0.7;
max_brk = 0.7;
max_cor = 1;
lean_spd = 40; % deg/s
%% 
%  =======================================================
%  =================Helper Functions======================
%  =======================================================

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

function LE = extract_LE(sma_slice_cln)

    % Leading edge point of CLN
    LE = [];
    for i = 1:length(sma_slice_cln)
        min_x = min(min(sma_slice_cln{i}.Ix_s));
        idx = find(sma_slice_cln{i}.Ix_s==min_x);
        LE(i,:) = [sma_slice_cln{i}.Ix_s(idx(1)),sma_slice_cln{i}.Iy_s(idx(1))];
    end

end

function LE_ice = useful_slice(sma_slice_ice,LE)

    % Extract useful part of ICE slice
    LE_ice = {};
    for i = 1:length(sma_slice_ice)
    % for i = 72
        fprintf('Processing %d/%d\n',[i,length(sma_slice_ice)])
        slice_temp = sma_slice_ice{i};
        pts = Ixy2pts(slice_temp);
        line_seg = [];
        for j = 1:size(pts,1)
            cur_pt = pts(j,:);
            if cur_pt(1) < LE(i,1)+1
                line_seg(end+1,:) = cur_pt;
            end
        end
        LE_ice{end+1} = line_seg;
    end

end

% return the horn height for the inquiry point relative to slice
% want: slice in Ixy format
function [hei,posid] = horn_height_pt(pt,slice_cln)

    height = [];
    for i = 1:size(slice_cln.Ix_s,1)
        pts = [slice_cln.Ix_s(i,:);slice_cln.Iy_s(i,:)].';
        if direct_above(pt,pts(1,:),pts(2,:))
            comp = (pt-pts(1,:))*(pts(2,:)-pts(1,:)).'/norm(pts(2,:)-pts(1,:));
            height(end+1,:) = [i,sqrt(norm(pt-pts(1,:))^2-comp^2)];
        end
    end
    if ~isempty(height)
        [~,minidx] = min(height(:,2));
        posid = height(minidx,1);
        hei = height(minidx,2);
    else
        posid = [];
        hei = [];
    end
end

