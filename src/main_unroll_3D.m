clc,clear;
% =========Slicer Set up==========
addpath(genpath('..'))
p0_total = []; n_total = [];
tot_idx = 1:1;
for i = tot_idx
% for i = 1
    p0_total(end+1,:) = [0,i,0];
    n_total(end+1, :) = [0,1,0];
end
display_fig = false;
save_video = false;

% stls_ice = ["Scallop-Max_Redo Root to Mid 200A stretch.stl",...
%             "Scallop-Max_Redo Root to Mid 200B stretch.stl",...
%             "Scallop-Max_Redo Mid to Tip 200A stretch.stl",...
%             "Scallop-Max_Redo Mid to Tip 200B stretch.stl"];
        
stls_ice = ["WB33_Root to Mid 200A stretch.stl",...
            "WB33_Root to Mid 200B stretch.stl",...
            "WB33_Mid to Tip 200A stretch.stl",...
            "WB33_Mid to Tip 200B stretch.stl"];

stls_cln = ["clean_leading_edge_surf.stl"];
vname = char(stls_ice(1));

%% ==========Slice=================
[tt,slc_ice,fv_ice] = STL_Main(stls_ice,p0_total,n_total,display_fig);
[tt,slc_cln,fv_cln] = STL_Main(stls_cln,p0_total,n_total,display_fig);

%% Generate unrolled STL
% (Lossless algo)
% Main idea: Find unrolled coordinate for each point, reconstruct the STL
% TODO: For all the points to the left of clean slice
% Notes: Ice STL : 8M points
%        CLN STL : 27K points

% Rotate STL (no need)
% fv_cln_backup = fv_cln;
% fv_cln.vertices = fv_cln.vertices*rotz(37.2);

facets_cln = newSTLfacets(fv_cln);
figure(1);clf;
display_downsample(facets_cln,50000,[0 .75 .75]);hold on;

facets_ice = newSTLfacets(fv_ice);
display_downsample(facets_ice,50000,[.5 .25 .25]);hold on;

function unr_pt = unroll_Coord(ice_pt, fv_cln)

    unr_pt = 0;

end

% Uncompress the facets from the fv struct to facets, an n pages facet
% book. Same as line 30 to 36
function facets =  newSTLfacets(fv)

    facets = zeros(3,3,size(fv.faces,1));
    temp = fv.faces';
    temp = temp(:);
    facets = fv.vertices(temp,:);
    facets = permute(reshape(facets',[3,3,size(facets,1)/3]),[2,1,3]);

end

function display_downsample(facets,sample_size,MKR_FC)

    disp_fv = [];
%     sample_size = 50000;
    sample_idx = randi(size(facets,3),min(size(facets,3),sample_size),1);
    sample_facets = facets(:,:,sample_idx);

    disp_fv.vertices = reshape(permute(sample_facets,[2,1,3]),3,[])';
    disp_fv.faces = reshape([1:size(disp_fv.vertices,1)],3,[])';

    scatter3(disp_fv.vertices(:,1),disp_fv.vertices(:,2),disp_fv.vertices(:,3),2,'filled','MarkerFaceColor',MKR_FC)
    xlabel('X');ylabel('Y');zlabel('Z');
    axis equal
    view([-3 1.5]);
end
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

