clc,clear;
% =========Slicer Set up==========
addpath(genpath('..'))
p0_total = []; n_total = [];
tot_idx = 1:1:9;
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

%% smash into 2D

[sma_slice_cln,sma_slice_ice] = smash_align(p0_total,n_total,slc_cln,slc_ice,display_fig);

%% Plot slices

if display_fig
    figure(3);clf
    plt_num = 1
    % plt_graph = [18]
    plt_graph = 1:min(length(tot_idx),25)
    for ff = plt_graph
        ff
        subplot(ceil(length(plt_graph)/5),min(length(plt_graph),5),plt_num)
        plt_num = plt_num + 1;
        plot(sma_slice_cln{ff}.Ix_s',sma_slice_cln{ff}.Iy_s','b');hold on;grid on;
        plot(sma_slice_ice{ff}.Ix_s',sma_slice_ice{ff}.Iy_s','r');hold on;axis equal
        title(['slice ',num2str(ff)])
    end
    drawnow
end

if save_video
    num_frame = length(tot_idx);
    pastnow = num2str(now,12);
    vname = char(stls_ice(1));
    mov_name = ['../results/','unroll_',pastnow,'t_',num2str(num_frame),'f.mp4'];
    delete(mov_name)
    myV = VideoWriter(mov_name,'MPEG-4');
    myV.Quality = 100;
    myV.FrameRate = length(tot_idx)/5;
    open(myV);
end

%% Extract LE and ices
% % Probably wont need this part, not sure yet.

% LE = extract_LE(sma_slice_cln);
% 
% if display_fig
%     figure(3);
%     for i = 1:20
%         subplot(ceil(length(plt_graph)/5),min(length(plt_graph),5),i)
%         scatter(LE(i,1),LE(i,2),'filled');hold on
%     end
%     figure(4);clf;
% end
% 
% LE_ice = useful_slice(sma_slice_ice,LE);
% 
% if display_fig
%     for i = 1:length(sma_slice_ice)
%         subplot(9,9,i)
%         plot(LE_ice{i}(:,1),LE_ice{i}(:,2),'r');hold on; grid on;axis equal
%         scatter(LE(i,1),LE(i,2),'filled');hold on
%         title(['slice ',num2str(i)])
%     end
% end

%% Horn height (New Lossless algo) 
% % BUG:point will show no found if two segments do not point directly to the
% % point.
% 
% % Find non-obtuse triangles (including right angle)
% test_slice = Ixy2pts(sma_slice_ice{1});
% ff = [];
% plt_idx = [];
% plotdata = [];cur_x = 0;
% for i = 1:size(test_slice,1)
%     [height,posid] = horn_height_pt(test_slice(i,:),sma_slice_cln{1});
%     if length(height) == 0
%         ff(end+1) = i;
%         fprintf('not found\n')
%     else
%         plt_idx(end+1) = posid;
%         plotdata(end+1,:) = [posid,height];
%     end
% end

%% Horn height (New Lossless algo V2) 
% Main idea: use dist2seg to find distance to every point on clean airfoil.
% Min distance would be the Min distance

LE = extract_LE(sma_slice_cln);
LE_ice = useful_slice(sma_slice_ice,LE);

for slcidx = 3
%     test_slice_i = Ixy2pts(sma_slice_ice{slcidx});
    test_slice_i = LE_ice{slcidx};
    test_slice_c = xy2combine(sma_slice_cln{slcidx});

    % Generate x axis tick scale 
    cur_length = 0;
    for i = 1:size(test_slice_c,1)
        cumu_length(i,:) = cur_length;
        cur_seg = squeeze(test_slice_c(i,:,:));
        cur_length = cur_length + norm(cur_seg(2,:)-cur_seg(1,:),2);
    end

    % Generate distance to slice (TODO:with shift)
    unroll_dist = [];
    rel_x = [];
    rel_x_idx = [];
    for j = 1:size(test_slice_i,1)
        temp_dist  = [];
        temp_shift = [];
        if mod(j,100) == 0
            fprintf('unrolling %d %d/%d \n',[slcidx,j,size(test_slice_i,1)])
        end
        for i = 1:size(test_slice_c,1)
            [temp_dist(i,:),temp_shift(i,:)] = dist2seg(test_slice_i(j,:),squeeze(test_slice_c(i,:,:)));
        end
        [unroll_dist(j,:),rel_x_idx(j,:)] = min(temp_dist);
        final_shift = temp_shift(rel_x_idx(j,:));
        rel_x(j,:) = cumu_length(rel_x_idx(j,:))+final_shift;
    end

    % plot
    figure(112);clf
    plot(sma_slice_cln{slcidx}.Ix_s',sma_slice_cln{slcidx}.Iy_s','b');hold on;grid on;
    plot(test_slice_i(:,1),test_slice_i(:,2),'r');hold on;grid on;
    axis equal
    xlim([min(test_slice_i(:,1))-0.2,max(test_slice_i(:,1))+0.2])
    ylim([min(test_slice_i(:,2))-0.2,max(test_slice_i(:,2))+0.2])
    figure(113);clf
    plot(rel_x,unroll_dist,'r');grid on; hold on;
    plot([-10,20],[0,0],'--b','linewidth',1);hold on;
    xlim([min(rel_x)-0.1,max(rel_x)+0.1])
    ylim([min(unroll_dist)-0.05,max(unroll_dist)+0.05])
    axis equal
    drawnow
    % (150:1300)

    if save_video
        frame = getframe(gcf);
        writeVideo(myV,frame);
    end

end

if save_video
    close(myV);
end


%%
figure(111);clf;
% plot([plt_idx,plt_idx+max(plt_idx)].',[plotdata(:,2);plotdata(:,2)]);grid on;hold on;
scatter([plt_idx].',[plotdata(:,2)]);grid on;hold on;
scatter([ff,ff+max(plt_idx)].',ones(1,length(ff)*2));hold on;
figure(112);clf;
plot(test_slice(:,1),test_slice(:,2));hold on;grid on;

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

