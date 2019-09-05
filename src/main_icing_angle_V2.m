clc,clear;
% =========Slicer Set up==========
addpath(genpath('..'))
p0_total = []; n_total = [];
tot_idx = 1:1:20;
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
    mov_name = ['../results/',vname(1:4),'_',pastnow,'t_',num2str(num_frame),'f.mp4'];
    delete(mov_name)
    myV = VideoWriter(mov_name,'MPEG-4');
    myV.Quality = 100;
    myV.FrameRate = length(tot_idx)/5;
    open(myV);
end

%% Extract LE and ices

LE = extract_LE(sma_slice_cln);

if display_fig
    figure(3);
    for i = 1:20
        subplot(ceil(length(plt_graph)/5),min(length(plt_graph),5),i)
        scatter(LE(i,1),LE(i,2),'filled');hold on
    end
    figure(4);clf;
end

LE_ice = useful_slice(sma_slice_ice,LE);

if display_fig
    for i = 1:length(sma_slice_ice)
        subplot(9,9,i)
        plot(line_seg(:,1),line_seg(:,2),'r');hold on; grid on;axis equal
        scatter(LE(i,1),LE(i,2),'filled');hold on
        title(['slice ',num2str(i)])
    end
end
%% Horn angle 

% Preperation
horn_angle_tot = [];
horn_height_tot = [];
% figure(7);clf;
    
for slc_idx = 1:length(tot_idx)
    tic
    
    %==========Horn Height===========
    temp_cln = xy2combine(sma_slice_cln{slc_idx});
    temp_ice = [];
    temp_ice(:,:,1) = LE_ice{slc_idx};
    temp_ice(:,:,2) = [LE_ice{slc_idx}(2:end,:);LE_ice{slc_idx}(1,:)];
    temp_ice = permute(temp_ice,[1,3,2]);

    % % Analyze angle, Horn angle definition
    figure(5);clf;
%     subplot(7,11,slc_idx)
    plot(LE_ice{slc_idx}(:,1),LE_ice{slc_idx}(:,2),'r');hold on; grid on;
    plot(sma_slice_cln{slc_idx}.Ix_s',sma_slice_cln{slc_idx}.Iy_s','b');hold on;
    axis equal

    % Use numerical method to find max distance point
    max_dist_pt = [];
    max_dist = 0;
    for i = 1:size(temp_cln,1)
        n = [squeeze(temp_cln(i,1,:))'-squeeze(temp_cln(i,2,:))',0];
        p0 = [squeeze(temp_cln(i,1,:))',0];
        if norm(n,2) > 0
            [IntP] = PlaneIntSlice(n,p0,temp_ice);
            if length(IntP) > 0
%                 scatter(IntP(:,1),IntP(:,2),'k');hold on
                dist_vec = IntP-p0;
                [max_val,max_idx] = max(vecnorm(dist_vec,2,2));
                if max_val > max_dist
                    max_dist = max_val;
                    max_dist_pt = IntP(max_idx,:);
                end
            end
        end
    end
    
    %==========Horn Height===========
    horn_height_tot(end+1) = max_dist;
    
    %==========Horn Angle===========
    coord_diff = max_dist_pt(1:2) - LE(slc_idx,:);
    horn_angle = rad2deg(atan2(coord_diff(2),coord_diff(1)));
    horn_angle = 180-mod(horn_angle+360,360);
    horn_angle_tot(end+1) = horn_angle;
    
    scatter(max_dist_pt(:,1),max_dist_pt(:,2),'c','filled');hold on
    scatter(LE(slc_idx,1),LE(slc_idx,2),'c','filled');hold on
    plot([max_dist_pt(:,1),LE(slc_idx,1)],[max_dist_pt(:,2),LE(slc_idx,2)],'k');hold on
    xlim([-0.6,0.6]);
    ylim([-0.8,0.6]);
    drawnow
    
    if save_video
        frame = getframe(gcf);
        writeVideo(myV,frame);
    end
    
    fprintf('Horn angle(deg): %6.2f | Time: %4.2f s | Progress: %2d/%d \n',[horn_angle,toc,slc_idx,length(tot_idx)]);
end

if save_video
    close(myV);
end

%% Fit line
syms x
figure(6);clf;
set(gcf,'Position', [60 580 1800 420])
scatter(tot_idx,horn_angle_tot,'r');hold on;grid on;
% plot(tot_idx,horn_angle_tot,'r');hold on;grid on;
ha_fit = fit(tot_idx',horn_angle_tot','poly2')
plot(ha_fit,tot_idx',horn_angle_tot','b--')
xlabel('Span');ylabel('Horn angle (deg)');title([vname(1:4),' Horn Angle with line fit'])
saveas(gcf,['..\results\',vname(1:4),'_ha_wfit.png'])

figure(7);clf;
set(gcf,'Position', [60 80 1800 420])
scatter(tot_idx,horn_height_tot,'r');hold on;grid on;
% plot(tot_idx,horn_height_tot,'b');hold on;grid on;
hh_fit = fit(tot_idx',horn_height_tot','poly2')
plot(hh_fit,tot_idx',horn_height_tot','b--')
xlabel('Span');ylabel('Horn height');title([vname(1:4),' Horn Height with line fit'])
saveas(gcf,['..\results\',vname(1:4),'_hh_wfit.png'])
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

