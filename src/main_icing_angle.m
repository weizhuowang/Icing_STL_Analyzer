clc,clear;
% =========Slicer Set up==========
addpath(genpath('..'))
p0_total = []; n_total = [];
tot_idx = 1:1:74;
for i = tot_idx
% for i = 1
    p0_total(end+1,:) = [0,i,0];
    n_total(end+1, :) = [0,1,0];
end
display_fig = false;
save_video = true;
% stls_ice = ["Scallop-Max_Redo Root to Mid 200A stretch.stl",...
%             "Scallop-Max_Redo Root to Mid 200B stretch.stl",...
%             "Scallop-Max_Redo Mid to Tip 200A stretch.stl",...
%             "Scallop-Max_Redo Mid to Tip 200B stretch.stl"];
        
stls_ice = ["WB33_Root to Mid 200A stretch.stl",...
            "WB33_Root to Mid 200B stretch.stl",...
            "WB33_Mid to Tip 200A stretch.stl",...
            "WB33_Mid to Tip 200B stretch.stl"];

stls_cln = ["clean_leading_edge_surf.stl"];


% ================================
%% ==========Slice=================
[tt,slc_ice,fv_ice] = STL_Main(stls_ice,p0_total,n_total,display_fig);
[tt,slc_cln,fv_cln] = STL_Main(stls_cln,p0_total,n_total,display_fig);


%% smash into 2D

tic

[sma_slice_cln] = smash_slice(p0_total,n_total,slc_cln,display_fig);
[sma_slice_ice] = smash_slice(p0_total,n_total,slc_ice,display_fig);

toc

% Take longest segment (ignore small ice cones)
tic

tic
for slc_idx = 1:size(sma_slice_ice,2)
    
    slc_idx
    seg = segmentation(sma_slice_ice{slc_idx})'
%     temp{slc_idx} = seg;
    seg = seg{1};
    sma_slice_ice{slc_idx}.Ix_s = [seg(:,1),[seg(end,1);seg(1:end-1,1)]];
    sma_slice_ice{slc_idx}.Iy_s = [seg(:,2),[seg(end,2);seg(1:end-1,2)]];

end

% Align all pairs of ice and clean slice
for i = 1:size(sma_slice_cln,2)
    i
    sma_slice_cln{i} = shift_slice_2d(sma_slice_cln{i},sma_slice_ice{i});
end

toc

%% Plot slices
tic

if display_fig
    figure(3);clf
    plt_num = 1
    % plt_graph = [18]
    plt_graph = 1:25
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
toc
%% Extract LE and ices

% Leading edge point of CLN
LE = [];
for i = 1:length(sma_slice_cln)
    min_x = min(min(sma_slice_cln{i}.Ix_s));
    idx = find(sma_slice_cln{i}.Ix_s==min_x);
    LE(i,:) = [sma_slice_cln{i}.Ix_s(idx(1)),sma_slice_cln{i}.Iy_s(idx(1))];
end

if display_fig
    figure(3);
    for i = 1:25
        subplot(ceil(length(plt_graph)/5),min(length(plt_graph),5),i)
        scatter(LE(i,1),LE(i,2),'filled');hold on
    end
end

% Extract useful part of ICE
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
    if display_fig
        figure(4);clf;
        subplot(9,9,i)
        plot(line_seg(:,1),line_seg(:,2),'r');hold on; grid on;axis equal
        scatter(LE(i,1),LE(i,2),'filled');hold on
        title(['slice ',num2str(i)])
    end
end

%% Horn angle 
% % Analyse angle (in deg) OLD
% LE_ang = {};
% figure(5);clf
% for i = 1:length(LE_ice)
%     cur_line = LE_ice{i};
%     tmp = [];
%     for j = 1:size(cur_line,1)
%         cur_pt = cur_line(j,:);
%         coord_diff = cur_pt - LE(i,:);
%         tmp(end+1) = rad2deg(atan2(coord_diff(2),coord_diff(1)));
%     end
%     tmp = mod(tmp+360,360);
%     LE_ang{end+1} = tmp;
%     subplot(5,5,i)
%     plot(1:length(tmp),tmp,'b');hold on; grid on;
% end

% Preperation
horn_angle_tot = [];
horn_height_tot = [];
% figure(7);clf;
if save_video
    num_frame = length(tot_idx);
    pastnow = num2str(now,12);
    vname = char(stls_ice(1));
%     datestr(pastnow)
    mov_name = ['../results/',vname(1:4),'_',pastnow,'t_',num2str(num_frame),'f.mp4'];
    delete(mov_name)
    myV = VideoWriter(mov_name,'MPEG-4');
    myV.Quality = 100;
    myV.FrameRate = length(tot_idx)/5;
    open(myV);
end
    
    
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
%     plt_circle(LE(slc_idx,1),LE(slc_idx,2),norm(coord_diff,2))
    
%     xlim([min(temp_ice(:,1,1)-0.2),max(temp_cln(:,1,1))+0.2]);
%     ylim([min(temp_cln(:,1,2)-0.2),max(temp_cln(:,1,2))+0.2]);

%     xlim([min(temp_ice(:,1,1)-0.3),max(temp_ice(:,1,1))+0.3]);
%     ylim([min(temp_ice(:,1,2)-0.3),max(temp_ice(:,1,2))+0.3]);
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

%% Show 3d plot
% figure(1);clf
% subplot(2,1,1);
% display_facets(fv_cln)
% title('clean')
% subplot(2,1,2);
% display_facets(fv_ice)
% title('icing')
% drawnow


%%
%  =======================================================
%  =================Helper Functions======================
%  =======================================================

function [pts] = Ixy2pts(slice)

    pts = [slice.Ix_s(:,1),slice.Iy_s(:,1)];
%     pts = zeros(2*size(slice.Ix_s,1),2);
%     for i = 1:size(slice.Ix_s,1)
%         pts(2*i,:) = [slice.Ix_s(i,1),slice.Iy_s(i,1)];
%         pts(2*i+1,:) = [slice.Ix_s(i,2),slice.Iy_s(i,2)];
%     end

end

% shift the clean slice to match the icing one
function [slice_cln] = shift_slice_orig(slice_cln,slice_ice)

    slice_cln_comb = xyz2combine(slice_cln);
    slice_ice_comb = xyz2combine(slice_ice);

    corr = slice_diff(slice_cln_comb,slice_ice_comb,[0,0,1],[0,0,0.25])

    % Shift the slice to fit iceing slice
    slice_cln.Ix = slice_cln.Ix+corr(1);
    slice_cln.Iy = slice_cln.Iy+corr(2);
    slice_cln.Iz = slice_cln.Iz+corr(3);
    
end

function slice_comb = xy2combine(slice)

    slice_comb = [];
    slice_comb(:,:,1) = slice.Ix_s;
    slice_comb(:,:,2) = slice.Iy_s;

end

function slice_comb = xyz2combine(slice)

    slice_comb = [];
    slice_comb(:,:,1) = slice.Ix;
    slice_comb(:,:,2) = slice.Iy;
    slice_comb(:,:,3) = slice.Iz;

end

function display_facets(fv)

    patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);grid on;
    xlabel('x');ylabel('y');zlabel('z');
    
    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');

    % Fix the axes scaling, and set a nice view angle
    axis('image');
    view([-90 10]);

end

% Determine the y=pt(2) line intersection with the segment, Int_y = pt(2)
function Int_x = seg_int_x(p1,p2,pt)

    p1x = p1(1);
    p2x = p2(1);
    p1y = p1(2);
    p2y = p2(2);
    
    Int_x = p2x - (p1x-p2x)/(p1y-p2y)*(p2y-pt(2));

end

function plt_circle(x,y,r)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit, yunit,'k--');
end