clc,clear;
% =========Slicer Set up==========
addpath(genpath('..'))
p0_total = [-0.2,0,0];
n_total  = [1,0,0];
display_fig = true;
stls = ["Scallop-Max_Redo Root to Mid 200A stretch.stl",...
        "Scallop-Max_Redo Root to Mid 200B stretch.stl",...
        "Scallop-Max_Redo Mid to Tip 200A stretch.stl",...
        "Scallop-Max_Redo Mid to Tip 200B stretch.stl"];
% ================================



% ==========Slice=================
[total_time,slices,fv] = STL_Main(stls,p0_total,n_total,display_fig);
% ================================




% =======Smash slice to 2D========
fprintf('\n=======\nSmashing 3d slice to 2d\n=======\n')

[smashed_slices] = smash_slice(p0_total,n_total,slices,display_fig)



%% =========Periodic Analysis=========
smashed_slices
tic
segment = segmentation(smashed_slices{1})';
size(segment)
toc

%% obtain center and segment length
seg_size = [];
filtered_seg = {};
figure(5);clf;
for i = 1:size(segment,1)
    ff = segment{i};
    if size(ff,1)>80
        plot(ff(:,1),ff(:,2),'linewidth',2);hold on;grid on;
        seg_size(end+1) = size(segment{i},1);
        filtered_seg{end+1} = ff;
    end
end

seg_center = [];
seg_len = [];
for i = 1:size(filtered_seg,2)
    ff = filtered_seg{i};
    seg_center(i,:) =  mean(ff);
    seg_len(i,:) = sum(vecnorm([ff(end,:);ff(1:end-1,:)]-ff,2,2));
end
scatter(seg_center(:,1),seg_center(:,2),'ro','filled');grid on;hold on;

total = sortrows([seg_center,seg_len],1);
seg_center = total(:,1:2);
seg_len = total(:,3);

figure(2);clf;
scatter(1:size(seg_size,2),seg_size);grid on;

figure(3);clf;
scatter3(seg_center(:,1),seg_center(:,2),seg_len);grid on;
axis equal


%% center distances
dist_center = vecnorm([seg_center(end,:);seg_center(1:end-1,:)]-seg_center,2,2);
figure(4);clf
histogram(dist_center(2:end),120);grid on;

figure(6);clf;
plot(seg_center(2:end,1),dist_center(2:end));grid on;
xlabel('y');ylabel('distance');
%%

figure(1);clf;
scatter3(p0_total(1)*ones(size(seg_center,1),1),-seg_center(:,1),seg_center(:,2),'ro','filled');grid on;hold on;
axis equal
patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);grid on; xlabel('x');ylabel('y');zlabel('z');
view([-90 0]);
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

for i = 1:size(segment,1)
    ff = segment{i};
    if size(ff,1)>80
        plot3(p0_total(1)*ones(size(ff,1),1),-ff(:,1),ff(:,2),'linewidth',2);hold on;grid on;
    end
end

%% test
[segment1] = iter_combine({segment{:},segment{:}});


%% helpers

function display_facets(fv)

    patch(fv,'FaceColor',       [0.8 0.8 1.0], ...
             'EdgeColor',       'none',        ...
             'FaceLighting',    'gouraud',     ...
             'AmbientStrength', 0.15);grid on;

    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');

    % Fix the axes scaling, and set a nice view angle
    axis('image');
    view([-135 35]);

end