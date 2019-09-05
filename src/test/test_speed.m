% ts = segment;
% ts = {ts{:},ts{:}};
% ts{817} = []
% 
% 
% tic
% for i = 1:10000
%     ff = ts(~cellfun('isempty',ts));
% end
% toc
% 
% 
% tic 
% for i = 1:10000
%     idx = false(size(ts,2),1);
%     for i = 1:size(ts,2)
%         sz = size(ts{i},1);
%         idx(i) = logical(sign(sz));
%     end
%     ff = ts(idx);
% end
% toc
% 
% 
% tic
% 
% for i = 1:10000
%     ff = ts;
%     ff(817) = [];
% end
% 
% toc

clc,clear;
% =========Slicer Set up==========
addpath(genpath('data'))
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