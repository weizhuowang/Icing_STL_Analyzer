clc,clear;
fprintf('=======\nRead STL\n=======\n');tic;

% Read in stl file and put the points in general format accepted by this
% program.

fv_stitched = [];

stls = ["Scallop-Max_Redo Root to Mid 200A stretch.stl",...
        "Scallop-Max_Redo Root to Mid 200B stretch.stl",...
        "Scallop-Max_Redo Mid to Tip 200A stretch.stl",...
        "Scallop-Max_Redo Mid to Tip 200B stretch.stl"];

idx = 1;
for stlname = stls
    tic
    fprintf('Loading : %s\n',stlname)
    tria = stlread(stlname);

    fv.vertices = tria.Points;
    fv.faces = tria.ConnectivityList;

    clear tria
    
    if idx == 1
        fv_stitched = fv;
    else
        shift = size(fv_stitched.vertices,1);
        fv_stitched.vertices = [fv_stitched.vertices;fv.vertices];
        fv_stitched.faces    = [fv_stitched.faces;fv.faces+shift];
    end

    clear fv
    idx = idx+1;
    
    toc
end

% facets = uncompress_fv(fv_stitched);
    
% figure(1);clf
% display_downsample(facets,50000)
% display_facets(fv)


%%
function display_downsample(facets,sample_size)

    disp_fv = [];
%     sample_size = 50000;
    sample_idx = randi(size(facets,3),min(size(facets,3),sample_size),1);
    sample_facets = facets(:,:,sample_idx);

    disp_fv.vertices = reshape(permute(sample_facets,[2,1,3]),3,[])';
    disp_fv.faces = reshape([1:size(disp_fv.vertices,1)],3,[])';

    scatter3(disp_fv.vertices(:,1),disp_fv.vertices(:,2),disp_fv.vertices(:,3),2,'filled','MarkerFaceColor',[0 .75 .75])
    xlabel('X');ylabel('Y');zlabel('Z');
    axis equal
    view([-3 1.5]);
end

function facets = uncompress_fv(fv)

    % % Uncompress the facets from the fv struct to facets, an n pages facet book

    tic
    % facets = newSTLfacets(fv);
    facets = zeros(3,3,size(fv.faces,1));
    temp = fv.faces';
    temp = temp(:);
    facets = fv.vertices(temp,:);
    facets = permute(reshape(facets',[3,3,size(facets,1)/3]),[2,1,3]);
    toc

end


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