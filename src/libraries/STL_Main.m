function [total_time,slices,fv] = STL_Main(stls,p0_total,n_total,display_figure)
% Last modified Feb 16. Functionalized main script
% single_slice = [Ix,Iy,Iz,centers,seg_len]

    % ===========
    % Control Panel

%     display_figure = false;
    start_t        = tic;
    % ===========

    fprintf('=======\nRead STL\n=======\n');tic;

    % Read in stl file and put the points in general format accepted by this
    % program.

    fv = stitch_fv(stls);

    toc

    
    fprintf('\n=======\nDisplay STL\n=======\n')

    % % Uncompress the facets from the fv struct to facets, an n pages facet book

    tic
    % facets = newSTLfacets(fv);
    facets = zeros(3,3,size(fv.faces,1));
    temp = fv.faces';
    temp = temp(:);
    facets = fv.vertices(temp,:);
    facets = permute(reshape(facets',[3,3,size(facets,1)/3]),[2,1,3]);
    toc

 

    fprintf('\n=======\nSetting Plane Definition\n=======\n')

    % Print range Info about the STL

    fprintf('Axis max [X|Y|Z]:');
    fprintf('| %2.2f |',max(fv.vertices))
    fprintf('\nAxis min [X|Y|Z]:');
    fprintf('| %2.2f |',min(fv.vertices))
    fprintf('\n');

    % Plane definition

    if isempty(p0_total)
        p0_total = (max(fv.vertices)+min(fv.vertices))/2
    end
    if isempty(n_total)
        n_total = [1,2,3];
    end 

    
    for slice_num = 1:size(n_total,1)

        if display_figure
            figure(slice_num);clf;
        end
    
        fprintf('\n=======\nPicking Relevant triangles\n=======\n')
        
        p0 = p0_total(slice_num,:);
        n  = n_total(slice_num,:);
        fprintf('Now working on: point(%2.2f,%2.2f,%2.2f)\n',p0)
        fprintf('Now working on: n    (%2.2f,%2.2f,%2.2f)\n',n)
        
        % % Display the STL in point cloud format, if desired
        
        if display_figure
            tic
            subplot(1,3,1);
            % Draw sampled points
            display_downsample(facets,50000);hold on;
            % Solve for plane definitions
            % Draw slicing plane
            model_range = max(fv.vertices)-min(fv.vertices);
            sz = abs(model_range(1))/2;
            if n(3) == 0
                [z1,z2,z3,z4] = deal(1,1,-1,-1);
                fill3(p0(1)+sz*[1 -1 -1 1], p0(2)+sz*[-1 1 1 -1], p0(3)+sz*[z1 z2 z4 z3],'r','FaceAlpha',0.7)
            else
                syms z
                z1 = double(solve(n*[1,1,z].'==0));
                z2 = double(solve(n*[1,-1,z].'==0));
                z3 = double(solve(n*[-1,1,z].'==0));
                z4 = double(solve(n*[-1,-1,z].'==0));
                fill3(p0(1)+sz*[1 1 -1 -1], p0(2)+sz*[1 -1 -1 1], p0(3)+sz*[z1 z2 z4 z3],'r','FaceAlpha',0.7)
            end
            % Draw norm line
            plot3(p0(1)+[0,n(1)],p0(2)+[0,n(2)],p0(3)+[0,n(3)],'k','LineWidth',3)
            drawnow;
            toc
        end
        


        % Pick out the relevant facets and put into useful_facets book

        tic
        % Find projected distance of each point in facet to p0 on direction n
        proj_dist2D = reshape(permute(facets-p0,[2,1,3]),3,[])'; % vec distance: flattened projection to speedup calculation
        proj_dist2D = reshape(proj_dist2D*n',3,[]);              % projected distance on n, has 3 entry, each is a vortex of triangle
        max_p = max(proj_dist2D);
        min_p = min(proj_dist2D);
        clear proj_dist2D
        useful_idx = (sign(max_p.*min_p) == -1);
        clear max_p min_p
        
        useful_facets = facets(:,:,useful_idx==1);
        toc

        % Display the relevant facets, if desired

        if display_figure
            subplot(1,3,2);
            disp_fv = [];
            disp_fv.vertices = reshape(permute(useful_facets,[2,1,3]),3,[])';
            disp_fv.faces = reshape([1:size(disp_fv.vertices,1)],3,[])';
            display_downsample(useful_facets,50000)
            drawnow;
        end


        fprintf('\n=======\nGenerating Slice\n=======\n')

        % Generate the slice using the plane definition above

        tic
        [Ix,Iy,Iz,centers,seg_len] = STL_Slice(n,p0,useful_facets);
        single_slice = [];
        single_slice.p0 = p0;
        single_slice.n  = n;
        single_slice.Ix = Ix;
        single_slice.Iy = Iy;
        single_slice.Iz = Iz;
        single_slice.centers = centers;
        single_slice.seg_len = seg_len;
        toc

        % Show the slices if desired

        if display_figure
            subplot(1,3,3);
            plot3(Ix',Iy',Iz','b');grid on;
            xlabel('X');ylabel('Y');zlabel('Z');
            axis equal
            drawnow;
        end
        
        slices{slice_num} = single_slice;
        
    end
    
    
    total_time = toc(start_t)

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

% Uncompress the facets from the fv struct to facets, an n pages facet
% book. Same as line 30 to 36
function facets =  newSTLfacets(fv)

    facets = zeros(3,3,size(fv.faces,1));
    temp = fv.faces';
    temp = temp(:);
    facets = fv.vertices(temp,:);
    facets = permute(reshape(facets',[3,3,size(facets,1)/3]),[2,1,3]);

end