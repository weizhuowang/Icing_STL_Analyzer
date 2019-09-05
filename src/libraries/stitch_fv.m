function [fv_stitched] = stitch_fv(stls)
%STITCH_FV stitch multiple fv object into one large fv
%   [fv_stitched] = stitch_fv(stls)

    fv_stitched = [];
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