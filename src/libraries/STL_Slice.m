function [Intersection_x,Intersection_y,Intersection_z,segment_center,seg_len] = STL_Slice(n,p0,facets)
%STL_Slice Return intersection coordinates from STL facets

    % Generate intersection point clouds
    
    Intersection_x = zeros(size(facets,3),2);
    Intersection_y = zeros(size(facets,3),2);
    Intersection_z = zeros(size(facets,3),2);
    id = 1;
    for i = 1:size(facets,3)
        [Intersect,IntP] = PlaneIntTrig(n,p0,facets(:,:,i));
        if Intersect
            Intersection_x(id,:) = IntP(:,1)';
            Intersection_y(id,:) = IntP(:,2)';
            Intersection_z(id,:) = IntP(:,3)';
            id = id+1;
        end
    end
    
    Intersection_x(id:end,:) = [];
    Intersection_y(id:end,:) = [];
    Intersection_z(id:end,:) = [];
    
    Segments = zeros(size(Intersection_x,1),size(Intersection_x,2),3);
    Segments(:,:,1) = Intersection_x;
    Segments(:,:,2) = Intersection_y;
    Segments(:,:,3) = Intersection_z;
    Segments = permute(Segments,[1,3,2]);
    
    segment_center = (Segments(:,:,1)+Segments(:,:,2))/2;
    seg_len = sqrt(sum((Segments(:,:,1)-Segments(:,:,2)).^2, 2));
    
    %**********Line segment identification UNFINISHED************
    % current_line = [Segments(1,:,1);Segments(1,:,2)];
    % all_lines = {};
    % 
    % % ****Loop over all segments****
    % 
    % for i = 2:size(Segments,1)
    % 
    %     curr_pt = Segments(i,:,1);
    %     if norm(curr_pt-current_line(end,:))<1e-5
    %         current_line(end+1,:) = Segments(i,:,2);
    %     else
    %         all_lines{end+1} = current_line;
    %         current_line = [Segments(i,:,1);Segments(i,:,2)];
    %     end
    %     prev_pt = curr_pt;
    %     
    % end
    
end

