function [slice_points,cutted_max,cutted_temp] = split_slices(max_slice,temp_slice)

    slice_points = [];

%     Find Intersection points
    for idx_max_1 = 1:size(max_slice,1)
        idx_max_2 = mod(idx_max_1,size(max_slice,1))+1;
        for idx_temp_1 = 1:size(temp_slice,1)
            idx_temp_2 = mod(idx_temp_1,size(temp_slice,1))+1;
            
            line1 = [max_slice(idx_max_1,:);max_slice(idx_max_2,:)];
            line2 = [temp_slice(idx_temp_1,:);temp_slice(idx_temp_2,:)];
            [Int,Int_log] = SegIntSeg(line1, line2);
            if Int_log
%                 slice_points(end+1,:) = [Int,idx_max_1,idx_temp_1];
                slice_points = [slice_points;[Int,idx_max_1,idx_temp_1]];
            end
            
        end
    end
    
%     Now split the slices into slice_pairs
    %     Add in the intersection points
%     tic
    for i = 1:size(slice_points,1)
        temp_pt = slice_points(i,1:2);
        max_slice  = insert_point(temp_pt,max_slice);
        temp_slice = insert_point(temp_pt,temp_slice);
    end
    
    
    %     debug plot
    figure(10);clf;
    plot3(max_slice(:,1),max_slice(:,2),1:size(max_slice,1),'k','linewidth',2);hold on;grid on;
    plot(temp_slice(:,1),temp_slice(:,2),'linewidth',2);hold on;grid on;
%     scatter(slice_points(:,1),slice_points(:,2))   %%
    
    %     Split open slices at intersection points
        %     Find cut points
    cutpt_max = [];
    cutpt_temp = [];
    for i = 1:size(max_slice,1)
        distan = [];
        for j = 1:size(slice_points,1)
            distan = [distan;norm(slice_points(j,1:2)-max_slice(i,:),2)];
        end
        if any(distan<1e-10)
%             cutpt_max(end+1) = i;
            cutpt_max = [cutpt_max,i];
        end
    end
    cutpt_max'
    
    for i = 1:size(temp_slice,1)
        distan = [];
        for j = 1:size(slice_points,1)
            distan = [distan;norm(slice_points(j,1:2)-temp_slice(i,:),2)];
        end
%         distan = vecnorm(slice_points(:,1:2)-temp_slice(i,:),2,2);
        if any(distan<1e-10)
%             cutpt_temp(end+1) = i;
            cutpt_temp = [cutpt_temp,i];
        end
    end
    cutpt_temp'
    
        %     Cut at cut points
    cutted_max = {};
    for i = 1:length(cutpt_max)
        i1 = cutpt_max(i);
        i2 = cutpt_max(mod(i,length(cutpt_max))+1);
%         fprintf('i1: %d i2: %d\n',[i1,i2])   %%
        if i2<i1
            cutted_max{end+1} = ([max_slice(i1:end,:);max_slice(1:i2,:)]);
        else
            cutted_max{end+1} = (max_slice(i1:i2,:));
        end
    end
    
    cutted_temp = {};
    for i = 1:length(cutpt_temp)
        i1 = cutpt_temp(i);
        i2 = cutpt_temp(mod(i,length(cutpt_temp))+1);
        if i2<i1
            cutted_temp{end+1} = ([temp_slice(i1:end,:);temp_slice(1:i2,:)]);
        else
            cutted_temp{end+1} = (temp_slice(i1:i2,:));
        end
    end
    
    
%     toc
end

function [slice] = insert_point(pt,slice)

    for i = 1:size(slice,1)
        p1 = slice(i,:);
        p2 = slice(mod(i,size(slice,1))+1,:);
        if onSegment(p1,pt,p2)
            slice = [slice(1:i,:);pt;slice(i+1:end,:)];
            break
        end
    end

end


% Given three colinear points p, q, r, the function checks if 
% point q lies on line segment 'pr' 

function [onSeg] = onSegment(p,q,r) 

    onSeg = false;
    if (q(1) <= max(p(1), r(1)) && q(1) >= min(p(1), r(1)) && q(2) <= max(p(2), r(2)) && q(2) >= min(p(2), r(2))) 
        onSeg = true; 
    end

end