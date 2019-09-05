% FIXED
% HTPOOLS structure:
%            3D                      3D
% [#1 start_x , #1 end_x] [#1 start_y , #1 end_x]
% [#2 start_x , #2 end_x] [#2 start_y , #2 end_x]
% [          ...        ] [          ...        ]

% VP is validity pool (flag denotes whether an element is deleted or not)

% squeeze's behavior is weird if there is only one row in the 3D array

% Potential bug: add a point into the segment
function [segment] = iter_combine(prev_segment)

    if size(prev_segment,1) > size(prev_segment,2)
        prev_segment = prev_segment.';
    end
    
    vp = 1:size(prev_segment,2);
    ht_pools = zeros(size(prev_segment,2),2,2);
    
    % pool all head and tail points
    for i = 1:size(prev_segment,2)
        cur_line = prev_segment{i};
        ht_pools(i,1,:) = cur_line(1,:);
        ht_pools(i,2,:) = cur_line(end,:);
    end
    % Add dummy data to avoid weird behavior of squeeze
    ht_pools(end+1,1,:) = [6666,6666];
    
    segment = {};
    
    tol = 1e-10;

    while size(ht_pools,1)>1
    
        % Extract data for this loop
        cur_seg = prev_segment{vp(1)};
%         size(cur_seg)
%         prev_segment{1} = [];
%         prev_segment = prev_segment(~cellfun('isempty',prev_segment));
        vp(1) = [];
        ht_pools(1,:,:) = [];
        
        
        if size(ht_pools,1) > 1
            % duplicate location at head (at other head, at other tail)
%            cur_seg(1,:)
%            cur_seg(end,:)
%            squeeze(ht_pools(:,1,:))
%            squeeze(ht_pools(:,2,:))
           for dup = 1:5
               switch dup
                   case 1
                       dup_loc = count_occur(cur_seg(1,:),squeeze(ht_pools(:,1,:)),tol,-1);
                       if dup_loc ~= -1
                           break
                       end
                   case 2
                       dup_loc = count_occur(cur_seg(1,:),squeeze(ht_pools(:,2,:)),tol,-1);
                       if dup_loc ~= -1
                           break
                       end
                   case 3
                       dup_loc = count_occur(cur_seg(end,:),squeeze(ht_pools(:,1,:)),tol,-1);
                       if dup_loc ~= -1
                           break
                       end
                   case 4
                       dup_loc = count_occur(cur_seg(end,:),squeeze(ht_pools(:,2,:)),tol,-1);
                       if dup_loc ~= -1
                           break
                       end
                   otherwise
               end
           end
        else
            dup = 5;
        end
        
%         BUG IS RIGHT HERE, SQUEEZE HTPOOLS should be SQUEEZE from Line

        % duplicate location at head (at other head)
        if dup == 1
            combined_seg = prev_segment{vp(dup_loc)};
            segment{end+1} = [flipud(cur_seg); combined_seg(2:end,:)];
            ht_pools(dup_loc,:,:) = [];
%             prev_segment{dup_loc} = [];
%             prev_segment = prev_segment(~cellfun('isempty',prev_segment));
            vp(dup_loc) = [];
        % duplicate location at head (at other tail)
        elseif dup == 2
            combined_seg = prev_segment{vp(dup_loc)};
            segment{end+1} = [combined_seg(1:end-1,:);cur_seg];
            ht_pools(dup_loc,:,:) = [];
%             prev_segment{dup_loc} = [];
%             prev_segment = prev_segment(~cellfun('isempty',prev_segment));
            vp(dup_loc) = [];
        % duplicate location at tail (at other head)
        elseif dup == 3
            combined_seg = prev_segment{vp(dup_loc)};
            segment{end+1} = [cur_seg;combined_seg(2:end,:)];
            ht_pools(dup_loc,:,:) = [];
%             prev_segment{dup_loc} = [];
%             prev_segment = prev_segment(~cellfun('isempty',prev_segment));
            vp(dup_loc) = [];
        % duplicate location at tail (at other tail)
        elseif dup == 4
            combined_seg = prev_segment{vp(dup_loc)};
            segment{end+1} = [cur_seg;flipud(combined_seg(1:end-1,:))];
            ht_pools(dup_loc,:,:) = [];
%             prev_segment{dup_loc} = [];
%             prev_segment = prev_segment(~cellfun('isempty',prev_segment));
            vp(dup_loc) = [];
        else
            segment{end+1} = cur_seg;
        end
    
    end

end

% Return locations of duplicate row in the array
function [occ] = count_occur(vec,total,tol,curidx)

    occ = find(sum(abs(vec-total),2)<tol);
    if isempty(occ)
        occ = -1;
    end
    
    occ = occ(1);
    
end



function [segment] = iter_combine_old(prev_segment)

    ht_pools = zeros(size(prev_segment,2),2,2);
    for i = 1:size(prev_segment,2)
        cur_line = prev_segment{i};
        ht_pools(i,1,:) = cur_line(1,:);
        ht_pools(i,2,:) = cur_line(end,:);
    end
    segment = {};
    
    tol = 1e-10;

    while size(ht_pools,1)>0
    
        % Extract data for this loop
        cur_seg = prev_segment{1};
%         prev_segment{1} = [];
%         prev_segment = prev_segment(~cellfun('isempty',prev_segment));
        prev_segment(1) = [];
        ht_pools(1,:,:) = [];
        
        
        if size(ht_pools,1) > 0
            % duplicate location at head (at other head, at other tail)
            
           for dup = 1:5
               switch dup
                   case 1
                       dup_loc = count_occur(cur_seg(1,:),squeeze(ht_pools(:,1,:)),tol,-1);
                       if dup_loc ~= -1
                           break
                       end
                   case 2
                       dup_loc = count_occur(cur_seg(1,:),squeeze(ht_pools(:,2,:)),tol,-1);
                       if dup_loc ~= -1
                           break
                       end
                   case 3
                       dup_loc = count_occur(cur_seg(2,:),squeeze(ht_pools(:,1,:)),tol,-1);
                       if dup_loc ~= -1
                           break
                       end
                   case 4
                       dup_loc = count_occur(cur_seg(2,:),squeeze(ht_pools(:,2,:)),tol,-1);
                       if dup_loc ~= -1
                           break
                       end
                   otherwise
               end
           end
        else
            dup = 5;
        end
        
%         BUG IS RIGHT HERE, SQUEEZE HTPOOLS should be SQUEEZE from Line

        % duplicate location at head (at other head)
        if dup == 1
            combined_seg = prev_segment{dup_loc};
            segment{end+1} = [flipud(cur_seg); combined_seg(2:end,:)];
            ht_pools(dup_loc,:,:) = [];
%             prev_segment{dup_loc} = [];
%             prev_segment = prev_segment(~cellfun('isempty',prev_segment));
            prev_segment(dup_loc) = [];
        % duplicate location at head (at other tail)
        elseif dup == 2
            combined_seg = prev_segment{dup_loc};
            segment{end+1} = [combined_seg(1:end-1,:);cur_seg];
            ht_pools(dup_loc,:,:) = [];
%             prev_segment{dup_loc} = [];
%             prev_segment = prev_segment(~cellfun('isempty',prev_segment));
            prev_segment(dup_loc) = [];
        % duplicate location at tail (at other head)
        elseif dup == 3
            combined_seg = prev_segment{dup_loc};
            segment{end+1} = [cur_seg;combined_seg(2:end,:)];
            ht_pools(dup_loc,:,:) = [];
%             prev_segment{dup_loc} = [];
%             prev_segment = prev_segment(~cellfun('isempty',prev_segment));
            prev_segment(dup_loc) = [];
        % duplicate location at tail (at other tail)
        elseif dup == 4
            combined_seg = prev_segment{dup_loc};
            segment{end+1} = [cur_seg;flipud(combined_seg(1:end-1,:))];
            ht_pools(dup_loc,:,:) = [];
%             prev_segment{dup_loc} = [];
%             prev_segment = prev_segment(~cellfun('isempty',prev_segment));
            prev_segment(dup_loc) = [];
        else
            segment{end+1} = cur_seg;
        end
    
    end

end
