function [smashed_slices] = smash_slice(p0_total,n_total,slices,display_fig)
%SMASH_SLICE Smash a 3D slice to 2D
%   [smashed_slices] = smash_slice(p0_total,n_total,slices,display_fig)


    % Define which plane to project to
    p0 = p0_total(1,:);
    n = n_total(1,:);

    % Define smashed 2D vector
    if isequal(n,[0,0,1])
        e1 = cross(n,[1,0,0])/norm(cross(n,[1,0,0]));
    else 
        e1 = cross(n,[0,0,1])/norm(cross(n,[0,0,1]));
    end
    e2 = cross(e1,n)/norm(cross(e1,n));

    % extract and smash
    if display_fig
        clf;
    end
    smashed_slices = {};
    for i = 1:size(slices,2)
        fprintf('Smashing slice %d\n',i)
        % Extract
        slice = slices{i};
        Ix = slice.Ix;
        Iy = slice.Iy;
        Iz = slice.Iz;
        centers = slice.centers;

        % Preview
        if display_fig
            subplot(1,2,1)
            plot3(Ix',Iy',Iz','b');grid on;hold on;
            xlabel('X');ylabel('Y');zlabel('Z');
            axis equal
        end

        % Smash
            % Prepare coordinates
        total = zeros([size(Ix),3]);
        total(:,:,1) = Ix;
        total(:,:,2) = Iy;
        total(:,:,3) = Iz;
        total = permute(total,[1,3,2]);

        total = proj2plane(p0,n,total);
        [Ix_s,Iy_s] = smash2D(total,e1,e2);
            % Preview smash
        if display_fig
            subplot(1,2,2)
            plot(Ix_s',Iy_s','r');grid on;hold on;
            axis equal
        end

        % Save in format
        temp_struct = [];
        temp_struct.Ix_s = Ix_s;
        temp_struct.Iy_s = Iy_s;
        smashed_slices{i} = temp_struct;
    end

    fprintf('Complete\n\n')

end

function [total] = proj2plane(p0,n,total)

    n = n/(norm(n));
    for i = 1:2
        dist = total(:,:,i)-p0;
        total(:,:,i) = total(:,:,i)-dist*n.'*n;
    end
    
end

function [Ix_s,Iy_s] = smash2D(total,e1,e2)
    
    Ix_s = zeros(size(total,1),size(total,3));
    Iy_s = zeros(size(total,1),size(total,3));
    
    for i = 1:2
        Ix_s(:,i) = total(:,:,i)*e1';
        Iy_s(:,i) = total(:,:,i)*e2';
    end

end