% =================================
% ==== Video processing script ====
% =================================

% Showcases generally how to save a video

figure(999);clf
mov_name = ['results/','testmov.mp4'];
myV = VideoWriter(mov_name,'MPEG-4');
myV.Quality = 100;
myV.FrameRate = 60;
open(myV)

for i = 1:1000
    scatter(rand(1),rand(1),'filled');hold on;grid on;
    frame = getframe(gcf);
    writeVideo(myV,frame);
end

close(myV)

%% frame rate modifier
% It will not generate new frame, it will only modify play speed

disp_vid = false;
mov_name = ['results/','maxarea_737597.476028.mp4'];
orig_V = VideoReader(mov_name)

modi_V = VideoWriter([mov_name(1:end-4),'_modFR.mp4'],'MPEG-4')
modi_V.Quality = 100;
modi_V.FrameRate = 3;
open(modi_V)

if disp_vid
    figure(999);clf;
    currAxes = axes;
end
while hasFrame(orig_V)
    vidFrame = readFrame(orig_V);
    if disp_vid
        image(vidFrame, 'Parent', currAxes);
        currAxes.Visible = 'off';
    end
    writeVideo(modi_V,vidFrame);
end

close(modi_V)

%% Frame interpolation
% It will try to interpolate between frames
disp_vid = false;
clc
fac = 4;
mov_name = ['results/','maxarea_737597.476028_modFR.mp4'];
% mov_name = ['results/','maxarea_737596.630126_modint.mp4'];
orig_V = VideoReader(mov_name)

modi_V = VideoWriter([mov_name(1:end-4),'_modFI.mp4'],'MPEG-4');
modi_V.Quality = 100;
modi_V.FrameRate = fac*orig_V.FrameRate
open(modi_V)

if disp_vid
    figure(999);clf;
    currAxes = axes;
end

prev_F = readFrame(orig_V);
i = 1;
while hasFrame(orig_V)
    curr_F = readFrame(orig_V);
    writeVideo(modi_V,prev_F);
    
    f1 = prev_F;
    f2 = curr_F;
    int_frames = interp_frame(f1,f2,fac);
    for j = 1:length(int_frames)
        writeVideo(modi_V,int_frames{j});
    end
    
    if disp_vid
        image(curr_F, 'Parent', currAxes);
        currAxes.Visible = 'off';
    end
    
    writeVideo(modi_V,curr_F);
    prev_F = curr_F;
    fprintf('Frame: %d \n',i)
    i = i+1;
end
writeVideo(modi_V,curr_F);
close(modi_V)




%% Helper Functions
function [int_frames] = interp_frame(f1,f2,fac)
    
    f1s = size(f1);
    f2s = size(f2);
    
    [X,Y,Z] = meshgrid(1:f1s(2),1:f1s(1),1:2);
    [Xq,Yq,Zq] = meshgrid(1:f1s(2),1:f1s(1),1:1/(fac+1):2);
    Vqs = {};
    for i = 1:3
        V = zeros(size(X));
        V(:,:,1) = f1(:,:,i); V(:,:,2) = f2(:,:,i);
        Vqs{end+1} = interp3(X,Y,Z,V,Xq,Yq,Zq);
    end

    int_frames = {};
    num_inter = size(Vqs{1},3)-2;
    for i = 2:num_inter+1
        temp_frame = zeros(f1s);
        temp_frame(:,:,1) = Vqs{1}(:,:,i);
        temp_frame(:,:,2) = Vqs{2}(:,:,i);
        temp_frame(:,:,3) = Vqs{3}(:,:,i);
        int_frames{end+1} = uint8(temp_frame);
    end

end