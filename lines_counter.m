count_dir = 'C:\Users\Wang Weizhuo\Documents\Summer 2019\BSW research\V3.00\src';
cd(count_dir)
files = dir([count_dir,'\**\*.m']);
n = [];
for i = 1:size(files,1)
    fid = fopen([files(i).folder,'\',files(i).name]);
    [nnn,ns,nall] = linecount(fid);
    n(end+1,:) = [nnn,ns,nall];
    fclose(fid);
end
n
sum(n)

function [nnn,ns,nall] = linecount(fid)
    nnn = 0; % no space + no comment
    ns = 0;  % no space
    nall = 0;
    
    tline = fgetl(fid);
    while ischar(tline)
      tline = fgetl(fid);
      if sum(~isspace(tline))
          if tline(1) ~= '%'
              nnn = nnn+1;
          end
          ns = ns+1;
      end
      nall = nall + 1;
    end
end