processforShearlets

Files = dir('*.*');

for k = 4:7
   FileNames = Files(k).name;
   [~,~,ext] = fileparts(FileNames);
   c = strcmp(ext,'.png');
   if c == 1
      [MapsX, Mapping, Xrec] = rundirwavLE(FileNames);
   else
   end
end