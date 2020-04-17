clear
path(path,'/home/yxzou/workfs/Matlab_code')
fileFolder = [pwd,'/frame_all'];
dirOutput=dir(fullfile(fileFolder,'*.sample'));
fileNames={dirOutput.name}';
file_prefix = 'voro_node_';
file_suffix = '.sample';
prefix_len = length(file_prefix);
suffix_len = length(file_suffix);
dump_frame = [];
for f = 1:1:length(fileNames);
    if fileNames{f}(1:prefix_len) == file_prefix;
        dump_frame = [dump_frame,str2num(fileNames{f}(prefix_len+1:end-suffix_len))];
    end
end
dump_frame = sort(dump_frame);
for f = 1:length(dump_frame);
    name1 = ['./cell_area_file/inf_frame_',num2str(dump_frame(f)),'.sh'];
    name2 = ['./cell_vectors_file/inf_frame_',num2str(dump_frame(f)),'.sh'];
    fid1 = fopen(name1,'r+');
    fid2 = fopen(name2,'w');
    text = fgetl(fid1);
    text = strrep(text,'%i %f','%l');
    fprintf(fid2,'%s\n',text);
    fclose all;
end
    