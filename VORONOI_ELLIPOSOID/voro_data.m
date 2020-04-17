clear
mkdir('surface_node_file');
mkdir('cell_area_file');
mkdir('cell_vol_file');
path(path,'/home/yxzou/workfs/Matlab_code')
fileFolder = [pwd,'/frame_all'];
dirOutput=dir(fullfile(fileFolder,'*.dat'));
fileNames={dirOutput.name}';
file_prefix = 'inf_frame_';
file_suffix = '.dat';
prefix_len = length(file_prefix);
suffix_len = length(file_suffix);
dump_frame = [];
for f = 1:1:length(fileNames);
    if fileNames{f}(1:prefix_len) == file_prefix;
        dump_frame = [dump_frame,str2num(fileNames{f}(prefix_len+1:end-suffix_len))];
    end
end
dump_frame = sort(dump_frame);
mark_end = dump_frame(end);
myOutputFile = fopen('words4voro.sh','w');
for f = 1:length(dump_frame);
    if dump_frame(f) == mark_end;
        fprintf(myOutputFile,'%s\n',['bash ',file_prefix,num2str(dump_frame(f)),'.sh']);
    else
        fprintf(myOutputFile,'%s\n',['bash ',file_prefix,num2str(dump_frame(f)),'.sh &']);
    end
end
fclose(myOutputFile);
clust = parcluster;
parpool(clust,16);
parfor f = 1:length(dump_frame);
    temp=[fileFolder,'/',file_prefix,num2str(dump_frame(f)),file_suffix];
    inf = dlmread(temp);
    [s1,s2]=size(inf);
    A=inf(:,1);
    c1 = 1;
    arrset = cell(0,0);
    while(c1<numel(A))
        c2 = 0;
        while (c1+c2+1<=numel(A)&&A(c1)+c2+1==A(c1+c2+1))
            c2 = c2+1;
        end
        if(c2>=1)
            arrset= [arrset;(A(c1:1:c1+c2))];
        end
        c1 = c1 + c2 +1;
    end
    d=1;
    surface_date=cell(0,0);
    for i=1:length(arrset);
        node_particle=inf(d:(d+length(arrset{i})-1),6:8);
        oo = [mean(node_particle(:,1)),mean(node_particle(:,2)),mean(node_particle(:,3))];
        K = convhulln(node_particle);
        single_surface_date =[];
        for j=1:length(node_particle);
            node_new = erosion(j,K,node_particle,oo);
            single_surface_date = [single_surface_date;node_new];
        end
        d=d+length(arrset{i});
        surface_date{i}=single_surface_date;
    end
    ff = fopen([fileFolder,'/voro_id_',num2str(dump_frame(f)),'.dat'],'w');
    for i = 1:length(surface_date);
        fprintf(ff,'%d\r\n',length(surface_date{i}));
    end
    ffid = fopen([fileFolder,'/voro_node_',num2str(dump_frame(f)),'.sample'],'w');
    num = 1;
    for i = 1:length(surface_date);
        for j = 1:length(surface_date{i});
            fprintf(ffid,'%d %12.6f %12.6f %12.6f\r\n',num,surface_date{i}(j,1),surface_date{i}(j,2),surface_date{i}(j,3));
            num = num + 1;
        end
    end
    boundary_file = fopen([fileFolder,'/boundary_node_',num2str(dump_frame(f)),'.dat'],'w');
    fprintf(boundary_file,'%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\r\n',min(inf(:,6)),max(inf(:,6)),min(inf(:,7)),max(inf(:,7)),min(inf(:,8)),max(inf(:,8))); 
  	myOutputFile2 = fopen(['surface_node_file/', file_prefix,num2str(dump_frame(f)),'.sh'],'w');
    myOutputFile3 = fopen(['cell_area_file/', file_prefix,num2str(dump_frame(f)),'.sh'],'w');
    myOutputFile4 = fopen(['cell_vol_file/', file_prefix,num2str(dump_frame(f)),'.sh'],'w');
    OutputFile_name = ['voro_node_',num2str(dump_frame(f)),'.sample'];
    fprintf(myOutputFile2,'%s %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %s\n','voro++ -c "%i %n" -o',min(inf(:,6)),max(inf(:,6)),min(inf(:,7)),max(inf(:,7)),min(inf(:,8)),max(inf(:,8)),OutputFile_name);
    fprintf(myOutputFile3,'%s %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %s\n','voro++ -c "%i %f" -o',min(inf(:,6)),max(inf(:,6)),min(inf(:,7)),max(inf(:,7)),min(inf(:,8)),max(inf(:,8)),OutputFile_name);
    fprintf(myOutputFile4,'%s %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %s\n','voro++ -c "%i %v" -o',min(inf(:,6)),max(inf(:,6)),min(inf(:,7)),max(inf(:,7)),min(inf(:,8)),max(inf(:,8)),OutputFile_name);
    fclose(myOutputFile2);
    fclose(myOutputFile3);
    fclose(myOutputFile4);
    fclose all;
end
