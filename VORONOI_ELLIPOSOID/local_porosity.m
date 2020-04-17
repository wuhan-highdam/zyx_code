clear
% path(path,'/home/zhouwei/workfs/zouyuxiong/Matlab_code')
mkdir('local_porosity');
dump_path = './vol_file';
dirOutput=dir(fullfile(dump_path,'*.dat'));
fileNames={dirOutput.name}';
file_prefix = 'cell_vol_';
file_suffix = '.dat';
prefix_len = length(file_prefix);
suffix_len = length(file_suffix);
frame_all = [];
for i = 1:length(fileNames)
    if fileNames{i}(end-suffix_len+1:end) == file_suffix;
        frame_all = [frame_all,str2num(fileNames{i}(prefix_len+1:end-suffix_len))];
    end
end
frame_all = sort(frame_all);
parfor i = 1:length(frame_all);
    frame = frame_all(i);
    inf = textread(['./frame_all/inf_frame_',num2str(frame),'.dat']);
    cell_vol = dlmread(['./vol_file/cell_vol_',num2str(frame),'.dat']);
    cell_vol = cell_vol(:,2);
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
    zx=[];
    d=1;
    jc=[];
    ox=0;
    oy=0;
    oz=0;
    oo=[];
    jiao=[];
    d=1;
    porosity = zeros(length(arrset),1);
    for j=1:length(arrset);
        sp=inf(d:(d+length(arrset{j})-1),1:12);
        node=sp(:,6:8);
        [K,v] = convhulln(node);
        porosity(j) = (cell_vol(j) - v)/cell_vol(j);
        dlmwrite(['./local_porosity/local_porosity_',num2str(frame),'.dat'],[j,porosity(j)'],'-append','precision',7);
        d=d+length(arrset{j});
    end
    fclose all;
end