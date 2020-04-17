clear
path(path,'/workfs/yxzou/Matlab_code')
mkdir('local_area');
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
    cell_area = dlmread(['./area_file/area_',num2str(frame),'.dat']);
    cell_area = sum(cell_area(:,2:end),2);
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
%         [K,v] = convhulln(node);
        [center, radii, evecs, v, chi2] = ellipsoid_fit(node, '' );
        a=radii(1);
        b=radii(2);
        c=radii(3);
        fun=@(A,B) a*b.*cos(A).*...
            sqrt(sin(A).^2+(c/a.*cos(A).*cos(B)).^2+(c/b.*cos(A).*sin(B)).^2);
        area=8*quad2d(fun,0,pi/2,0,pi/2);
        porosity(j) = (cell_area(j) - area)/cell_area(j);
        dlmwrite(['./local_area/local_area_',num2str(frame),'.dat'],[j,porosity(j)'],'-append','precision',7);
        d=d+length(arrset{j});
    end
    fclose all;
end