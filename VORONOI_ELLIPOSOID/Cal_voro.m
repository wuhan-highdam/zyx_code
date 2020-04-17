clear
mkdir('neighbors_file4area');
mkdir('area_file');
mkdir('vol_file')
dump_path = './surface_node_file';
dirOutput=dir(fullfile(dump_path,'*.vol'));
fileNames={dirOutput.name}';
file_prefix = 'voro_node_';
file_suffix = '.sample.vol';
prefix_len = length(file_prefix);
suffix_len = length(file_suffix);
frame_all = [];
for i = 1:length(fileNames)
    if fileNames{i}(end-suffix_len+1:end) == file_suffix;
        frame_all = [frame_all,str2num(fileNames{i}(prefix_len+1:end-suffix_len))];
    end
end
frame_all = sort(frame_all);
clust = parcluster;
parpool(clust,16);
parfor i = 1:length(frame_all);
    frame = frame_all(i);
    data_dump = textread(['./surface_node_file/voro_node_',num2str(frame),'.sample.vol']);
    data_dump2 = textread(['./cell_area_file/voro_node_',num2str(frame),'.sample.vol']);
    data_dump3 = textread(['./cell_vol_file/voro_node_',num2str(frame),'.sample.vol']);
    voro_id = dlmread(['./frame_all/voro_id_',num2str(frame),'.dat']);
    cell_vol =[];
    acc = 1;
    par_id = zeros(length(data_dump),1);
    for j = 1:length(voro_id);
        par_id(acc:(acc+voro_id(j)-1)) = j;
        acc = acc + voro_id(j);
    end
    acc = 1;
    for j = 1:length(voro_id);
        vol = data_dump3(acc:(acc+voro_id(j)-1),2);
        cell_vol = sum(vol);
        dlmwrite(['./vol_file/cell_vol_',num2str(frame),'.dat'],[j,cell_vol],'-append');
        neighbors = data_dump(acc:(acc+voro_id(j)-1),2:end);
        area = data_dump2(acc:(acc+voro_id(j)-1),2:end);
        single_particle_neig = [];
        [s1,s2] = size(neighbors);
        for m = 1:s1;
            for n = 1:s2;
                if neighbors(m,n)==0;
                    continue
                elseif neighbors(m,n)<0;
                    ell_id = neighbors(m,n);
                else
                    ell_id = par_id(neighbors(m,n));
                end
                if ell_id==j;
                    continue
                else
                    single_particle_neig = [single_particle_neig;ell_id,area(m,n)];
                end
            end
        end
        ell_id_sort = sort(unique(single_particle_neig(:,1)));
        area_cluster = [];
        for m = ell_id_sort';
            pos_cluster = find(single_particle_neig(:,1)==m);
            area_cluster = [area_cluster;m,sum(single_particle_neig(pos_cluster,2))];
        end
        dlmwrite(['./neighbors_file4area/neighbor_',num2str(frame),'.dat'],[j,area_cluster(:,1)'],'-append','precision',7);
        dlmwrite(['./area_file/area_',num2str(frame),'.dat'],[j,area_cluster(:,2)'],'-append','precision',7);
        acc = acc + voro_id(j);
    end
    fclose all;
end
delete(gcp('nocreate'))
