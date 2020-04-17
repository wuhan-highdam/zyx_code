mkdir('vectors_file')
data_dump = textread(['./surface_node_file/voro_node_',num2str(frame),'.sample.vol']);
data_dump2 = textread(['./cell_area_file/voro_node_',num2str(frame),'.sample.vol']);

%     data_dump3 = textread(['./cell_vectors_file/voro_node_',num2str(frame),'.sample.vol']);
voro_id = dlmread('./frame_all/voro_id_0.dat');
acc = 1;
par_id = zeros(length(data_dump),1);
for j = 1:length(voro_id);
    par_id(acc:(acc+voro_id(j)-1)) = j;
    acc = acc + voro_id(j);
end
acc = 1;
fid1 = fopen(['./cell_vectors_file/voro_node_',num2str(frame),'.sample.vol']);
%     eval(['data_area',num2str(i),'=cell(1,1)']);
%     eval(['data_vectors',num2str(i),'=cell(1,1)']);
data_area = cell(1,1);
data_vectors = cell(1,1);
for j = 1:length(voro_id);
    vectors = cell(1,1);
    js = 0;
    for ff = acc:(acc+voro_id(j)-1);
        js = js + 1;
        tline1 = fgetl(fid1);
        tline1(1) = '';
        tline1(end) = '';
        text_split = strsplit(tline1,') (');
        for k = 1:length(text_split);
            text_split{k} = str2num(text_split{k});
        end
        vectors{js} = text_split;
    end
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
                single_particle_neig = [single_particle_neig;ell_id,area(m,n),vectors{m}{n}];
            end
        end
    end
%         ell_id_sort = sort(unique(single_particle_neig(:,1)));
%         area_cluster = [];
%         for m = ell_id_sort';
%             pos_cluster = find(single_particle_neig(:,1)==m);
%             area_cluster = [area_cluster;m,sum(single_particle_neig(pos_cluster,2))];
%         end
%         dlmwrite(['./area_file4vectors/neighbor_',num2str(frame),'.dat'],[j,area_cluster(:,1)'],'-append','precision',7);
%         dlmwrite(['./area_file/area_',num2str(frame),'.dat'],[j,area_cluster(:,2)'],'-append','precision',7);

%         eval(['data_area',num2str(i),'{',num2str(j),'}=','single_particle_neig(:,2)']);
%         eval(['data_vectors',num2str(i),'{',num2str(j),'}=','single_particle_neig(:,3:end))']);
    data_area{j} = single_particle_neig(:,2);
    data_vectors{j} = single_particle_neig(:,3:end);
    acc = acc + voro_id(j);
end
save(['./vectors_file/vectors_',num2str(frame),'.mat'],'data_area','data_vectors');
