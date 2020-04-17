clear
tic
ffid = dlmread('inf_frame_50.dat');
ffid = ffid(:,[1,6,7,8]);
[s1,s2]=size(ffid);
A=ffid(:,1);
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
bin = 5e-4;  % Minimum radius divided by 15
x_min = floor(min(ffid(:,2))/bin)*bin;
x_max = ceil(max(ffid(:,2))/bin)*bin;
y_min = floor(min(ffid(:,3))/bin)*bin;
y_max = ceil(max(ffid(:,3))/bin)*bin;
z_min = floor(min(ffid(:,4))/bin)*bin;
z_max = ceil(max(ffid(:,4))/bin)*bin;
% domain = [-0.145,0.145, -0.145, 0.0.145, -0.145, 0.145];
% voxel_coor = [-0.0035:5e-4:0.2965];
% voxel_coor = [0.1:5e-4:0.2];
x_voxel_coor = [x_min:bin:x_max];
y_voxel_coor = [y_min:bin:y_max];
z_voxel_coor = [z_min:bin:z_max];

% x_voxel_coor = x_voxel_coor(ceil(0.65*length(x_voxel_coor)):end);
% y_voxel_coor = y_voxel_coor(ceil(0.65*length(y_voxel_coor)):end);
% z_voxel_coor = z_voxel_coor(ceil(0.65*length(z_voxel_coor)):end);
% voxel_coor3 = cell(length(voxel_coor),length(voxel_coor),length(voxel_coor));
num1 = length(x_voxel_coor);
num2 = length(y_voxel_coor);
num3 = length(z_voxel_coor);
num = [num1,num2,num3];
voxel_coor3 = zeros(num1*num2*num3,3);
for i = 1:num1;
    for j = 1:num2;
        for k = 1:num3;
            index = (k-1)*num1*num2 + (j-1)*num1 + i;
            voxel_coor3(index,:) = [x_voxel_coor(i),y_voxel_coor(j),z_voxel_coor(k)];
        end
    end
end

d=1;
oo=[];
voids = true(num1,num2,num3);
for f = 1:length(arrset);
    sp=ffid(d:(d+length(arrset{f})-1),2:4);
    ox=mean(sp(:,1));
    oy=mean(sp(:,2));
    oz=mean(sp(:,3));
    oo=[oo;ox,oy,oz];
    d=d+length(arrset{f});
    box_boundary1 = floor((min(sp)-[min(x_voxel_coor),min(y_voxel_coor),min(z_voxel_coor)])./bin);
    box_boundary2 = ceil((max(sp)-[min(x_voxel_coor),min(y_voxel_coor),min(z_voxel_coor)])./bin);
    for i = 1:3;
        if box_boundary1(i) < 1;
            box_boundary1(i) = 1;
        end
        if box_boundary2(i) > num(i);
            box_boundary2(i) = num(i);
        end
    end
    if box_boundary1(1) == box_boundary2(1) || box_boundary1(2) == box_boundary2(2) || box_boundary1(3) == box_boundary2(3);
        continue
    end
%     tri = convhulln(sp);
%     coor = [sp;oo(f,:)];
%     tri = [tri,length(coor)*ones(length(tri),1)];
%     tri = triangulation(tri,coor(:,1),coor(:,2),coor(:,3));
    [center, radii, evecs, v, chi2] = ellipsoid_fit(sp, '' );
    for i = box_boundary1(1):box_boundary2(1);
        for j = box_boundary1(2):box_boundary2(2);
            for k = box_boundary1(3):box_boundary2(3);
%                 in_particle = pointLocation(tri,[voxel_coor(i),voxel_coor(j),voxel_coor(k)]);
                in_particle = cal_in_particle(v,x_voxel_coor(i),y_voxel_coor(j),z_voxel_coor(k));
                if ~in_particle;
                    voids(i,j,k) = in_particle;
                end
            end
        end
    end
end

adress = find(voids == true);
adress3 = fix((adress-1)/(num1*num2))+1;
adress2 = fix(rem((adress-1),(num1*num2))/num1)+1;
adress1 = rem(rem((adress-1),num1*num2),num1)+1;
points = voxel_coor3(adress,:);
edt = bwdist(~voids); 
% Search cubic block of odd size 7*7*7
block_size = 7;
upper = block_size - 1;
shift = upper/2;
void_center_7 = false(num1,num2,num3);
background = zeros(num1+upper,num2+upper,num3+upper);
background(4:num1+3, 4:num2+3, 4:num3+3) = edt;

for f = 1:length(adress1);
    [i,j,k] = deal(adress1(f)+shift,adress2(f)+shift,adress3(f)+shift);
    search = background((i-shift):(i+shift), (j-shift):(j+shift), (k-shift):(k+shift));
    max_adress = find(search == max(max(max(search))));
    if max_adress(1) == 172;
        void_center_7(i-3,j-3,k-3) = true;
    end
end

center_adress = find(void_center_7 == true);
center_adress = center_adress(randperm(numel(center_adress)));  %´òÂÒcenterµÄÎ»ÖÃÅÅÐò
adress3_c = fix((center_adress-1)/(num1*num2))+1;
adress2_c = fix(rem((center_adress-1),(num1*num2))/num1)+1;
adress1_c = rem(rem((center_adress-1),num1*num2),num1)+1;
adress_c = [adress1_c,adress2_c,adress3_c];
voxel_type = zeros(length(adress1),1);
voxel_type_num = zeros(length(center_adress),1);

for f = 1:length(adress1);
    [i,j,k] = deal(adress1(f),adress2(f),adress3(f));
    box_voxel1 = [i,j,k] - 50;
    box_voxel2 = [i,j,k] + 50;
    center_adress_local_index = find(adress1_c> box_voxel1(1) & adress1_c <box_voxel2(1) & adress2_c > box_voxel1(2) & adress2_c < box_voxel2(2) & adress3_c > box_voxel1(3) & adress3_c < box_voxel2(3));
    center_adress_local = adress_c(center_adress_local_index,:);
    voxel_type(f) = center_adress_local_index(near_center([i,j,k],center_adress_local));
    voxel_type_num(voxel_type(f)) = voxel_type_num(voxel_type(f)) + 1;
end

myOutputFile1 = fopen('test-ellipsoid51.dump','w');
fprintf(myOutputFile1,'%s\r\n','ITEM: TIMESTEP');
fprintf(myOutputFile1,'%s\r\n','35800000');
fprintf(myOutputFile1,'%s\r\n','ITEM: NUMBER OF ATOMS');
fprintf(myOutputFile1,'%d\r\n',length(points));
fprintf(myOutputFile1,'%s\r\n','ITEM: BOX BOUNDS ff ff ff');
fprintf(myOutputFile1,'%4.4f %4.4f\r\n',x_min,x_max);
fprintf(myOutputFile1,'%4.4f %4.4f\r\n',y_min,y_max);
fprintf(myOutputFile1,'%4.4f %4.4f\r\n',z_min,z_max);
fprintf(myOutputFile1,'%s\r\n','ITEM: ATOMS id type x y z edt cluster_num');

id = [1:length(points)]';
edt1 = find(edt>0);
edt1 = double(edt(edt1));
voxel_type_num4id = zeros(length(points),1);
for i = 1 : length(points);
    voxel_type_num4id(i) = voxel_type_num(voxel_type(i));
end
data4print = [id,voxel_type,points(:,1),points(:,2),points(:,3),edt1,voxel_type_num4id]';
fprintf(myOutputFile1,'%d %d %6.4f %6.4f %6.4f %6.2f %d\r\n',data4print);
fclose all;
toc

% function  voxel_type = near_center(coor,center_pos);
%     dis_voxel_center = (coor - center_pos).^2;
%     dis_voxel_center = dis_voxel_center(:,1) + dis_voxel_center(:,2) +dis_voxel_center(:,3);
%     min_dis_adress = find(dis_voxel_center == min(dis_voxel_center));
%     voxel_type = min_dis_adress(1);
% end
% function result = is_closest(coor,coor_center);
%     dis_voxel_center = (coor - coor_center).^2;
%     dis_voxel_center = dis_voxel_center(:,1) + dis_voxel_center(:,2) +dis_voxel_center(:,3);
%     min_dis_adress = find(dis_voxel_center == min(dis_voxel_center));
%     result = min_dis_adress(1);
% end
