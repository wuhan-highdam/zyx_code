function  voxel_type = near_center(coor,center_pos);
    dis_voxel_center = (coor - center_pos).^2;
    dis_voxel_center = dis_voxel_center(:,1) + dis_voxel_center(:,2) +dis_voxel_center(:,3);
    min_dis_adress = find(dis_voxel_center == min(dis_voxel_center));
    voxel_type = min_dis_adress(1);
end