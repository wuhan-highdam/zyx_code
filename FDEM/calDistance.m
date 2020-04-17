%% 计算矩阵中点与点之间的距离
function [ dis ] = calDistance( x )
    [m,n]=size(x);
    dis = zeros(m,m);
    x=x';
    for i = 1:m;
        for j = m:-1:i;
            dis(i, j)=pdist(x(:, [i j])');
            dis(j, i) = dis(i, j);
        end
     end
end