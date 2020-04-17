clear
path(path,'F:\桌面_9.6\后处理代码\FDEM')
inf = dlmread('inf_frame_0.dat');
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

myOutputFile1 = fopen('FDEM_ellipsoid_OVITO.dump','w');
fprintf(myOutputFile1,'%s\r\n','ITEM: TIMESTEP');
fprintf(myOutputFile1,'%s\r\n','35800000');
fprintf(myOutputFile1,'%s\r\n','ITEM: NUMBER OF ATOMS');
fprintf(myOutputFile1,'%d\r\n',length(arrset));
fprintf(myOutputFile1,'%s\r\n','ITEM: BOX BOUNDS ff ff ff');
fprintf(myOutputFile1,'%6.6f %6.6f\r\n',min(inf(:,6)),max(inf(:,6)));
fprintf(myOutputFile1,'%6.6f %6.6f\r\n',min(inf(:,7)),max(inf(:,7)));
fprintf(myOutputFile1,'%6.6f %6.6f\r\n',min(inf(:,8)),max(inf(:,8)));
fprintf(myOutputFile1,'%s\r\n','ITEM: ATOMS id x y z shapea shapeb shapec quat1 quat2 quat3 quat4 radius');
d=1;
chi = [];
tic
for i= 1:length(arrset);
    sp=inf(d:(d+length(arrset{i})-1),6:8);
    [center, radii, evecs, v, chi2] = ellipsoid_fit(sp, '' );
%     plot_Ellipsoid(sp,v);
%     hold on
    eul=rotm2eul(evecs);
    quat = eul2quat(eul);
%     quat = rotm2quat(evecs);
%     quat2 = dcm2quat(inv(evecs));
%     quat = rotm2quat_zyx(evecs);
    radius = (radii(1)*radii(2)*radii(3)).^(1/3);
    fprintf(myOutputFile1,'%d %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f %6.8f\r\n',i,center(1),center(2),center(3),radii(1),radii(2),radii(3),quat(1),quat(2),quat(3),quat(4),radius);
    d=d+length(arrset{i});
%     chi = [chi;chi2];
end
toc
fclose all;

function plot_Ellipsoid(coor,v)
    mind = min( coor );
    maxd = max( coor );
    nsteps = 10;
    step = ( maxd - mind ) / nsteps;
    [ x, y, z ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
    Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
          2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
          2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
    p = patch( isosurface( x, y, z, Ellipsoid, -v(10) ) );
    set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
end

function quat = rotm2quat_zyx(rotm)
    tr = [rotm(1,1),rotm(2,2),rotm(3,3)];
    if 1 + sum(tr) > 0.0004;
        q0 = sqrt(1 + rotm(1,1) + rotm(2,2) + rotm(3,3))/2;
        q1 = (rotm(3,2) - rotm(2,3))/(4*q0);
        q2 = (rotm(1,3) - rotm(3,1))/(4*q0);
        q3 = (rotm(2,1) - rotm(1,2))/(4*q0);
    elseif find(tr == max(tr)) == 1;
        t = sqrt(1 + rotm(1,1) - rotm(2,2) - rotm(3,3));
        q0 = (rotm(3,2) - rotm(2,3))/t;
        q1 = t/4;
        q2 = (rotm(1,3) + rotm(3,1))/t;
        q3 = (rotm(1,2) + rotm(2,1))/t;
    elseif find(tr == max(tr)) == 2;
        t = sqrt(1 - rotm(1,1) + rotm(2,2) - rotm(3,3));
        q0 = (rotm(1,3) - rotm(3,1))/t;
        q1 = (rotm(1,2) + rotm(2,1))/t;
        q2 = t/4;
        q3 = (rotm(3,2) + rotm(2,3))/t;
    elseif find(tr == max(tr)) == 3;
        t = sqrt(1 - rotm(1,1) - rotm(2,2) + rotm(3,3));
        q0 = (rotm(2,1) - rotm(1,2))/t;
        q1 = (rotm(1,3) + rotm(3,1))/t;
        q2 = (rotm(2,3) - rotm(3,2))/t;
        q3 = t/4;
    end
    quat = [q0,q1,q2,q3];
end
    
    