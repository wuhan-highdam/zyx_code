clear
dirOutput=dir(fullfile(pwd,'*.prop_cforce'));
fileNames={dirOutput.name}';
temp=[pwd,'\',fileNames{1}];
ffid = textread(temp,'','headerlines',9);
zhi = [ffid(:,1),ffid(:,3),ffid(:,2)] - [ffid(:,4),ffid(:,6),ffid(:,5)];
zhi = [zhi;-zhi];
myOutputFile7=fopen('Branch1-Tecplot.dat','w');
[s1,s2]=size(zhi);
angle = 10;
ncol = 180/angle;
nrow = 360/angle;
cNum=0;
zhiD=zeros(ncol,nrow);
zhiF=zeros(ncol,nrow);
[s1,s2]=size(zhi);
for k=1:s1;
    zhim=sqrt(zhi(k,1)^2+zhi(k,2)^2+zhi(k,3)^2);
    zhi1=zhi(k,1);
    zhi2=zhi(k,2);
    zhi3=zhi(k,3);
    zhi0 = sqrt(zhi1^2 + zhi3^2);
    if zhi0 > 0;
        cNum=cNum+1;
        theta1 = acos(zhi2/zhim)*180/pi;
        if theta1 < 0;
            theta1 = 0;
        end
        if theta1 >= 180;
            theta1 = 180 - 0.1;
        end
        phi1 = zhi3/zhi0;
        if zhi1 > 0;
            phi1 = acos(phi1)*180/pi;
        else
            phi1 = 360 - acos(phi1)*180/pi;
        end
        if phi1 < 0;
             phi1 = 0;
        end
        if phi1 >= 360;
            phi1 = 360 - 0.1;
        end
        i = floor(theta1/angle)+1;
        j = floor(phi1/angle)+1;
        if i >= (ncol + 1);
            i = ncol - 1;
        end
        if j >= (nrow + 1);
            j = nrow - 1;
        end
        zhiD(i,j)=zhiD(i,j)+1;
        zhiF(i,j)=zhiF(i,j)+zhim;
    end
end
ang=12;
for i =1:ang;
    zhiD(1,nrow/ang*(i-1)+1:nrow/ang*i)=mean(zhiD(1,nrow/ang*(i-1)+1:nrow/ang*i));
    zhiD(ncol,nrow/ang*(i-1)+1:nrow/ang*i)=mean(zhiD(ncol,nrow/ang*(i-1)+1:nrow/ang*i));
    zhiF(1,nrow/ang*(i-1)+1:nrow/ang*i)=mean(zhiF(1,nrow/ang*(i-1)+1:nrow/ang*i));
    zhiF(ncol,nrow/ang*(i-1)+1:nrow/ang*i)=mean(zhiF(ncol,nrow/ang*(i-1)+1:nrow/ang*i));
end
for i = 1:ncol;
    for j = 1:nrow;
        if zhiD(i,j)>0;
            zhiF(i,j)=zhiF(i,j)/zhiD(i,j);
        else
            zhiF(i,j)=0;
        end
    end
end
origin = [0,0,0];
angle_deg = angle*pi/180;
fprintf(myOutputFile7,'%s\r\n','TITLE=Rose Diagram');
fprintf(myOutputFile7,'%s\r\n','VARIABLES=''X'',''Y'',''Z'',''Branch-size'',''Branch-density''');
fprintf(myOutputFile7,'%s %8.0f, %s %8.0f\r\n','ZONE N=    ',(ncol*nrow*4+1),'E=   ',(ncol*nrow));
fprintf(myOutputFile7,'%s\r\n','DATAPACKING=BLOCK,ZONETYPE=FEBRICK');
fprintf(myOutputFile7,'%s\r\n','VARLOCATION=([4-5]=CELLCENTERED)');
fprintf(myOutputFile7,'%12.4f\r\n',origin(1));
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = zhiD(i,j)/abs(cos(theta - 0.5*angle_deg)-cos(theta + 0.5*angle_deg));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*cos(phi - 0.5*angle_deg)));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*cos(phi + 0.5*angle_deg)));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*cos(phi + 0.5*angle_deg)));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*cos(phi - 0.5*angle_deg)));
    end
end
fprintf(myOutputFile7,'%12.4f\r\n','');
fprintf(myOutputFile7,'%12.4f\r\n',origin(2));
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = zhiD(i,j)/abs(cos(theta - 0.5*angle_deg)-cos(theta + 0.5*angle_deg));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*sin(phi - 0.5*angle_deg)));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*sin(phi + 0.5*angle_deg)));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*sin(phi + 0.5*angle_deg)));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*sin(phi - 0.5*angle_deg)));
    end
end
fprintf(myOutputFile7,'%12.4f\r\n','');
fprintf(myOutputFile7,'%12.4f\r\n',origin(3));
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = zhiD(i,j)/abs(cos(theta - 0.5*angle_deg)-cos(theta + 0.5*angle_deg));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*cos(theta - 0.5*angle_deg)));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*cos(theta - 0.5*angle_deg)));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*cos(theta + 0.5*angle_deg)));
        fprintf(myOutputFile7,'%12.4f\r\n',(rho1*cos(theta + 0.5*angle_deg)));
    end
end
fprintf(myOutputFile7,'%12.4f\r\n','');
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = zhiF(i,j);
        fprintf(myOutputFile7,'%12.4f\r\n',rho1);
    end
end
fprintf(myOutputFile7,'%12.4f\r\n','');
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = zhiD(i,j)/abs(cos(theta - 0.5*angle_deg)-cos(theta + 0.5*angle_deg));
        fprintf(myOutputFile7,'%12.4f\r\n',rho1);
    end
end
fprintf(myOutputFile7,'%12.4f\r\n','');
oid = 1;
eid = 1;
for i = 1:ncol;
    for j = 1:nrow;
        node1 = (eid-1)*4 + 2;
        node2 = (eid-1)*4 + 3;
        node3 = (eid-1)*4 + 4;
        node4 = (eid-1)*4 + 5;
        eid = eid + 1;
        fprintf(myOutputFile7,'%8.0f, %8.0f, %8.0f, %8.0f, %8.0f, %8.0f, %8.0f, %8.0f\r\n',oid, oid, oid, oid, node1, node2, node3, node4);
    end
end
fclose all;
