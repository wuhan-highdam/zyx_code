clear
dirOutput=dir(fullfile(pwd,'*.prop_cforce'));
fileNames={dirOutput.name};
temp = fileNames{1};
ffid = textread(temp,'','headerlines',9);
zhi = [ffid(:,1),ffid(:,2),ffid(:,3)] - [ffid(:,4),ffid(:,5),ffid(:,6)];
force = ffid(:,10:12);
CN = [];
CN0 = [];
for i = 1:length(force);
    CN(i,1:3) = (sum(force(i,:).*zhi(i,:))/(zhi(i,1).^2+zhi(i,2).^2+zhi(i,3).^2))*zhi(i,:);
    CN(i+length(force),1:3) = -CN(i,1:3);
    CN0(i,1) = (sum(force(i,:).*zhi(i,:))/sqrt(zhi(i,1).^2+zhi(i,2).^2+zhi(i,3).^2));
    CN0(i+length(force),1) = CN0(i,1);
end
%CN = force.*zhi./((zhi(:,1).^2 + zhi(:,2).^2 + zhi(:,3).^2).^0.5);
force = [force;-force];
CS = force - CN;
CS0 =(CS(:,1).^2 + CS(:,2).^2 + CS(:,3).^2).^0.5;
inf = [CN0,CN(:,1),CN(:,3),CN(:,2),CS0,CS(:,1),CS(:,3),CS(:,2)];
myOutputFile1=fopen('CNormF.dat','w');
myOutputFile2=fopen('CNormF-Tecplot.dat','w');
myOutputFile3=fopen('CShearF.dat','w');
myOutputFile4=fopen('CShearF-Tecplot.dat','w');
myOutputFile5=fopen('ContactForce.dat','w');
myOutputFile6=fopen('RoseDiagrm.dat','w');
[s1,s2]=size(inf);
fric = 0.2;
angle = 10;
ncol = 180/angle;
nrow = 360/angle;
cNum = 0;
cNormD=zeros(ncol,nrow);
cNormF=zeros(ncol,nrow);
cShearD=zeros(ncol,nrow);
cShearF=zeros(ncol,nrow);
cNormD1=zeros(60,1);
cNormF1=zeros(60,1);
cShearD1=zeros(60,1);
cShearF1=zeros(60,1);
cFricI=zeros(60,1);
cNormF2=[];
cShearF2=[];
cFricMob=[];
for k=1:s1;
    cnormfm = inf(k,1);
    cnormf1 = inf(k,2);
    cnormf2 = inf(k,3);
    cnormf3 = inf(k,4);
    cnormf0 = sqrt(cnormf1^2 + cnormf3^2);
    cshearfm = inf(k,5);
    cshearf1 = inf(k,6);
    cshearf2 = inf(k,7);
    cshearf3 = inf(k,8);
    cshearf0 = sqrt(cshearf1^2 + cshearf3^2);
    if cnormf0 > 0;
        mobilization = abs(cshearfm)/(fric*cnormfm);
        cNum=cNum+1;
        cNormF2(cNum)=cnormfm;
        cShearF2(cNum)=cshearfm;
        cFricMob(cNum)=mobilization;
        theta1 = acos(cnormf2/cnormfm)*180/pi;
        if theta1 < 0;
            theta1 = 0;
        end
        if theta1 >= 180;
            theta1 = 180 - 0.1;
        end
        phi1 = cnormf3/cnormf0;
        if cnormf1 > 0;
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
        cNormD(i,j)=cNormD(i,j)+1;
        cNormF(i,j)=cNormF(i,j)+cnormfm;
    end
    if cshearfm > 0 && cshearf0 > 0 
        theta2 = acos(cshearf2/cshearfm)*180/pi;
        if theta2 < 0;
            theta2 = 0;
        end
        if theta2 >= 180;
            theta2 = 180 - 0.1;
        end
        phi2 = cshearf3/cshearf0;
        if cshearf1 > 0;
            phi2 = acos(phi2)*180/pi;
        else
            phi2 = 360 - acos(phi2)*180/pi;
        end
        if phi2 < 0;
             phi2 = 0;
        end
        if phi2 >= 360;
            phi2 = 360 - 0.1;
        end
        m = floor(theta2/angle)+1;
        n = floor(phi2/angle)+1;
        if m >= (ncol + 1);
            m = ncol - 1;
        end
        if n >= (nrow + 1);
            n = nrow - 1;
        end
        cShearD(m,n)=cShearD(m,n)+1;
        cShearF(m,n)=cShearF(m,n)+cshearfm;
    end
    if cshearfm > 0
        rad = cnormf0/cnormfm;
        if rad > 1.0
            rad = 1;
        end
        if cnormf1 > 0
            if cnormf2 > 0
                ang = acos(rad)*180/pi;
            else
                ang =  360 - acos(rad)*180/pi;
            end
        else
            if cnormf2 > 0
                ang = 180 - acos(rad)*180/pi;
            else
                ang = 180 + acos(rad)*180/pi;
            end
        end
        o = floor((ang + 3)/6.0)+1;
        if o >= 2 && o < 61;
            cNormD1(o)=cNormD1(o)+1;
            cNormF1(o)=cNormF1(o)+cnormfm;
            cFricI(o)=cFricI(o)+mobilization;
        else
            cNormD1(1)=cNormD1(1)+1;
            cNormF1(1)=cNormF1(1)+cnormfm;
            cFricI(1)=cFricI(1)+mobilization;
        end
        rad = cshearf0/cshearfm;
        if rad > 1.0
            rad = 1;
        end
        if cshearf1 > 0
            if cshearf2 > 0
                ang = acos(rad)*180/pi;
            else
                ang =  360 - acos(rad)*180/pi;
            end
        else
            if cshearf2 > 0
                ang = 180 - acos(rad)*180/pi;
            else
                ang = 180 + acos(rad)*180/pi;
            end
        end
        o = floor((ang + 3)/6.0)+1;
        if o >= 2 && o < 61;
            cShearD1(o)=cShearD1(o)+1;
            cShearF1(o)=cShearF1(o)+cshearfm;
        else
            cShearD1(1)=cShearD1(1)+1;
            cShearF1(1)=cShearF1(1)+cshearfm;
        end
    end
end
for i = 1:ncol;
    for j = 1:nrow;
        if cNormD(i,j)>0;
            cNormF(i,j)=cNormF(i,j)/cNormD(i,j);
        else
            cNormF(i,j)=0;
        end
    end
end
for i = 1:ncol;
    for j = 1:nrow;
        if cShearD(i,j)>0;
            cShearF(i,j)=cShearF(i,j)/cShearD(i,j);
        else
            cShearF(i,j)=0;
        end
    end
end
for i = 1:60;
    if cNormD1(i)>0
        cNormF1(i) = cNormF1(i)/cNormD1(i);
        cFricI(i) = cFricI(i)/cNormD1(i);
    else
        cNormF1(i) = 0;
        cFricI(i) = 0;
    end
    if cShearD1(i) > 0
        cShearF1(i) = cShearF1(i)/cShearD1(i);
    else
        cShearF1(i) = 0;
    end
end
for i = 1:cNum
    fprintf(myOutputFile5,'%12.4f,  %12.4f,  %12.4f\r\n',cNormF2(i),cShearF2(i),cFricMob(i));
end
for i = 1:length(cNormD1)
    fprintf(myOutputFile6,'%12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f\r\n',6*(i-1), cNormD1(i), cNormF1(i), cFricI(i), cShearD1(i), cShearF1(i));
end
for i = 1:ncol;
    for j = 1:nrow;
        fprintf(myOutputFile1,'%12.4f,  %12.4f,  %12.4f,  %12.4f\r\n',angle*(i-0.5),angle*(j-0.5),cNormD(i,j),cNormF(i,j));
        fprintf(myOutputFile3,'%12.4f,  %12.4f,  %12.4f,  %12.4f\r\n',angle*(i-0.5),angle*(j-0.5),cShearD(i,j),cShearF(i,j));
    end
end
ang=12;
for i =1:ang;
    cNormD(1,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cNormD(1,nrow/ang*(i-1)+1:nrow/ang*i));
    cNormD(ncol,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cNormD(ncol,nrow/ang*(i-1)+1:nrow/ang*i));
    cNormF(1,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cNormF(1,nrow/ang*(i-1)+1:nrow/ang*i));
    cNormF(ncol,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cNormF(ncol,nrow/ang*(i-1)+1:nrow/ang*i));
    cShearD(1,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cShearD(1,nrow/ang*(i-1)+1:nrow/ang*i));
    cShearD(ncol,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cShearD(ncol,nrow/ang*(i-1)+1:nrow/ang*i));
    cShearF(1,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cShearF(1,nrow/ang*(i-1)+1:nrow/ang*i));
    cShearF(ncol,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cShearF(ncol,nrow/ang*(i-1)+1:nrow/ang*i));
end
origin = [0,0,0];
angle_deg = angle*pi/180;
fprintf(myOutputFile2,'%s\r\n','TITLE=Rose Diagram');
fprintf(myOutputFile2,'%s\r\n','VARIABLES=''X'',''Y'',''Z'',''Fn''');
fprintf(myOutputFile2,'%s %8.0f, %s %8.0f\r\n','ZONE N=    ',(ncol*nrow*4+1),'E=   ',(ncol*nrow));
fprintf(myOutputFile2,'%s\r\n','DATAPACKING=BLOCK,ZONETYPE=FEBRICK');
fprintf(myOutputFile2,'%s\r\n','VARLOCATION=([4]=CELLCENTERED)');
fprintf(myOutputFile2,'%12.4f\r\n',origin(1));
fprintf(myOutputFile4,'%s\r\n','TITLE=Rose Diagram');
fprintf(myOutputFile4,'%s\r\n','VARIABLES=''X'',''Y'',''Z'',''Ft''');
fprintf(myOutputFile4,'%s %8.0f, %s %8.0f\r\n','ZONE N=    ',(ncol*nrow*4+1),'E=   ',(ncol*nrow));
fprintf(myOutputFile4,'%s\r\n','DATAPACKING=BLOCK,ZONETYPE=FEBRICK');
fprintf(myOutputFile4,'%s\r\n','VARLOCATION=([4]=CELLCENTERED)');
fprintf(myOutputFile4,'%12.4f\r\n',origin(1));
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = cNormF(i,j);
        rho2 = cShearF(i,j);
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*cos(phi - 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*cos(phi + 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*cos(phi + 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*cos(phi - 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*sin(theta - 0.5*angle_deg)*cos(phi - 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*sin(theta - 0.5*angle_deg)*cos(phi + 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*sin(theta + 0.5*angle_deg)*cos(phi + 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*sin(theta + 0.5*angle_deg)*cos(phi - 0.5*angle_deg)));
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
fprintf(myOutputFile2,'%12.4f\r\n',origin(2));
fprintf(myOutputFile4,'%12.4f\r\n','');
fprintf(myOutputFile4,'%12.4f\r\n',origin(2));
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = cNormF(i,j);
        rho2 = cShearF(i,j);
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*sin(phi - 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*sin(phi + 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*sin(phi + 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*sin(phi - 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*sin(theta - 0.5*angle_deg)*sin(phi - 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*sin(theta - 0.5*angle_deg)*sin(phi + 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*sin(theta + 0.5*angle_deg)*sin(phi + 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*sin(theta + 0.5*angle_deg)*sin(phi - 0.5*angle_deg)));
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
fprintf(myOutputFile2,'%12.4f\r\n',origin(3));
fprintf(myOutputFile4,'%12.4f\r\n','');
fprintf(myOutputFile4,'%12.4f\r\n',origin(3));
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = cNormF(i,j);
        rho2 = cShearF(i,j);
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*cos(theta - 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*cos(theta - 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*cos(theta + 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*cos(theta + 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*cos(theta - 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*cos(theta - 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*cos(theta + 0.5*angle_deg)));
        fprintf(myOutputFile4,'%12.4f\r\n',(rho2*cos(theta + 0.5*angle_deg)));
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
fprintf(myOutputFile4,'%12.4f\r\n','');
for i = 1:ncol;
    for j = 1:nrow;
        rho1 = cNormF(i,j);
        rho2 = cShearF(i,j);
        fprintf(myOutputFile2,'%12.4f\r\n',rho1);
        fprintf(myOutputFile4,'%12.4f\r\n',rho2);
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
fprintf(myOutputFile4,'%12.4f\r\n','');
oid = 1;
eid = 1;
for i = 1:ncol;
    for j = 1:nrow;
        node1 = (eid-1)*4 + 2;
        node2 = (eid-1)*4 + 3;
        node3 = (eid-1)*4 + 4;
        node4 = (eid-1)*4 + 5;
        eid = eid + 1;
        fprintf(myOutputFile2,'%8.0f, %8.0f, %8.0f, %8.0f, %8.0f, %8.0f, %8.0f, %8.0f\r\n',oid, oid, oid, oid, node1, node2, node3, node4);
        fprintf(myOutputFile4,'%8.0f, %8.0f, %8.0f, %8.0f, %8.0f, %8.0f, %8.0f, %8.0f\r\n',oid, oid, oid, oid, node1, node2, node3, node4);
    end
end
fclose all;
