clear
dirOutput=dir(fullfile(pwd,'*.dat'));
fileNames={dirOutput.name}';
for f = 1:length(fileNames);
    if isempty(strfind(fileNames{f},'CNS_frame_'));
        continue
    end
    inf=dlmread(fileNames{f});
end
myOutputFile2=fopen('Comtact_Normal-Tecplot.dat','w');
[s1,s2]=size(inf);
fric = 0.5;
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
ang=12;
for i =1:ang;
    cNormD(1,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cNormD(1,nrow/ang*(i-1)+1:nrow/ang*i));
    cNormD(ncol,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cNormD(ncol,nrow/ang*(i-1)+1:nrow/ang*i));
    cNormF(1,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cNormF(1,nrow/ang*(i-1)+1:nrow/ang*i));
    cNormF(ncol,nrow/ang*(i-1)+1:nrow/ang*i)=mean(cNormF(ncol,nrow/ang*(i-1)+1:nrow/ang*i));
end
origin = [0,0,0];
angle_deg = angle*pi/180;
fprintf(myOutputFile2,'%s\r\n','TITLE=Rose Diagram');
fprintf(myOutputFile2,'%s\r\n','VARIABLES=''X'',''Y'',''Z'',''Fn'',''Contact_Normal''');
fprintf(myOutputFile2,'%s %8.0f, %s %8.0f\r\n','ZONE N=    ',(ncol*nrow*4+1),'E=   ',(ncol*nrow));
fprintf(myOutputFile2,'%s\r\n','DATAPACKING=BLOCK,ZONETYPE=FEBRICK');
fprintf(myOutputFile2,'%s\r\n','VARLOCATION=([4-5]=CELLCENTERED)');
fprintf(myOutputFile2,'%12.4f\r\n',origin(1));
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = cNormD(i,j)/abs(cos(theta - 0.5*angle_deg)-cos(theta + 0.5*angle_deg));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*cos(phi - 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*cos(phi + 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*cos(phi + 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*cos(phi - 0.5*angle_deg)));
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
fprintf(myOutputFile2,'%12.4f\r\n',origin(2));
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = cNormD(i,j)/abs(cos(theta - 0.5*angle_deg)-cos(theta + 0.5*angle_deg));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*sin(phi - 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta - 0.5*angle_deg)*sin(phi + 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*sin(phi + 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*sin(theta + 0.5*angle_deg)*sin(phi - 0.5*angle_deg)));
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
fprintf(myOutputFile2,'%12.4f\r\n',origin(3));
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = cNormD(i,j)/abs(cos(theta - 0.5*angle_deg)-cos(theta + 0.5*angle_deg));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*cos(theta - 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*cos(theta - 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*cos(theta + 0.5*angle_deg)));
        fprintf(myOutputFile2,'%12.4f\r\n',(rho1*cos(theta + 0.5*angle_deg)));
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
for i = 1:ncol;
    for j = 1:nrow;
        rho1 = cNormF(i,j);
        fprintf(myOutputFile2,'%12.4f\r\n',rho1);
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
for i = 1:ncol;
    for j = 1:nrow;
        theta = angle_deg*(i-0.5);
        phi = angle_deg*(j-0.5);
        rho1 = cNormD(i,j)/abs(cos(theta - 0.5*angle_deg)-cos(theta + 0.5*angle_deg));
        fprintf(myOutputFile2,'%12.4f\r\n',rho1);
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
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
    end
end
fclose all;
