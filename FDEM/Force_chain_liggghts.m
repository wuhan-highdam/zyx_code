tic
clear
path(path,'C:\Users\Gao\Desktop\后处理代码\FDEM');
dirOutput=dir(fullfile(pwd,'*.dat'));
adress=pwd;
shape=adress(12);
if isnan(str2double(adress(end-1)));
    frame=adress(end);
else
    frame=adress(end-1:end);
end
fileNames={dirOutput.name}';
for f = 1:length(fileNames);
    if isempty(strfind(fileNames{f},'inf_frame_'));
        continue
    end
    inf=dlmread(fileNames{f});
end
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
d=1;
jc=[];
ox=0;
oy=0;
oz=0;
oo=[];
for i=1:length(arrset);
    sp=inf(d:(d+length(arrset{i})-1),1:12);
    ox=mean(sp(:,6));
    oy=mean(sp(:,7));
    oz=mean(sp(:,8));
    oo=[oo;ox,oy,oz];
    x=[];
    y=[];
    z=[];
    CF=[];
    cnx=[];
    cny=[];
    cnz=[];
    cn=[];
    csx=[];
    csy=[];
    csz=[];
    cs=[];
    for j = 1 : length(arrset{i});
        if sp(j,2)==0;
            continue
        else
            x(end+1)=sp(j,6);
            y(end+1)=sp(j,7);
            z(end+1)=sp(j,8);
            CF=[CF;sp(j,:)];
        end
    end
    if length(x)==0;
        jc(i)=0;
        d=d+length(arrset{i});
    elseif length(x)==1
        jc(i)=1;
        d=d+length(arrset{i});
        cnx=[cnx;CF(3)];
        cny=[cny;CF(4)];
        cnz=[cnz;CF(5)];
        cn=[cn;CF(2)];
        csx=[csx;CF(10)];
        csy=[csy;CF(11)];
        csz=[csz;CF(12)];
        cs=[cs;CF(9)];
        ocn{i}=[cn,cnx,cny,cnz,cs,csx,csy,csz];
        os{i}=[x,y,z];
    else
        A=[x',y',z'];
        dis=calDistance(A);
        T=DPSCAN(0.003,dis);
        osx=[];
        osy=[];
        osz=[];
        for k=1:max(T);
            cla=find(T==k);
            cnx=[cnx;sum(CF(cla,3))];
            cny=[cny;sum(CF(cla,4))];
            cnz=[cnz;sum(CF(cla,5))];
            cn=[cn;sqrt((sum(CF(cla,3)))^2+(sum(CF(cla,4)))^2+(sum(CF(cla,5)))^2)];
            csx=[csx;sum(CF(cla,10))];
            csy=[csy;sum(CF(cla,11))];
            csz=[csz;sum(CF(cla,12))];
            cs=[cs;sqrt((sum(CF(cla,10)))^2+(sum(CF(cla,11)))^2+(sum(CF(cla,12)))^2)];
            osx=[osx;mean(x(cla))];
            osy=[osy;mean(y(cla))];
            osz=[osz;mean(z(cla))];
       end
        jc(i)=max(T);
        os{i}=[osx,osy,osz];
        ocn{i}=[cn,cnx,cny,cnz,cs,csx,csy,csz];
        d=d+length(arrset{i});
    end
end
ld=[];
lcn=[];
id=[];
pp=zeros(length(arrset),1);
for i=1:length(os);
    if isempty(os{i});
        continue
    end
    for m=1:length(os);
        if m==i || isempty(os{m});
           continue
        end
        if sqrt(sum((oo(i,:)-oo(m,:)).^2))<0.05;
            [s1,s2]=size(os{i});
            [s3,s4]=size(os{m});
            for j=1:s1;
                for n=1:s3;
                    if sqrt(sum((os{i}(j,:)-os{m}(n,:)).^2))<0.003;
                        ld=[ld;(oo(i,:)-oo(m,:))];
                        id=[id;oo(i,:),oo(m,:),i,m,ocn{i}(j,2)+ocn{i}(j,6),ocn{i}(j,3)+ocn{i}(j,7),ocn{i}(j,4)+ocn{i}(j,8),ocn{i}(j,1)];
                        lcn=[lcn;ocn{i}(j,:)];
                        break
                    end
                end
            end
        end
    end
end
lcn=sortrows(lcn,-1);
id=sortrows(id,-12);
a=[];
s1=length(lcn);
for i =2 : s1-1;
    if abs(lcn(i,1)-lcn(i-1,1)) > 0.01 && abs(lcn(i,1)-lcn(i+1,1)) > 0.01;
        a=[a;i];
    end
end
lcn(a,:)=[];
id(a,:)=[];
s1=length(lcn);
b=[];
for i =2 : s1-1;
    if abs(lcn(i,3)+lcn(i-1,3)) > 0.01 && abs(lcn(i,2)+lcn(i+1,2)) > 0.01;
        b=[b;i];
    end
end
lcn(b,:)=[];
id(b,:)=[];
s1=length(lcn);
c=[];
for i =2 : s1-1;
    if abs(lcn(i,5)-lcn(i-1,5)) > 0.01 && abs(lcn(i,5)-lcn(i+1,5)) > 0.01;
        c=[c;i];
    end
end
lcn(c,:)=[];
id(c,:)=[];
s1=length(lcn);
d=[];
for i =2 : s1-1;
    if abs(lcn(i,6)+lcn(i-1,6)) > 0.01 && abs(lcn(i,6)+lcn(i+1,6)) > 0.01;
        d=[d;i];
    end
end
lcn(d,:)=[];
id(d,:)=[];
s1=length(lcn);
d=[];
for i =2 : s1;
    if id(i,7)==id(i-1,8) && id(i,8)==id(i-1,7);
        d=[d;i];
    end
end
lcn(d,:)=[];
id(d,:)=[];

board_x1=min(inf(:,6));
board_x2=max(inf(:,6));
board_y1=min(inf(:,7));
board_y2=max(inf(:,7));
board_z1=min(inf(:,8));
board_z2=max(inf(:,8));
filename=['Force_chain_',num2str(shape),'_frame',num2str(frame),'.dat'];
ffid = fopen(filename,'w');
fprintf(ffid,'%s\r\n','ITEM: TIMESTEP');
fprintf(ffid,'%d\r\n',str2num(frame)*100000);
fprintf(ffid,'%s\r\n','ITEM: NUMBER OF ENTRIES');
fprintf(ffid,'%d\r\n',length(id));
fprintf(ffid,'%s\r\n','ITEM: BOX BOUNDS pp pp pp');
fprintf(ffid,'%3.7f %3.7f\r\n',board_x1,board_x2);
fprintf(ffid,'%3.7f %3.7f\r\n',board_y1,board_y2);
fprintf(ffid,'%3.7f %3.7f\r\n',board_z1,board_z2);
fprintf(ffid,'%s\r\n','ITEM: ENTRIES c_cforce[1] c_cforce[2] c_cforce[3] c_cforce[4] c_cforce[5] c_cforce[6] c_cforce[7] c_cforce[8] c_cforce[9] c_cforce[10] c_cforce[11] c_cforce[12] ');
for i = 1:length(id);
    fprintf(ffid,'%5.7f %5.7f %5.7f %5.7f %5.7f %5.7f %d %d %d %5.7f %5.7f %5.7f\r\n',id(i,1),id(i,2),id(i,3),id(i,4),id(i,5),id(i,6),id(i,7),id(i,8),0,id(i,9),id(i,10),id(i,11));
end
fclose all
toc
