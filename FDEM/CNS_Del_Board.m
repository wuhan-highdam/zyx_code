clear
inf=dlmread('inf_frame_5.dat');
loc=dlmread('board.dat');
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
cnx=[];
cny=[];
cnz=[];
cn=[];
csx=[];
csy=[];
csz=[];
cs=[];
for i=1:length(arrset);
    sp=inf(d:(d+length(arrset{i})-1),1:12);
    x=[];
    y=[];
    z=[];
    CF=[];
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
    if isempty(x);
        jc(i)=0;
        d=d+length(arrset{i});
    elseif length(x)==1;
        if x<(loc(3)-0.0005) || x>(loc(4)+0.0005) || y<(loc(2)-0.0005) || y>(loc(1)+0.0005) || z<(loc(5)-0.0005) || z>(loc(6)+0.0005);
            jc(i)=0;
        else
            jc(i)=1; 
            cnx=[cnx,CF(3)];
            cny=[cny,CF(4)];
            cnz=[cnz,CF(5)];
            cn=[cn,CF(2)];
            csx=[csx,CF(10)];
            csy=[csy,CF(11)];
            csz=[csz,CF(12)];
            cs=[cs,CF(9)];
        end
        d=d+length(arrset{i});
    else
        A=[x',y',z'];
        dis=calDistance(A);
        T=DPSCAN(0.003,dis);
        for k=1:max(T);
            cla=find(T==k);
            for g=1:length(cla);
                if x(cla(g))<(loc(3)-0.0005) || x(cla(g))>(loc(4)+0.0005) || y(cla(g))<(loc(2)-0.0005) || y(cla(g))>(loc(1)+0.0005) || z(cla(g))<(loc(5)-0.0005) || z(cla(g))>(loc(6)+0.0005);
                    T(cla)=0;
                    break
                end
            end
        end
        if ~isempty(T);
            kz=0;
            for k=1:max(T);
                cla=find(T==k);
                if isempty(cla);
                    kz=kz+1;
                    continue
                end
                cnx=[cnx,sum(CF(cla,3))];
                cny=[cny,sum(CF(cla,4))];
                cnz=[cnz,sum(CF(cla,5))];
                cn=[cn,sqrt((sum(CF(cla,3)))^2+(sum(CF(cla,4)))^2+(sum(CF(cla,5)))^2)];
                csx=[csx,sum(CF(cla,10))];
                csy=[csy,sum(CF(cla,11))];
                csz=[csz,sum(CF(cla,12))];
                cs=[cs,sqrt((sum(CF(cla,10)))^2+(sum(CF(cla,11)))^2+(sum(CF(cla,12)))^2)];
            end
            jc(i)=max(T)-kz;
        end
        d=d+length(arrset{i});
    end
end
ff = fopen('CNS-3.dat','w');
for i=1:length(cn);
    fprintf(ff,'%12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f\r\n',cn(i),cnx(i),cny(i),cnz(i),cs(i),csx(i),csy(i),csz(i));
end
fclose all;
