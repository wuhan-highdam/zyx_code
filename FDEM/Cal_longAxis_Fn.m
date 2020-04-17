clear
path(path,'C:\Users\Gao\Desktop\后处理代码\FDEM')
dirOutput=dir(fullfile(pwd,'*.dat'));
fileNames={dirOutput.name}';
for f = 1:length(fileNames);
    if isempty(strfind(fileNames{f},'inf_frame_'));
        continue
    end
    inf=dlmread(fileNames{f});
end
fclose all;
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
zx=[];
d=1;
jc=[];
ox=0;
oy=0;
oz=0;
oo=[];
jiao=[];
for i=1:length(arrset);
    sp=inf(d:(d+length(arrset{i})-1),1:12);
    ox=mean(sp(:,6));
    oy=mean(sp(:,7));
    oz=mean(sp(:,8));
    oo=[oo;ox,oy,oz];
    dis=[];
    for j=1:length(sp);
        dis(j)=sqrt(sum((sp(j,6:8)-[ox,oy,oz]).^2));
    end
    [dis1,index1]=sort(dis,'descend');
    index2=index1(1:floor(length(sp)/20));
    sp0=sp(index2,6:8);
    dis0=calDistance(sp0);
    [row,col]=find(dis0==max(max(dis0)));
    B = sp0(row(1),:)-sp0(row(2),:);
    zx=[zx;B];
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
                        lcn=[lcn;ocn{i}(j,:)];
                        xl1=zx(i,:);
                        xl2=ocn{i}(j,2:4);
                        jiajiao=acos(dot(xl1,xl2)/(norm(xl1)*norm(xl2)))*180/pi;
                        if jiajiao>90;
                            jiajiao=180-jiajiao;
                        end
                        jiao=[jiao;jiajiao];
                        break
                    end
                end
            end
        end
    end
end
tj=zeros(90,1);
for i=1:length(jiao);
    for j=1:90;
        if jiao(i)>=j-1 && jiao(i)<j;
            tj(j)=tj(j)+1;
        end
    end
end