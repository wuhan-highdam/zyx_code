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
    for j = 1 : length(arrset{i});
        if sp(j,2)==0;
            continue
        else
            x(end+1)=sp(j,6);
            y(end+1)=sp(j,7);
            z(end+1)=sp(j,8);
        end
    end
    if length(x)==0;
        jc(i)=0;
        d=d+length(arrset{i});
    elseif length(x)==1
        jc(i)=1;
        d=d+length(arrset{i});
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
            osx=[osx;mean(x(cla))];
            osy=[osy;mean(y(cla))];
            osz=[osz;mean(z(cla))];
       end
        jc(i)=max(T);
        os{i}=[osx,osy,osz];
        d=d+length(arrset{i});
    end
end
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
                        pp(i)=pp(i)+1;
                        break
                    end
                end
            end
        end
    end
end
pc = zeros(max(pp),1);
for i = 1 : length(pp);
    for j = 1:max(pp);
        if pp(i)==j;
            pc(j) = pc(j) + 1;
        end
    end
end
pc=[length(find(pp==0));pc];