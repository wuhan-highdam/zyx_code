clear
path(path,'C:\Users\Gao\Desktop\后处理代码\FDEM')
GI=[];
for f=0:50;
    filename1=['frame_all\CNS_frame_',num2str(f),'.dat'];
    inf=dlmread(filename1);
    inf=sortrows(inf,-1);
    a=[];
    s1=length(inf);
    for i =2 : s1-1;
        if abs(inf(i,1)-inf(i-1,1)) > 0.01 && abs(inf(i,1)-inf(i+1,1)) > 0.01;
            a=[a;i];
        end
    end
    inf(a,:)=[];
    s1=length(inf);
    b=[];
    for i =2 : s1-1;
        if abs(inf(i,3)+inf(i-1,3)) > 0.01 && abs(inf(i,2)+inf(i+1,2)) > 0.01;
            b=[b;i];
        end
    end
    inf(b,:)=[];
    s1=length(inf);
    c=[];
    for i =2 : s1-1;
        if abs(inf(i,5)-inf(i-1,5)) > 0.01 && abs(inf(i,5)-inf(i+1,5)) > 0.01;
            c=[c;i];
        end
    end
    inf(c,:)=[];
    s1=length(inf);
    d=[];
    for i =2 : s1-1;
        if abs(inf(i,6)+inf(i-1,6)) > 0.01 && abs(inf(i,6)+inf(i+1,6)) > 0.01;
            d=[d;i];
        end
    end
    inf(d,:)=[];
    force=inf(:,1);
    %force=force(2:2:end);
    G=Gini(force);
    GI=[GI;G];
end
