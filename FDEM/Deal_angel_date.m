clear
bin=5;
jiao=xlsread('Angle_date2.xlsx');
angle_1=zeros(90/bin,1);
Fn_1=zeros(90/bin,1);
angle_2=zeros(90/bin,1);
Fn_2=zeros(90/bin,1);
angle_3=zeros(90/bin,1);
Fn_3=zeros(90/bin,1);
for i = 1:length(jiao);
    for j=1 : 90/bin;
        if jiao(i,1)>(j-1)*bin && jiao(i,1)<=j*bin;
            angle_1(j)=angle_1(j)+1;
            Fn_1(j)=Fn_1(j)+jiao(i,5);
        end
    end
    for j=1 : 90/bin;
        if jiao(i,2)>(j-1)*bin && jiao(i,2)<=j*bin;
            angle_2(j)=angle_2(j)+1;
            Fn_2(j)=Fn_2(j)+jiao(i,5);
        end
    end
    if jiao(i,2)==0;
        angle_2(1)=angle_2(1)+1;
        Fn_2(1)=Fn_2(1)+jiao(i,5);
    end
    for j=1 : 90/bin;
        if jiao(i,3)>(j-1)*bin && jiao(i,3)<=j*bin;
            angle_3(j)=angle_3(j)+1;
            Fn_3(j)=Fn_3(j)+jiao(i,5);
        end
    end
    if jiao(i,3)==0;
        angle_3(1)=angle_3(1)+1;
        Fn_3(1)=Fn_3(1)+jiao(i,5);
    end
end
Fn_1=Fn_1./angle_1;
Fn_2=Fn_2./angle_2;
Fn_3=Fn_3./angle_3;
pdf2=angle_2./sum(angle_2)./(bin*pi/180);
pdf3=angle_3./sum(angle_3)./(bin*pi/180);
