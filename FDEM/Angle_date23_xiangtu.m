clear
bin=5;
jiao=xlsread('Angle_date4.xlsx');
ffid = dlmread('ELE_NOD2.inp');
node=ffid(1:361,2:3);
ele=ffid(362:end,2:5);
angle_23=zeros(90/bin,90/bin);
Fn_23=zeros(90/bin,90/bin);
for i = 1:length(jiao);
    for j=1 : 90/bin;
        for k=1 : 90/bin;
            if jiao(i,2)>(j-1)*bin && jiao(i,2)<=j*bin && jiao(i,3)>(k-1)*bin && jiao(i,3)<=k*bin;
                angle_23(j,k)=angle_23(j,k)+1;
                Fn_23(j,k)=Fn_23(j,k)+jiao(i,6);
            end
            if jiao(i,2)==0 && jiao(i,3)>(k-1)*bin && jiao(i,3)<=k*bin && j==1;
                angle_23(1,k)=angle_23(1,k)+1;
                Fn_23(1,k)=Fn_23(1,k)+jiao(i,6);
            end
            if jiao(i,2)>(j-1)*bin && jiao(i,2)<=j*bin && jiao(i,3)==0 && k==1;
                angle_23(j,1)=angle_23(j,1)+1;
                Fn_23(j,1)=Fn_23(j,1)+jiao(i,6);
            end
        end
    end
end
for j = 1:90/bin;
    for k = 1:90/bin;
        if angle_23(j,k)>0;
            Fn_23(j,k)=Fn_23(j,k)/angle_23(j,k);
        else
            Fn_23(j,k)=0;
        end
    end
end
myOutputFile2=fopen('loc_num3.dat','w');
fprintf(myOutputFile2,'%s\r\n','TITLE=Rose Diagram');
fprintf(myOutputFile2,'%s\r\n','VARIABLES=''X'',''Y'',''Nun'',''Fn''');
fprintf(myOutputFile2,'%s %8.0f, %s %8.0f\r\n','ZONE N=    ',length(node) ,'E=   ',(90/bin)^2);
fprintf(myOutputFile2,'%s\r\n','DATAPACKING=BLOCK,ZONETYPE=FEQUADRILATERAL');
fprintf(myOutputFile2,'%s\r\n','VARLOCATION=([3-4]=CELLCENTERED)');
for j =1 : length(node);
    fprintf(myOutputFile2,'%12.4f\r\n',node(j,1));
end
fprintf(myOutputFile2,'%12.4f\r\n','');
for j =1 : length(node);
    fprintf(myOutputFile2,'%12.4f\r\n',node(j,2));
end
fprintf(myOutputFile2,'%12.4f\r\n','');
for k= 1:90/bin;
    for j = 1:90/bin;
        fprintf(myOutputFile2,'%12.4f\r\n',angle_23(j,k));
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
for k = 1:90/bin;
    for j = 1:90/bin;
        fprintf(myOutputFile2,'%12.4f\r\n',Fn_23(j,k));
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
eid = 1;
for j = 1 : length(ele);
    fprintf(myOutputFile2,'%8.0f, %8.0f, %8.0f, %8.0f\r\n',ele(j,1), ele(j,2), ele(j,3),ele(j,4));
end
fclose all;