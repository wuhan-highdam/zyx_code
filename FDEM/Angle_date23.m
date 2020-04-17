bin=5;
jiao=xlsread('Angle_date.xlsx');
s=[0.01222	0.03668	0.06015	0.08048	0.09598	0.1063	0.1121	0.1146	0.1146	0.1133	0.1111	0.1087	0.1062	0.104	0.1023	0.1018	0.1019	0.1037]';
angle_23=zeros(90/bin,90/bin);
Fn_23=zeros(90/bin,90/bin);
for i = 1:length(jiao);
    for j=1 : 90/bin;
        for k=1 : 90/bin;
            if jiao(i,2)>(j-1)*bin && jiao(i,2)<=j*bin && jiao(i,3)>(k-1)*bin && jiao(i,3)<=k*bin;
                angle_23(j,k)=angle_23(j,k)+1;
                Fn_23(j,k)=Fn_23(j,k)+jiao(i,5);
            end
            if jiao(i,2)==0 && jiao(i,3)>(k-1)*bin && jiao(i,3)<=k*bin && j==1;
                angle_23(1,k)=angle_23(1,k)+1;
                Fn_23(1,k)=Fn_23(1,k)+jiao(i,5);
            end
            if jiao(i,2)>(j-1)*bin && jiao(i,2)<=j*bin && jiao(i,3)==0 && k==1;
                angle_23(j,1)=angle_23(j,1)+1;
                Fn_23(j,1)=Fn_23(j,1)+jiao(i,5);
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
for j = 1:90/bin;
    angle_23(j,:)=angle_23(j,:)./s(j).*mean(s);
end
myOutputFile2=fopen('loc_num.dat','w');
fprintf(myOutputFile2,'%s\r\n','TITLE=Rose Diagram');
fprintf(myOutputFile2,'%s\r\n','VARIABLES=''X'',''Y'',''Nun'',''Fn''');
fprintf(myOutputFile2,'%s %8.0f, %s %8.0f\r\n','ZONE N=    ',4*(90/bin)^2 ,'E=   ',(90/bin)^2);
fprintf(myOutputFile2,'%s\r\n','DATAPACKING=BLOCK,ZONETYPE=FEQUADRILATERAL');
fprintf(myOutputFile2,'%s\r\n','VARLOCATION=([3-4]=CELLCENTERED)');
for j = 1:90/bin;
    for k = 1:90/bin;
        fprintf(myOutputFile2,'%12.4f\r\n',(j-1)*bin);
        fprintf(myOutputFile2,'%12.4f\r\n',j*bin);
        fprintf(myOutputFile2,'%12.4f\r\n',j*bin);
        fprintf(myOutputFile2,'%12.4f\r\n',(j-1)*bin);
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
for j = 1:90/bin;
    for k = 1:90/bin;
        fprintf(myOutputFile2,'%12.4f\r\n',(k-1)*bin);
        fprintf(myOutputFile2,'%12.4f\r\n',(k-1)*bin);
        fprintf(myOutputFile2,'%12.4f\r\n',k*bin);
        fprintf(myOutputFile2,'%12.4f\r\n',k*bin);
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
for j = 1:90/bin;
    for k = 1:90/bin;
        fprintf(myOutputFile2,'%12.4f\r\n',angle_23(j,k));
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
for j = 1:90/bin;
    for k = 1:90/bin;
        fprintf(myOutputFile2,'%12.4f\r\n',Fn_23(j,k));
    end
end
fprintf(myOutputFile2,'%12.4f\r\n','');
eid = 1;
for j = 1:90/bin;
    for k = 1:90/bin;
        node1 = (eid-1)*4 + 1;
        node2 = (eid-1)*4 + 2;
        node3 = (eid-1)*4 + 3;
        node4 = (eid-1)*4 + 4;
        eid = eid + 1;
        fprintf(myOutputFile2,'%8.0f, %8.0f, %8.0f, %8.0f\r\n',node1, node2, node3, node4);
    end
end
fclose all;