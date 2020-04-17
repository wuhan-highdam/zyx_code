clear
dirOutput=dir(fullfile(pwd,'*.dat'));
fileNames={dirOutput.name}';
for f = 1:length(fileNames);
    if isempty(strfind(fileNames{f},'CNS_frame_'));
        continue
    end
    inf=dlmread(fileNames{f});
end
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
zNormF=sum(sum(cNormF))/sum(sum(cNormD))
