clear
path(path,'C:\Users\Gao\Desktop\后处理代码\LIGGGHTS')
cNum = [];
cNum_f = [];
cNum_s = [];
cNum_w = [];
a_c = [];
a_fn = [];
a_ft = [];
a_d = [];
cNormF =  [];
cShearF = [];
cNormF_s =  [];
cShearF_s =  [];
cNormF_w =  [];
cShearF_w =  [];
fricIndex = [];
fricIndex_s = [];
fricIndex_w = [];
slideFrac = [];
p = zeros(51,1);
q = zeros(51,1);
fileFolder = [pwd,'\post1'];
dirOutput=dir(fullfile(fileFolder,'*.prop_cforce'));
fileNames={dirOutput.name}';
js = 0;
temp = [pwd,'\output.dat'];
Vol=textread(temp,'','headerlines',1);
Vol = [Vol(1,:);Vol];
for f = 1:length(fileNames);
    if isempty(strfind(fileNames{f},'30000.prop'));
        continue
    end
    js = js + 1;
    temp=[fileFolder,'\',fileNames{f}];
    ffid = textread(temp,'','headerlines',9);
    zhi = [ffid(:,1),ffid(:,2),ffid(:,3)] - [ffid(:,4),ffid(:,5),ffid(:,6)];
    zhi = [zhi;-zhi];
    force = ffid(:,10:12);
    CN = [];
    CN0 = [];
    for j = 1:length(force);
        CN(j,1:3) = (sum(force(j,:).*zhi(j,:))/(zhi(j,1).^2+zhi(j,2).^2+zhi(j,3).^2))*zhi(j,:);
        CN(j+length(force),1:3) = -CN(j,1:3);
        CN0(j,1) = (sum(force(j,:).*zhi(j,:))/sqrt(zhi(j,1).^2+zhi(j,2).^2+zhi(j,3).^2));
        CN0(j+length(force),1) = CN0(j,1);
    end
    force = [force;-force];
    CS = force - CN;
    CS0 =(CS(:,1).^2 + CS(:,2).^2 + CS(:,3).^2).^0.5;
    inf = [CN0,CN(:,1:3),CS0,CS(:,1:3)];
    fric = 0.5;
    cNum2 = 0;
    cNum1_s = 0;
    cNum1_w = 0;
    cNormF1 = 0;
    cShearF1 = 0;
    cNormF1_s = 0;
    cShearF1_s = 0;
    cNormF1_w = 0;
    cShearF1_w = 0;
    fricIndex1 = 0;
    fricIndex1_s = 0;
    fricIndex1_w = 0;
    slideFrac1 = 0;
    a_c1 = 0;
    a_fn1 = 0;
    a_ft1 = 0;
    a_d1 = 0;
    cFabric = zeros(3,3);
    fnFabric = zeros(3,3);
    ftFabric = zeros(3,3);
    dFabric = zeros(3,3);
    sigema = zeros(3,3);
    [s1,s2]=size(inf);
    cNum1 = s1;
    for k=1:s1;
        cnormfm = inf(k,1);
        cnormf1 = inf(k,2);
        cnormf2 = inf(k,3);
        cnormf3 = inf(k,4);
        cshearfm = inf(k,5);
        cshearf1 = inf(k,6);
        cshearf2 = inf(k,7);
        cshearf3 = inf(k,8);
        cNormF1 = cNormF1 + cnormfm;
        cShearF1 = cShearF1 + cshearfm;
        temp = cshearfm/(fric*cnormfm);
        if temp > 1.0;
            temp = 1.0;
        end
        if (cshearfm > cnormfm*fric);
            cNum2 = cNum2 + 1.0;
        end
        fricIndex1 =fricIndex1 + temp;
        lmn = [cnormf1/cnormfm, cnormf2/cnormfm, cnormf3/cnormfm];
        for i=1:3;
            for j=1:3;
                cFabric(i,j) = cFabric(i,j) + lmn(i)*lmn(j);
            end
        end
    end
    cFabric = cFabric/cNum1;
    cFabric_d = cFabric - eye(3)/3.0;
    aFabric_c = 15./2*cFabric_d;
    cNormF1 = cNormF1/cNum1;
    cShearF1 = cShearF1/cNum1;
    fricIndex1 = fricIndex1/cNum1;
    slideFrac1 = cNum2/cNum1;
    for k=1:s1;
        cnormfm = inf(k,1);
        cnormf1 = inf(k,2);
        cnormf2 = inf(k,3);
        cnormf3 = inf(k,4);
        cshearfm = inf(k,5);
        cshearf1 = inf(k,6);
        cshearf2 = inf(k,7);
        cshearf3 = inf(k,8);
        if cnormfm > cNormF1;
            temp = cshearfm/(fric*cnormfm);
            if temp > 1.0;
                temp = 1.0;
            end
            cNormF1_s = cNormF1_s + cnormfm;
            cShearF1_s = cShearF1_s + cshearfm;
            fricIndex1_s = fricIndex1_s + temp;
            cNum1_s = cNum1_s + 1;
        else
            cNormF1_w = cNormF1_w + cnormfm;
            cShearF1_w = cShearF1_w + cshearfm;
            fricIndex1_w = fricIndex1_w + temp;
            cNum1_w = cNum1_w + 1;
        end
        lmn = [cnormf1/cnormfm, cnormf2/cnormfm, cnormf3/cnormfm];
        lmn1 = lmn';
        lmn2 = lmn;
        deta = 1 + lmn2*(aFabric_c*lmn1);
        for i=1:3;
            for j=1:3;
                fnFabric(i,j) = fnFabric(i,j) + cnormfm*lmn(i)*lmn(j)/deta;
                sigema(i,j) = sigema(i,j) + (inf(k,i+1)+inf(k,i+5))*(zhi(k,j));
            end
        end
        if cshearfm >0 && (cnormf1*cshearf1 + cnormf2*cshearf2 + cnormf3*cshearf3)/(cnormfm*cshearfm)<0.034899 && (cnormf1*cshearf1 + cnormf2*cshearf2 + cnormf3*cshearf3)/(cnormfm*cshearfm)>-0.034899;
            lmt = [cshearf1/cshearfm, cshearf2/cshearfm, cshearf3/cshearfm];
            for i=1:3;
                for j=1:3;
                    ftFabric(i,j) = ftFabric(i,j) + cshearfm*lmt(i)*lmn(j)/deta;
                end
            end
        end
    end
    [s1,s2]=size(zhi);
    for k=1:s1;
        zhim=sqrt(zhi(k,1)^2+zhi(k,2)^2+zhi(k,3)^2);
        zhix=zhi(k,1);
        zhiy=zhi(k,2);
        zhiz=zhi(k,3);
        lmd = [zhix/zhim, zhiy/zhim, zhiz/zhim];
        lmn = [inf(k,2)/inf(k,1), inf(k,3)/inf(k,1), inf(k,4)/inf(k,1)];
        lmn1 = lmn';
        lmn2 = lmn;
        deta = 1 + lmn2*(aFabric_c*lmn1);
        for i=1:3;
            for j=1:3;
                dFabric(i,j) = dFabric(i,j) + zhim*lmd(i)*lmd(j)/deta;
            end
        end
    end
    %sig=dlmread('Job-ApplyLoad-Strain.dat');
    %sigema = [0.5,0,0;0,sig(f+3,4),0;0,0,0.5];
    sigema = sigema/(2*Vol(f,2));
    sigema0 = (sigema(1,1) + sigema(2,2) + sigema(3,3))/3;
    sigema_d = sigema - sigema0*eye(3);
    p(js) = sigema0/10^6;
    q(js) = sqrt(1.5*sum(sum(sigema_d.^2)))/10^6;
    fricIndex1_s = fricIndex1_s/cNum1_s;
    fricIndex1_w = fricIndex1_w/cNum1_w;
    cNormF1_s = cNormF1_s/cNum1_s;
    cShearF1_s = cShearF1_s/cNum1_s;
    cNormF1_w = cNormF1_w/cNum1_w;
    cShearF1_w = cShearF1_w/cNum1_w;
    fnFabric = fnFabric/cNum1;
    ftFabric = ftFabric/cNum1;
    dFabric = dFabric/s1;
    fn0 = fnFabric(1,1) + fnFabric(2,2) + fnFabric(3,3);
    ft0 = ftFabric(1,1) + ftFabric(2,2) + ftFabric(3,3);
    d0 = dFabric(1,1) + dFabric(2,2) + dFabric(3,3);
    fnFabric_d = fnFabric - fn0/3*eye(3);
    ftFabric_d = ftFabric - ft0/3*eye(3);
    dFabric_d = dFabric - d0/3*eye(3);
    aFabric_fn = 15./2*fnFabric_d/fn0;
    aFabric_ft = 15./3*ftFabric_d/fn0;
    aFabric_d = 15./2*dFabric_d/d0;
        
    a_c1 = asign(aFabric_c,sigema_d)*sqrt(1.5*sum(sum(aFabric_c.^2)));
    a_fn1 = asign(aFabric_fn,sigema_d)*sqrt(1.5*sum(sum(aFabric_fn.^2)));
    a_ft1 = asign(aFabric_ft,sigema_d)*sqrt(1.5*sum(sum(aFabric_ft.^2)));
    a_d1 = asign(aFabric_d,sigema_d)*sqrt(1.5*sum(sum(aFabric_d.^2)));
    
    cNum=[cNum,cNum1];
    cNum_f=[cNum_f,cNum2];
    cNum_s=[cNum_s,cNum1_s];
    cNum_w=[cNum_w,cNum1_w];
    cNormF=[cNormF,cNormF1];
    cShearF=[cShearF,cShearF1];
    cNormF_s=[cNormF_s,cNormF1_s];
    cShearF_s=[cShearF_s,cShearF1_s];
    cNormF_w=[cNormF_w,cNormF1_w];
    cShearF_w=[cShearF_w,cShearF1_w];
    fricIndex=[fricIndex,fricIndex1];
    fricIndex_s=[fricIndex_s,fricIndex1_s];
    fricIndex_w=[fricIndex_w,fricIndex1_w];
    slideFrac=[slideFrac,slideFrac1];
    a_c=[a_c,a_c1];
    a_fn=[a_fn,a_fn1];
    a_ft=[a_ft,a_ft1];
    a_d=[a_d,a_d1];
end    
ffid = fopen('Anisotropy.dat','w');
fprintf(ffid,'%s\r\n','Anisotropy evolution during loading precess');
fprintf(ffid,'%s\r\n','~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i=1:length(cNum);
    fprintf(ffid,'%d, %d, %d, %d, %d\r\n',i, cNum(i), cNum_s(i), cNum_w(i), cNum_f(i));
end
fprintf(ffid,'%s\r\n','~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i=1:length(cNum);
    fprintf(ffid,'%d, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f, %12.4f\r\n',i, cNormF(i), cShearF(i), cNormF_s(i), cShearF_s(i), cNormF_w(i), cShearF_w(i));
end
fprintf(ffid,'%s\r\n','~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i=1:length(cNum);
    fprintf(ffid,'%d, %12.4f, %12.4f, %12.4f, %12.4f\r\n',i, fricIndex(i), fricIndex_s(i), fricIndex_w(i), slideFrac(i));
end
fprintf(ffid,'%s\r\n','~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i=1:length(cNum);
    fprintf(ffid,'%d, %12.6f, %12.6f, %12.6f, %12.6f\r\n',i, a_c(i), a_fn(i), abs(a_ft(i)),a_d(i));
end
fclose all;