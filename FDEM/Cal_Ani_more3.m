clear
cNum = [];
cNum_f = [];
cNum_s = [];
cNum_w = [];
a_c=[];
a_fn=[];
a_ft=[];
a_cf=[];
a_bn=[];
a_bt=[];
a_cb=[];
a_fb=[];
a_s=[];
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
for f=0:50;
    filename1=['frame_all\CNS_frame_',num2str(f),'.dat'];
    inf=dlmread(filename1);
    filename2=['frame_all\ZHI_frame_',num2str(f),'.dat'];
    zhi=dlmread(filename2);
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
    a_bn1 = 0;
    a_bt1 = 0;
    cFabric = zeros(3,3);
    fnFabric = zeros(3,3);
    ftFabric = zeros(3,3);
    bnFabric = zeros(3,3);
    btFabric = zeros(3,3);
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
        lmt = [cshearf1/cshearfm, cshearf2/cshearfm, cshearf3/cshearfm];
        for i=1:3;
            for j=1:3;
                ftFabric(i,j) = ftFabric(i,j) + cshearfm*lmt(i)*lmn(j)/deta;
            end
        end
        %if cshearfm >0 && (cnormf1*cshearf1 + cnormf2*cshearf2 + cnormf3*cshearf3)/(cnormfm*cshearfm)<0.0087265 && (cnormf1*cshearf1 + cnormf2*cshearf2 + cnormf3*cshearf3)/(cnormfm*cshearfm)>-0.0087265;
         %   lmt = [cshearf1/cshearfm, cshearf2/cshearfm, cshearf3/cshearfm];
          %  for i=1:3;
           %     for j=1:3;
            %        sigema(i,j) = sigema(i,j) + (inf(k,i+1)+inf(k,i+5))*(zhi(k,j));
             %   end
           % end
        %end
    end
    [s1,s2]=size(zhi);
    for k=1:s1;
        zhim=sqrt(zhi(k,1)^2+zhi(k,2)^2+zhi(k,3)^2);
        zhix=zhi(k,1);
        zhiy=zhi(k,2);
        zhiz=zhi(k,3);
        lmd = [zhix/zhim, zhiy/zhim, zhiz/zhim];
        lmn = [zhi(k,5)/zhi(k,4), zhi(k,6)/zhi(k,4), zhi(k,7)/zhi(k,4)];
        bn = zhix*lmn(1)+zhiy*lmn(2)+zhiz*lmn(3);
        dmn = bn*lmn;
        dmt = [zhix,zhiy,zhiz]-dmn;
        bt = sqrt(sum(dmt.^2));
        lmdt = dmt./bt;
        lmn1 = lmn';
        lmn2 = lmn;
        deta = 1 + lmn2*(aFabric_c*lmn1);
        for i=1:3;
            for j=1:3;
                bnFabric(i,j) = bnFabric(i,j) + bn*lmn(i)*lmn(j)/deta;
                btFabric(i,j) = btFabric(i,j) + bt*lmdt(i)*lmn(j)/deta;
            end
        end
    end
    Vol=dlmread('Job-ApplyLoad-Volume.dat');
    %sig=dlmread('Job-ApplyLoad-Strain.dat');
    %sigema = [0.5,0,0;0,sig(f+3,4),0;0,0,0.5];
    sigema = sigema/Vol(f+3);
    sigema0 = (sigema(1,1) + sigema(2,2) + sigema(3,3))/3;
    sigema_d = sigema - sigema0*eye(3);
    sigema_d1 = sigema_d./sqrt(sum(sum(sigema_d.^2)));
    p(f+1) = sigema0/10^6;
    q(f+1) = sqrt(1.5*sum(sum(sigema_d.^2)))/10^6;
    fricIndex1_s = fricIndex1_s/cNum1_s;
    fricIndex1_w = fricIndex1_w/cNum1_w;
    cNormF1_s = cNormF1_s/cNum1_s;
    cShearF1_s = cShearF1_s/cNum1_s;
    cNormF1_w = cNormF1_w/cNum1_w;
    cShearF1_w = cShearF1_w/cNum1_w;
    fnFabric = fnFabric/cNum1;
    ftFabric = ftFabric/cNum1;
    bnFabric = bnFabric/s1;
    btFabric = btFabric/s1;
    fn0 = fnFabric(1,1) + fnFabric(2,2) + fnFabric(3,3);
    ft0 = ftFabric(1,1) + ftFabric(2,2) + ftFabric(3,3);
    bn0 = bnFabric(1,1) + bnFabric(2,2) + bnFabric(3,3);
    bt0 = btFabric(1,1) + btFabric(2,2) + btFabric(3,3);
    fnFabric_d = fnFabric - fn0/3*eye(3);
    ftFabric_d = ftFabric - ft0/3*eye(3);
    bnFabric_d = bnFabric - bn0/3*eye(3);
    btFabric_d = btFabric - bt0/3*eye(3);
    aFabric_fn = 15./2*fnFabric_d/fn0;
    aFabric_ft = 15./3*ftFabric_d/fn0;
    aFabric_bn = 15./2*bnFabric_d/bn0;
    aFabric_bt = 15./3*btFabric_d/bn0;
    aFabric_f = aFabric_fn - aFabric_ft;
    aFabric_b = aFabric_bn - aFabric_bt;
    aFabric_ftc = aFabric_ft*aFabric_c - sum(sum(aFabric_ft.*(aFabric_c')))*eye(3)/3;
    aFabric_cbt = aFabric_c*aFabric_bt - sum(sum(aFabric_c.*(aFabric_bt')))*eye(3)/3;
    aFabric_ftbt = aFabric_ft*aFabric_bt - sum(sum(aFabric_ft.*(aFabric_bt')))*eye(3)/3;
    aFabric_ftb = aFabric_ft*aFabric_b - sum(sum(aFabric_ft.*(aFabric_b')))*eye(3)/3;
    aFabric_fbt = aFabric_f*aFabric_bt - sum(sum(aFabric_f.*(aFabric_bt')))*eye(3)/3;
    aFabric_cf = (aFabric_c*aFabric_f + (aFabric_c*aFabric_f)')./2 - sum(sum(aFabric_c.*(aFabric_f')))*eye(3)/3;
    aFabric_cb = (aFabric_c*aFabric_b + (aFabric_c*aFabric_b)')./2 - sum(sum(aFabric_c.*(aFabric_b')))*eye(3)/3;
    aFabric_bf = (aFabric_b*aFabric_f + (aFabric_b*aFabric_f)')./2 - sum(sum(aFabric_b.*(aFabric_f')))*eye(3)/3;
    A_c = sqrt(1.5*sum(sum(aFabric_c.^2)));
    A_fn = sign(sum(sum(aFabric_c.*aFabric_fn)))*sqrt(1.5*sum(sum(aFabric_fn.^2)));
    A_ft = sign(sum(sum(aFabric_c.*aFabric_ft)))*sqrt(1.5*sum(sum(aFabric_ft.^2)));
    A_bn = sign(sum(sum(aFabric_c.*aFabric_bn)))*sqrt(1.5*sum(sum(aFabric_bn.^2)));
    A_bt = sign(sum(sum(aFabric_c.*aFabric_bt)))*sqrt(1.5*sum(sum(aFabric_bt.^2)));
    A_ftc = sign(sum(sum(aFabric_c.*aFabric_ftc)))*sqrt(1.5*sum(sum(aFabric_ftc.^2)));
    A_cbt = sign(sum(sum(aFabric_c.*aFabric_cbt)))*sqrt(1.5*sum(sum(aFabric_cbt.^2)));
    A_ftbt = sign(sum(sum(aFabric_c.*aFabric_ftbt)))*sqrt(1.5*sum(sum(aFabric_ftbt.^2)));
    A_ftb = sign(sum(sum(aFabric_c.*aFabric_ftb)))*sqrt(1.5*sum(sum(aFabric_ftb.^2)));
    A_fbt = sign(sum(sum(aFabric_c.*aFabric_fbt)))*sqrt(1.5*sum(sum(aFabric_fbt.^2)));
    A_cf = sign(sum(sum(aFabric_c.*aFabric_cf)))*sqrt(1.5*sum(sum(aFabric_cf.^2)));
    A_cb = sign(sum(sum(aFabric_c.*aFabric_cb)))*sqrt(1.5*sum(sum(aFabric_cb.^2)));
    A_bf = sign(sum(sum(aFabric_c.*aFabric_bf)))*sqrt(1.5*sum(sum(aFabric_bf.^2)));
    D0 = 1+4/45*(sign(sum(sum(aFabric_c.*aFabric_fn)))*abs(A_c*A_fn)+sign(sum(sum(aFabric_c.*aFabric_bn)))*abs(A_c*A_bn)+sign(sum(sum(aFabric_bn.*aFabric_fn)))*abs(A_bn*A_fn)+1.5*sign(sum(sum(aFabric_bt.*aFabric_ft)))*abs(A_bt*A_ft));
    
    a_c1 = 0.4*A_c/D0;
    a_fn1 = 0.4*A_fn/D0;
    a_ft1 = 0.6*A_ft/D0;
    a_cf1 = 2/35*(7*A_ftc+4*A_cf)/D0;
    a_bn1 = 0.4*A_bn/D0;
    a_bt1 = 0.6*A_bt/D0;
    a_cb1 = 2/35*(7*A_cbt+4*A_cb)/D0;
    a_fb1 = (A_ftbt + 0.4*(A_ftb +A_fbt)+8/35*A_bf)/D0;
    a_s1 = a_bn1 + a_bt1 + a_cb1 + a_fb1;
    
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
    a_cf=[a_cf,a_cf1];
    a_bn=[a_bn,a_bn1];
    a_bt=[a_bt,a_bt1];
    a_cb=[a_cb,a_cb1];
    a_fb=[a_fb,a_fb1];
    a_s=[a_s,a_s1];
end    
ffid = fopen('Anisotropy_new.dat','w');
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
    fprintf(ffid,'%d, %12.6f, %12.6f, %12.6f, %12.6f\r\n',i, a_bn(i), a_bt(i), a_cb(i),a_fb(i));
end
fprintf(ffid,'%s\r\n','~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for i=1:length(cNum);
    fprintf(ffid,'%d, %12.6f, %12.6f, %12.6f, %12.6f, %12.6f\r\n',i, a_c(i), a_fn(i), a_ft(i),a_cf(i),a_s(i));
end
fclose all;
a_c=a_c';
a_fn=a_fn';
a_ft=a_ft';
a_cf=a_cf';
a_bn=a_bn';
a_bt=a_bt';
a_cb=a_cb';
a_fb=a_fb';