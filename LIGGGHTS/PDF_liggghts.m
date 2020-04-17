clear
dirOutput=dir(fullfile(pwd,'*.prop_cforce'));
fileNames={dirOutput.name};
temp = fileNames{1};
ffid = textread(temp,'','headerlines',9);
zhi = [ffid(:,1),ffid(:,2),ffid(:,3)] - [ffid(:,4),ffid(:,5),ffid(:,6)];
force = ffid(:,10:12);
CN = [];
CN0 = [];
for i = 1:length(force);
    CN(i,1:3) = (sum(force(i,:).*zhi(i,:))/(zhi(i,1).^2+zhi(i,2).^2+zhi(i,3).^2))*zhi(i,:);
    CN(i+length(force),1:3) = -CN(i,1:3);
    CN0(i,1) = (sum(force(i,:).*zhi(i,:))/sqrt(zhi(i,1).^2+zhi(i,2).^2+zhi(i,3).^2));
    CN0(i+length(force),1) = CN0(i,1);
end
%CN = force.*zhi./((zhi(:,1).^2 + zhi(:,2).^2 + zhi(:,3).^2).^0.5);
force = [force;-force];
CS = force - CN;
CS0 =(CS(:,1).^2 + CS(:,2).^2 + CS(:,3).^2).^0.5;
inf = [CN0,CN(:,1),CN(:,3),CN(:,2),CS0,CS(:,1),CS(:,3),CS(:,2)];
force=inf(:,1);
s=length(force);
f=force/mean(force);
bin=0.1;
fmax=floor(10*max(force)/mean(force))/10;
num=roundn(fmax/bin,0);
pdf=zeros(num,1);
for i=1 : s;
    for j = 1:num;
        if f(i)>(j-1)*bin && f(i)<=j*bin
            pdf(j)=pdf(j)+1;
        end
    end
end
pdf=pdf/(s*bin);
cf=bin/2:bin:fmax-bin/2;
cf=cf';
kz=find(pdf==0);
pdf(kz)=[];
cf(kz)=[];
p=fittype('k.*x.^(c-1).*exp(-b^a*(abs(x-d)).^a)','independent','x','coefficients',{'a','c','k','b','d'});
st_=[1.5,0.7,1,0.7,0.1]; %设置初始值
lo_=[0.5,0,0,0,0]; %设置各参数下限
up_=[2.5,1.5,10,2,2]; %设置各参数上限
[cfun,gof]=fit(cf,pdf,p,'startpoint',st_,'lower',lo_,'upper',up_) %显示拟合函数
%pdf=pdf(1:100);
%cf=cf(1:100);
xi=cf;
yi=cfun(xi);
plot(cf,pdf,'r*',xi,yi,'b-');
