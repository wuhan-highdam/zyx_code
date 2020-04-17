clear
path(path,'C:\Users\Gao\Desktop\后处理代码\LIGGGHTS')
fileFolder = [pwd,'\post1'];
dirOutput=dir(fullfile(fileFolder,'*.prop_cforce'));
fileNames={dirOutput.name}';
js = 0;
temp = [pwd,'\output.dat'];
fn_av = [];
ft_av = [];
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
    fn_av = [fn_av;mean(CN0)];
    ft_av = [ft_av;mean(CS0)];
end
