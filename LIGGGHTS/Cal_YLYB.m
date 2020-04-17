clear
fileFolder = [pwd,'\post'];
dirOutput=dir(fullfile(fileFolder,'*.prop'));
fileNames={dirOutput.name}';
contnum = [];
js = 0;
temp = [pwd,'\output.dat'];
inf = textread(temp,'','headerlines',1);
inf = [inf(1,:);inf];
inf0 = inf(1,:);
for i = 1:length(fileNames);
    if isempty(strfind(fileNames{i},'30000.prop'));
        continue
    end
    js = js + 1;
    z_yb(js) = -log(1-(inf0(5)-inf(i,5))/inf0(5));
    sigema_d(js) = (inf(i,8)-inf(i,6))/10^6;
    v_yb(js) = log(1+(inf(i,2)-inf0(2))/inf0(2));
end
z_yb = z_yb';
sigema_d = sigema_d';
v_yb = v_yb';
