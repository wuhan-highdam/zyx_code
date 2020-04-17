clear
fileFolder = [pwd,'\post'];
dirOutput=dir(fullfile(fileFolder,'*.prop'));
fileNames={dirOutput.name}';
contnum = [];
js = 0;
for i = 1:length(fileNames);
    if isempty(strfind(fileNames{i},'30000.prop'));
        continue
    end
    js = js + 1;
    temp=[fileFolder,'\',fileNames{i}];
    ffid = textread(temp,'','headerlines',9);
    jc = ffid(:,18);
    z1 = length(find(jc==1));
    z0 = length(find(jc==0));
    contnum(js) = (sum(jc)-z1)/(length(jc)-z0-z1);
end
contnum = contnum';