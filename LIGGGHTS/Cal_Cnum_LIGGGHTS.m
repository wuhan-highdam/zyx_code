clear
dirOutput=dir(fullfile(pwd,'*.prop'));
fileNames={dirOutput.name};
temp = fileNames{1};
pp = textread(temp,'','headerlines',9);
pc = zeros(max(pp),1);
for i = 1 : length(pp);
    for j = 1:max(pp);
        if pp(i)==j;
            pc(j) = pc(j) + 1;
        end
    end
end
pc=[length(find(pp==0));pc];
