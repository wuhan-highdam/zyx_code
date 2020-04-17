clear
dirOutput=dir(fullfile(pwd,'*.prop'));
fileNames={dirOutput.name};
temp = fileNames{1};
ffid = textread(temp,'','headerlines',9);
board_x1=min(ffid(:,3))-0.0034;
board_x2=max(ffid(:,3))+0.0034;
board_y1=min(ffid(:,4))-0.0034;
board_y2=max(ffid(:,4))+0.0034;
board_z1=min(ffid(:,5))-0.0034;
board_z2=max(ffid(:,5))+0.0034;
dirOutput=dir(fullfile(pwd,'*.prop_cforce'));
fileNames={dirOutput.name};
temp = fileNames{1};
ffid = fopen(fileNames{1},'r+');
dirOutput=dir(fullfile(pwd,'*.dat'));
adress=pwd;
rolling_faction=adress(17:end);
rolling_faction=strrep(rolling_faction,'\','_');
filename=['Force_chain_',rolling_faction,'.dat'];
fid = fopen(filename,'w');
while feof(ffid) == 0;
    tline1 = fgetl(ffid);
    if ~isempty(strfind(tline1,'ITEM: BOX BOUNDS mm mm mm'));
        fprintf(fid,'%s\r\n','ITEM: BOX BOUNDS pp pp pp');
        fprintf(fid,'%3.7f %3.7f\r\n',board_x1,board_x2);
        fprintf(fid,'%3.7f %3.7f\r\n',board_y1,board_y2);
        fprintf(fid,'%3.7f %3.7f\r\n',board_z1,board_z2);
        tline1 = fgetl(ffid);
        tline1 = fgetl(ffid);
        tline1 = fgetl(ffid);
    else
        fprintf(fid,'%s\r\n',tline1);
    end
end
fclose all;
