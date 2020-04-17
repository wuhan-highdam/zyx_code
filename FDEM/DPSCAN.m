%% æ€¿‡À„∑®
function [class] = DPSCAN(r,dis)
x=dis(:,1);
dealed=zeros(length(x),1);
num=1;
class=zeros(length(x),1);
for i=1:length(x)
    if dealed(i) == 0
        dealed(i)=1;
        class(i)=num;
        D=dis(i,:);
    	ind=find(D<r);
        for j=1:length(ind);
            if dealed(ind(j))==1;
                ind(j)=0;
            end
        end
        ind(ind==0)=[];
        if ~isempty(ind);
            for k=1:length(ind);
                dealed(ind(k))=1;
                class(ind(k))=num;
            end
            while ~isempty(ind)
                D=dis(ind(1),:);
                ind_1=find(D<r);
                if ~isempty(ind_1);
                    for j=1:length(ind_1);
                        if dealed(ind_1(j))==1;
                           ind_1(j)=0;
                        end
                    end
                end
                ind_1(ind_1==0)=[];
                if ~isempty(ind_1);
                    for k=1:length(ind_1);
                        dealed(ind_1(k))=1;
                        class(ind_1(k))=num;
                    end
                    ind=[ind ind_1];
                    ind(1)=[];
                else
                    ind(1)=[];
                end
            end
            num=num+1;
        else
            num=num+1;
        end
    end
end
end

