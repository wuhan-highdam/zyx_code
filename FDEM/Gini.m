% Contributing author and copyright for this file:
% Wang Di 
% E-mail:
% wangdi1010@whu.edu.cn
% Reference:
% Hurley, R. C., Hall, S. A., & Wright, J. P. (2017). Multi-scale mechanics of granular solids from grain-resolved X-ray measurements. Proc. R. Soc. A, 473(2207), 20170491.
% Parameter:
% fn is a colume vector
% function:
% calculate Gini coefficient
function G=Gini(fn)
number = length(fn);
force= fn;
force_sort=sort(force);
sum1=0;sum2=0;
for i=1:number
    sum1 = (number+1-i)*force_sort(i)+sum1;
    sum2 = force_sort(i)+sum2;
end
G = (number+1-2*(sum1/sum2))/number;
end