function [y,z]=ell(m,n,a,b,c);
    y=((- a^2*b^2*c^2 + a^2*b^2*n^2 + b^2*c^2*m^2)/(a^2*(b^2 - c^2)))^(1/2);
    z=(-(c^2*(- a^2*b^2 + a^2*n^2 + b^2*m^2))/(a^2*(b^2 - c^2)))^(1/2);
    y=real(y);
    z=real(z);
end