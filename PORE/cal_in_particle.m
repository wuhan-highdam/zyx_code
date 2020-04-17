function result = cal_in_particle(v,x,y,z);
    result = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
              2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
              2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z + v(10);
    if v(1) > 0;
        if result < 0;
            result = false;
        else
            result = true;
        end
    else
        if result >= 0;
            result = false;
        else
            result = true;
        end
    end
end
