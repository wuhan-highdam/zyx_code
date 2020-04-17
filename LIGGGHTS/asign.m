function Sr = asign(x,y)
    [evecs,evals] =  eig(x);
    evals = sum(evals);
    T = find(evals==max(evals));
    x1 = evecs(:,T)';
    [evecs,evals] =  eig(y);
    evals = sum(evals);
    T = find(evals==max(evals));
    y1 = evecs(:,T);
    if x1*y1<1/3^0.5 && x1*y1>-1/3^0.5;
        Sr=-1;
    else
        Sr=1;
    end
end
