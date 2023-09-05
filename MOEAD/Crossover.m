function y = Crossover(x1, x2, params)

    gamma = params.gamma;
    VarMin = params.VarMin;
    VarMax = params.VarMax;
    
    alpha = unifrnd(-gamma, 1+gamma, size(x1));
    
    y = alpha.*x1+(1-alpha).*x2;

    y = min(max(y, VarMin), VarMax);
    
end