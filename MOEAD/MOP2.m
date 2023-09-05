function z = MOP2(x)

    n = numel(x);
    
    z = [1-exp(-sum((x-1/sqrt(n)).^2))
       1-exp(-sum((x+1/sqrt(n)).^2))];
    
end