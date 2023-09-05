function z = ZDT(x)

    n = numel(x);

    f1 = x(1);
    
    g = 1+9/(n-1)*sum(x(2:end));
    
    h = 1-sqrt(f1/g);
    
    f2 = g*h;
    
    z = [f1
       f2];

end