function b = Dominates(x, y)

    if isfield(x, 'Cost')
        x = x.Cost;
    end

    if isfield(y, 'Cost')
        y = y.Cost;
    end
    
    b = all(x <= y) && any(x<y);

end