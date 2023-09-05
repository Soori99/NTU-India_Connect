function g = DecomposedCost(individual, z, lambda)

    if isfield(individual, 'Cost')
        fx = individual.Cost;
    else
        fx = individual;
    end
    
    g = max(lambda.*abs(fx-z));

end