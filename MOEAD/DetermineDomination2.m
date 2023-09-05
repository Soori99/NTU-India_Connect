function pop = DetermineDomination2(pop)

    nPop = numel(pop);

    for i = 1:nPop
        pop(i).IsDominated = false;
    end
    
    for i = 1:nPop
        for j = i+1:nPop
            if Dominates(pop(i), pop(j))
                pop(j).IsDominated = true;
                
            elseif Dominates(pop(j), pop(i))
                pop(i).IsDominated = true;
                
            end
        end
    end

end