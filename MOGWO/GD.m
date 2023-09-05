function Score = GD(PopObj,PF)
    Distance = min(pdist2(PopObj,PF),[],2);
    Score    = norm(Distance) / length(Distance);
end