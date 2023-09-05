function Score = Spacing(PopObj,PF)
    Distance = pdist2(PopObj,PopObj,'cityblock');
    Distance(logical(eye(size(Distance,1)))) = inf;
    Score    = std(min(Distance,[],2));
end