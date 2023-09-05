function Score = Spread(PopObj,PF)
    Dis1  = pdist2(PopObj,PopObj);
    Dis1(logical(eye(size(Dis1,1)))) = inf;
    [~,E] = max(PF,[],1);
    Dis2  = pdist2(PF(E,:),PopObj);
    d1    = sum(min(Dis2,[],2));
    d2    = mean(min(Dis1,[],2));
    Score = (d1+sum(abs(min(Dis1,[],2)-d2))) / (d1+(size(PopObj,1)-size(PopObj,2))*d2);
end