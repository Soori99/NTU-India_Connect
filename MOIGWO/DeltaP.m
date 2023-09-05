function Score = DeltaP(PopObj,PF)
    IGDp  = mean(min(pdist2(PF,PopObj),[],2));
    GDp   = mean(min(pdist2(PopObj,PF),[],2));
    Score = max(IGDp,GDp);
end