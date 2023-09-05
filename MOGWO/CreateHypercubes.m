function G=CreateHypercubes(costs,ngrid,alpha)

    nobj=size(costs,1);
    
    empty_grid.Lower=[];
    empty_grid.Upper=[];
    G=repmat(empty_grid,nobj,1);
    
    for j=1:nobj
        
        min_cj=min(costs(j,:));
        max_cj=max(costs(j,:));
        
        dcj=alpha*(max_cj-min_cj);
        
        min_cj=min_cj-dcj;
        max_cj=max_cj+dcj;
        
        gx=linspace(min_cj,max_cj,ngrid-1);
        
        G(j).Lower=[-inf gx];
        G(j).Upper=[gx inf];
        
    end

end