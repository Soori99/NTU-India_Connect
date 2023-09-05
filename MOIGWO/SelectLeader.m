function rep_h=SelectLeader(rep,beta)
    if nargin<2
        beta=1;
    end

    [occ_cell_index, occ_cell_member_count]=GetOccupiedCells(rep);
    
    p=occ_cell_member_count.^(-beta);
    p=p/sum(p);
    
    selected_cell_index=occ_cell_index(RouletteWheelSelection(p));
    
    GridIndices=[rep.GridIndex];
    
    selected_cell_members=find(GridIndices==selected_cell_index);
    
    n=numel(selected_cell_members);
    
    selected_member_index=randi([1 n]);
    
    h=selected_cell_members(selected_member_index);
    
    rep_h=rep(h);
end