function rep=DeleteFromRep(rep,EXTRA,gamma)

    if nargin<3
        gamma=1;
    end

    for k=1:EXTRA
        [occ_cell_index occ_cell_member_count]=GetOccupiedCells(rep);

        p=occ_cell_member_count.^gamma;
        p=p/sum(p);

        selected_cell_index=occ_cell_index(RouletteWheelSelection(p));

        GridIndices=[rep.GridIndex];

        selected_cell_members=find(GridIndices==selected_cell_index);

        n=numel(selected_cell_members);

        selected_memebr_index=randi([1 n]);

        j=selected_cell_members(selected_memebr_index);
        
        rep=[rep(1:j-1); rep(j+1:end)];
    end
    
end