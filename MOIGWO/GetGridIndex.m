function [Index SubIndex]=GetGridIndex(particle,G)

    c=particle.Cost;
    
    nobj=numel(c);
    ngrid=numel(G(1).Upper);
    
    str=['sub2ind(' mat2str(ones(1,nobj)*ngrid)];

    SubIndex=zeros(1,nobj);
    for j=1:nobj
        
        U=G(j).Upper;
        
        i=find(c(j)<U,1,'first');
        
        SubIndex(j)=i;
        
        str=[str ',' num2str(i)];
    end
    
    str=[str ');'];
    
    Index=eval(str);
    
end