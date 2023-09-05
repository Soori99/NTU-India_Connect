function fobj = USVobj(name)
    switch name
        case 'F1'
            fobj = @F1;
        case 'F2'
            fobj = @F2; 
        case 'F3'
            fobj = @F3;  
        case 'F1F2'
            fobj = @F1F2;
        case 'F2F3'
            fobj = @F2F3; 
        case 'F1F3'
            fobj = @F1F3;
        case 'F1F2F3'
            fobj = @F1F2F3;      
        otherwise
            fobj = @F1;
    end
end
%% Objective 1 - Shortest path length
function minlength = F1(sp,ep,points,nVar,entgreen,entred,exitgreen,exitred,finalobsselection)
    minlength =0;
    minlength = minlength + sqrt((sp(1,1)-points(1,1))^2+(sp(1,2)-points(1,2))^2);
    for i=2:nVar
    minlength = minlength + sqrt((points(i,1)-points(i-1,1))^2+(points(i,2)-points(i-1,2))^2);     
    end
    minlength = minlength + sqrt((ep(1,1)-points(nVar,1))^2+(ep(1,2)-points(nVar,2))^2);
end 
%% Objective 2 - Path Smoothness 
function mintheta = F2(sp,ep,points,nVar,entgreen,entred,exitgreen,exitred,finalobsselection)
    mintheta =0;
    mintheta = mintheta + abs(atan2(points(2,2)-points(1,2),points(2,1)-points(1,1))-atan2(points(1,2)-sp(1,2),points(1,1)-sp(1,1)))*180/pi;
    for i=2:nVar-1
    mintheta = mintheta + abs(atan2(points(i+1,2)-points(i,2),points(i+1,1)-points(i,1))-atan2(points(i,2)-points(i-1,2),points(i,1)-points(i-1,1)))*180/pi;        
    end
    mintheta = mintheta + abs(atan2(ep(1,2)-points(nVar,2),ep(1,1)-points(nVar,1))-atan2(points(nVar,2)-points(nVar-1,2),points(nVar,1)-points(nVar-1,1)))*180/pi;        
end
%% Objective 3 - Path Safety
function minsafety = F3(sp,ep,points,nVar,entgreen,entred,exitgreen,exitred,finalobsselection)
    minsafety = 0;
    noofobs = numel(finalobsselection)/2;
    tmp = noofobs+4;
    dreal = zeros(nVar,tmp);
    dsafemin = ones(nVar,tmp)*3.25;
    dsafemax = ones(nVar,tmp)*27.5;
    for i=1:nVar
        dsafemax(nVar,tmp) = dsafemax(nvar,tmp) - 0.9;
    end
    normaliseddist = zeros(nVar,1);
    for i=1:nVar
        j=1;
        while (j<6)
            if j==1
                dreal(i,j) = dreal(i,j) + sqrt((entgreen(1,1)-points(i,1))^2+(entgreen(1,2)-points(i,2))^2);
            elseif j==2
                dreal(i,j) = dreal(i,j) + sqrt((entred(1,1)-points(i,1))^2+(entred(1,2)-points(i,2))^2);
            elseif j==3
                dreal(i,j) = dreal(i,j) + sqrt((exitgreen(1,1)-points(i,1))^2+(exitgreen(1,2)-points(i,2))^2);
            elseif j==4
                dreal(i,j) = dreal(i,j) + sqrt((exitred(1,1)-points(i,1))^2+(exitred(1,2)-points(i,2))^2);
            else
                for k=1:noofobs
                    dreal(i,j) = dreal(i,j) + sqrt((finalobsselection(k,1)-points(i,1))^2+(finalobsselection(k,2)-points(i,2))^2);
                    dsafemin(i,j) = 2.35;
                end      
            end
            j=j+1;
        end
    end
    for i=1:nVar
        for j=1:5
            if dreal(i,j)<=dsafemin(i,j)
                normaliseddist(i) = normaliseddist(i) + 1;
            elseif dreal(i,j)>=dsafemax(i,j)
                normaliseddist(i) = normaliseddist(i) + 0;
            else
                normaliseddist(i) = normaliseddist(i) + (dreal(i,j)-dsafemin(i,j))/(dsafemax(i,j)-dsafemin(i,j));
            end
        end
    end
    for i=1:tmp
        minsafety = minsafety + normaliseddist(i);
    end   
end
%% Objective 1 and 2
function minlengththeta = F1F2(sp,ep,points,nVar,entgreen,entred,exitgreen,exitred,finalobsselection)
    minlengththeta = zeros(1,2);

    minlength =0;
    minlength = minlength + sqrt((sp(1,1)-points(1,1))^2+(sp(1,2)-points(1,2))^2);
    for i=2:nVar
    minlength = minlength + sqrt((points(i,1)-points(i-1,1))^2+(points(i,2)-points(i-1,2))^2);     
    end
    minlength = minlength + sqrt((ep(1,1)-points(nVar,1))^2+(ep(1,2)-points(nVar,2))^2);
    
    mintheta =0;
    mintheta = mintheta + abs(atan2(points(2,2)-points(1,2),points(2,1)-points(1,1))-atan2(points(1,2)-sp(1,2),points(1,1)-sp(1,1)))*180/pi;
    for i=2:nVar-1
    mintheta = mintheta + abs(atan2(points(i+1,2)-points(i,2),points(i+1,1)-points(i,1))-atan2(points(i,2)-points(i-1,2),points(i,1)-points(i-1,1)))*180/pi;        
    end
    mintheta = mintheta + abs(atan2(ep(1,2)-points(nVar,2),ep(1,1)-points(nVar,1))-atan2(points(nVar,2)-points(nVar-1,2),points(nVar,1)-points(nVar-1,1)))*180/pi;
    
    minlengththeta(1,1) = minlength;
    minlengththeta(1,2) = mintheta;
end
%% Objective 2 and 3
function minthetasafety = F2F3(sp,ep,points,nVar,entgreen,entred,exitgreen,exitred,finalobsselection)
    minthetasafety = zeros(1,2);
    
    mintheta =0;
    mintheta = mintheta + abs(atan2(points(2,2)-points(1,2),points(2,1)-points(1,1))-atan2(points(1,2)-sp(1,2),points(1,1)-sp(1,1)))*180/pi;
    for i=2:nVar-1
    mintheta = mintheta + abs(atan2(points(i+1,2)-points(i,2),points(i+1,1)-points(i,1))-atan2(points(i,2)-points(i-1,2),points(i,1)-points(i-1,1)))*180/pi;  
    end
    mintheta = mintheta + abs(atan2(ep(1,2)-points(nVar,2),ep(1,1)-points(nVar,1))-atan2(points(nVar,2)-points(nVar-1,2),points(nVar,1)-points(nVar-1,1)))*180/pi;      
    
    minsafety = 0;
    noofobs = numel(finalobsselection)/2;
    tmp = 5;
    dreal = zeros(nVar,tmp);
    dsafemin = ones(nVar,tmp)*3.25;
    dsafemax = ones(nVar,tmp)*27.5;
    for i=1:nVar
        dsafemax(i,tmp) = 26.6;
    end
    normaliseddist = zeros(nVar,1);
    for i=1:nVar
        j=1;
        while (j<6)
            if j==1
                dreal(i,j) = dreal(i,j) + sqrt((entgreen(1,1)-points(i,1))^2+(entgreen(1,2)-points(i,2))^2);
            elseif j==2
                dreal(i,j) = dreal(i,j) + sqrt((entred(1,1)-points(i,1))^2+(entred(1,2)-points(i,2))^2);
            elseif j==3
                dreal(i,j) = dreal(i,j) + sqrt((exitgreen(1,1)-points(i,1))^2+(exitgreen(1,2)-points(i,2))^2);
            elseif j==4
                dreal(i,j) = dreal(i,j) + sqrt((exitred(1,1)-points(i,1))^2+(exitred(1,2)-points(i,2))^2);
            else
                for k=1:noofobs
                    dreal(i,j) = dreal(i,j) + sqrt((finalobsselection(k,1)-points(i,1))^2+(finalobsselection(k,2)-points(i,2))^2);
                    dsafemin(i,j) = 2.35;
                end      
            end
            j=j+1;
        end
    end
    for i=1:nVar
        for j=1:5
            if dreal(i,j)<=dsafemin(i,j)
                normaliseddist(i) = normaliseddist(i) + 1;
            elseif dreal(i,j)>=dsafemax(i,j)
                normaliseddist(i) = normaliseddist(i) + 0;
            else
                normaliseddist(i) = normaliseddist(i) + (dreal(i,j)-dsafemin(i,j))/(dsafemax(i,j)-dsafemin(i,j));
            end
        end
    end
    for i=1:nVar
        minsafety = minsafety + normaliseddist(i);
    end 
    
    minthetasafety(1,1) = mintheta;
    minthetasafety(1,2) = minsafety;
end
%% Objective 1 and 3
function minlengthsafety = F1F3(sp,ep,points,nVar,entgreen,entred,exitgreen,exitred,finalobsselection)
    minlengthsafety = zeros(1,2);
    
    minlength =0;
    minlength = minlength + sqrt((sp(1,1)-points(1,1))^2+(sp(1,2)-points(1,2))^2);
    for i=2:nVar
    minlength = minlength + sqrt((points(i,1)-points(i-1,1))^2+(points(i,2)-points(i-1,2))^2);        
    end
    minlength = minlength + sqrt((ep(1,1)-points(nVar,1))^2+(ep(1,2)-points(nVar,2))^2);
    
    minsafety = 0;
    noofobs = numel(finalobsselection)/2;
    tmp = 5;
    dreal = zeros(nVar,tmp);
    dsafemin = ones(nVar,tmp)*3.25;
    dsafemax = ones(nVar,tmp)*27.5;
    for i=1:nVar
        dsafemax(i,tmp) = 26.6;
    end
    normaliseddist = zeros(nVar,1);
    for i=1:nVar
        j=1;
        while (j<6)
            if j==1
                dreal(i,j) = dreal(i,j) + sqrt((entgreen(1,1)-points(i,1))^2+(entgreen(1,2)-points(i,2))^2);
            elseif j==2
                dreal(i,j) = dreal(i,j) + sqrt((entred(1,1)-points(i,1))^2+(entred(1,2)-points(i,2))^2);
            elseif j==3
                dreal(i,j) = dreal(i,j) + sqrt((exitgreen(1,1)-points(i,1))^2+(exitgreen(1,2)-points(i,2))^2);
            elseif j==4
                dreal(i,j) = dreal(i,j) + sqrt((exitred(1,1)-points(i,1))^2+(exitred(1,2)-points(i,2))^2);
            else
                for k=1:noofobs
                    dreal(i,j) = dreal(i,j) + sqrt((finalobsselection(k,1)-points(i,1))^2+(finalobsselection(k,2)-points(i,2))^2);
                    dsafemin(i,j) = 2.35;
                end      
            end
            j=j+1;
        end
    end
    for i=1:nVar
        for j=1:5
            if dreal(i,j)<=dsafemin(i,j)
                normaliseddist(i) = normaliseddist(i) + 1;
            elseif dreal(i,j)>=dsafemax(i,j)
                normaliseddist(i) = normaliseddist(i) + 0;
            else
                normaliseddist(i) = normaliseddist(i) + (dreal(i,j)-dsafemin(i,j))/(dsafemax(i,j)-dsafemin(i,j));
            end
        end
    end
    for i=1:nVar
        minsafety = minsafety + normaliseddist(i);
    end
    
    minlengthsafety(1,1) = minlength;
    minlengthsafety(1,2) = minsafety;    
end

%% Objective 1,2 and 3
function minlengththetasafety = F1F2F3(sp,ep,points,nVar,entgreen,entred,exitgreen,exitred,finalobsselection)
    minlengththetasafety = zeros(1,3);

    minlength =0;
    minlength = minlength + sqrt((sp(1,1)-points(1,1))^2+(sp(1,2)-points(1,2))^2);
    for i=2:nVar
    minlength = minlength + sqrt((points(i,1)-points(i-1,1))^2+(points(i,2)-points(i-1,2))^2);        
    end
    minlength = minlength + sqrt((ep(1,1)-points(nVar,1))^2+(ep(1,2)-points(nVar,2))^2);
    
    mintheta =0;
    mintheta = mintheta + abs(atan2(points(2,2)-points(1,2),points(2,1)-points(1,1))-atan2(points(1,2)-sp(1,2),points(1,1)-sp(1,1)))*180/pi;
    for i=2:nVar-1
    mintheta = mintheta + abs(atan2(points(i+1,2)-points(i,2),points(i+1,1)-points(i,1))-atan2(points(i,2)-points(i-1,2),points(i,1)-points(i-1,1)))*180/pi;  
    end
    mintheta = mintheta + abs(atan2(ep(1,2)-points(nVar,2),ep(1,1)-points(nVar,1))-atan2(points(nVar,2)-points(nVar-1,2),points(nVar,1)-points(nVar-1,1)))*180/pi;      
    
    minsafety = 0;
    noofobs = numel(finalobsselection)/2;
    tmp = 5;
    dreal = zeros(nVar,tmp);
    dsafemin = ones(nVar,tmp)*3.25;
    dsafemax = ones(nVar,tmp)*27.5;
    for i=1:nVar
        dsafemax(i,tmp) = 26.6;
    end
    normaliseddist = zeros(nVar,1);
    for i=1:nVar
        j=1;
        while (j<6)
            if j==1
                dreal(i,j) = dreal(i,j) + sqrt((entgreen(1,1)-points(i,1))^2+(entgreen(1,2)-points(i,2))^2);
            elseif j==2
                dreal(i,j) = dreal(i,j) + sqrt((entred(1,1)-points(i,1))^2+(entred(1,2)-points(i,2))^2);
            elseif j==3
                dreal(i,j) = dreal(i,j) + sqrt((exitgreen(1,1)-points(i,1))^2+(exitgreen(1,2)-points(i,2))^2);
            elseif j==4
                dreal(i,j) = dreal(i,j) + sqrt((exitred(1,1)-points(i,1))^2+(exitred(1,2)-points(i,2))^2);
            else
                for k=1:noofobs
                    dreal(i,j) = dreal(i,j) + sqrt((finalobsselection(k,1)-points(i,1))^2+(finalobsselection(k,2)-points(i,2))^2);
                    dsafemin(i,j) = 2.35;
                end      
            end
            j=j+1;
        end
    end
    for i=1:nVar
        for j=1:5
            if dreal(i,j)<=dsafemin(i,j)
                normaliseddist(i) = normaliseddist(i) + 1;
            elseif dreal(i,j)>=dsafemax(i,j)
                normaliseddist(i) = normaliseddist(i) + 0;
            else
                normaliseddist(i) = normaliseddist(i) + (dreal(i,j)-dsafemin(i,j))/(dsafemax(i,j)-dsafemin(i,j));
            end
        end
    end
    for i=1:nVar
        minsafety = minsafety + normaliseddist(i);
    end   
    
    minlengththetasafety(1,1) = minlength;
    minlengththetasafety(1,2) = mintheta;
    minlengththetasafety(1,3) = minsafety; 
end


