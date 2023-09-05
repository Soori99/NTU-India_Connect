function [FinalBestPath, FinalBestCost] = Traj_MOGWO(startp,endp,startpgreen,startpred,endpgreen,endpred,obs,USVvelocity,option,headingangle,nVar)
    tic
    drawing_flag = 0;
    Lusv = 24.25;
    Pset = [0.3,0.3,0.4];  
    VarSize=[2 nVar];
    TestProblem='F1F2F3';
    fobj = USVobj(TestProblem);
    %% Boundary Setting
    xmax = endp(1,1); %positive x axis boundary
    xmin = startp(1,1); %negative x axis boundary
    smalld = (xmax-xmin)/(nVar+1);
    if startp(1,2) < endp(1,2)
        ymax = endp(1,2)+ Lusv; %positive y axis boundary
        ymin = startp(1,2)-Lusv; %negative y axis boundary
    else
        ymax = startp(1,2)+ Lusv; %positive y axis boundary
        ymin = endp(1,2) - Lusv; %negative y axis boundary
    end

    GreyWolves_num=30;
    Tmax=50;  % Maximum Number of Iterations
    Archive_size=20;   % Repository Size

    alpha=0.1;  % Grid Inflation Parameter
    nGrid=10;   % Number of Grids per each Dimension
    beta=4;   % Leader Selection Pressure Parameter
    gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure
    
    
    %% Grey Wolves Intialisation
    GreyWolves=CreateEmptyParticle(GreyWolves_num);
    noofobs = numel(obs)/2;
    tmp = noofobs+4;      
    
 for i=1:GreyWolves_num
    GreyWolves(i).Velocity=USVvelocity;
    GreyWolves(i).Position=zeros(2,nVar);       
    j=1;      
    while(j<nVar+1)
        check =0;
        t=1;
        if j==1
            GreyWolves(i).Position(1,j) = startp(1,1)+smalld;
            GreyWolves(i).Position(2,j) = unifrnd(startp(1,2)-Lusv/4,startp(1,2)+Lusv/4);
        else     
            GreyWolves(i).Position(1,j)= unifrnd(GreyWolves(i).Position(1,j-1),(GreyWolves(i).Position(1,j-1)+smalld));
            GreyWolves(i).Position(2,j)= unifrnd(ymin,ymax);
        end

        %% Safe Area Constraint
        while (t<tmp+1)
            if t==1
                if sqrt((startpgreen(1,1)-GreyWolves(i).Position(1,j))^2+(startpgreen(1,2)-GreyWolves(i).Position(2,j))^2)<(Lusv+3.25)
                    check = 1;
                end
            elseif t==2
                if sqrt((startpred(1,1)-GreyWolves(i).Position(1,j))^2+(startpred(1,2)-GreyWolves(i).Position(2,j))^2)<(Lusv+3.25)
                    check = 2;
                end
            elseif t==3
                if sqrt((endpgreen(1,1)-GreyWolves(i).Position(1,j))^2+(endpgreen(1,2)-GreyWolves(i).Position(2,j))^2)<(Lusv+3.25)
                    check = 3;
                end
            elseif t==4
                if sqrt((endpred(1,1)-GreyWolves(i).Position(1,j))^2+(endpred(1,2)-GreyWolves(i).Position(2,j))^2)<(Lusv+3.25)
                    check = 4;
                end
            else
                for k=1:noofobs
                    if sqrt((obs(k,1)-GreyWolves(i).Position(1,j))^2+(obs(k,2)-GreyWolves(i).Position(2,j))^2) <(Lusv+2.35)
                        check = 5;
                    end
                end      
            end
            t=t+1;
        end 
        %% Angle Constraint
        thetamax = 60;
        theta = 0;
        theta1 = 0;
        if j==2
            theta = theta + abs(atan2(GreyWolves(i).Position(2,j)-GreyWolves(i).Position(2,j-1),GreyWolves(i).Position(1,j)-GreyWolves(i).Position(1,j-1))-atan2(GreyWolves(i).Position(2,j-1)-startp(1,2),GreyWolves(i).Position(1,j-1)-startp(1,1)))*180/pi;
        elseif j==1
            theta =0;
        else
            theta = theta + abs(atan2(GreyWolves(i).Position(2,j)-GreyWolves(i).Position(2,j-1),GreyWolves(i).Position(1,j)-GreyWolves(i).Position(1,j-1))-atan2(GreyWolves(i).Position(2,j-1)-GreyWolves(i).Position(2,j-2),GreyWolves(i).Position(1,j-1)-GreyWolves(i).Position(1,j-2)))*180/pi;
            if j==nVar
                theta1 = theta1 + abs(atan2(endp(1,2)-GreyWolves(i).Position(2,j),endp(1,1)-GreyWolves(i).Position(1,j))-atan2(GreyWolves(i).Position(2,j)-GreyWolves(i).Position(2,j-1),GreyWolves(i).Position(1,j)-GreyWolves(i).Position(1,j-1)))*180/pi;
            end
        end
        if theta>thetamax
            check = 6;
        end
        if theta<thetamax       
            if j==nVar
                if theta1>thetamax
                    check =6; 
                end
            end
        end
            %% If all constraints are satisified
        if check ==0
            j=j+1;
        elseif check ==6
            j=1;
        end
    end
    GreyWolves(i).Position = sortrows(GreyWolves(i).Position');
    GreyWolves(i).Cost=fobj(startp,endp,GreyWolves(i).Position,nVar,startpgreen,startpred,endpgreen,endpred,obs);
    GreyWolves(i).Best.Position=GreyWolves(i).Position;
    GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
 end
    %% Intialising Archive and Storing Non Dominated Solutions 

    GreyWolves=DetermineDomination(GreyWolves);

    Archive=GetNonDominatedParticles(GreyWolves);

    Archive_costs=GetCosts(Archive);

    G=CreateHypercubes(Archive_costs,nGrid,alpha);

    for i=1:numel(Archive)
        [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
    end

    %% MOGWO Main Loop
    
    for ti=1:Tmax
    a=2-ti*((2)/Tmax);
        constraints = zeros(1,3);
        for i=1:GreyWolves_num
%             Choose the alpha, beta, and delta grey wolves
            Delta=SelectLeader(Archive,beta);
            Beta=SelectLeader(Archive,beta);
            Alpha=SelectLeader(Archive,beta);
            
            c=2.*rand(1, nVar);
            D=abs(c.*Delta.Position' -GreyWolves(i).Position');
            A=2.*a.*rand(1, nVar)-a;
            X1=Delta.Position'-A.*abs(D);

            c=2.*rand(1, nVar);
            D=abs(c.*Beta.Position'-GreyWolves(i).Position');
            A=2.*a.*rand(1, nVar)-a;
            X2=Beta.Position'-A.*abs(D);

            c=2.*rand(1, nVar);
            D=abs(c.*Alpha.Position'-GreyWolves(i).Position');
            A=2.*a.*rand(1, nVar)-a;
            X3=Alpha.Position'-A.*abs(D);

            GreyWolves(i).Position= (X1+X2+X3)./3;
            %% Checking of constraints

%         Boundary Contraint
            for k=1:nVar
                GreyWolves(i).Position(1,k) = min(max(GreyWolves(i).Position(1,k),xmin),xmax);
                GreyWolves(i).Position(2,k) = min(max(GreyWolves(i).Position(2,k),ymin),ymax);
            end

%           Safe Area Constraint
            j=1;
            while(j<nVar+1)
                while (t<tmp+1)
                    if t==1
                        if sqrt((startpgreen(1,1)-GreyWolves(i).Position(1,j))^2+(startpgreen(1,2)-GreyWolves(i).Position(2,j))^2)<3.25
                            constraints(1,1) = 1;
                        end
                    elseif t==2
                        if sqrt((startpred(1,1)-GreyWolves(i).Position(1,j))^2+(startpred(1,2)-GreyWolves(i).Position(2,j))^2)<3.25
                            constraints(1,1) = 1;
                        end
                    elseif t==3
                        if sqrt((endpgreen(1,1)-GreyWolves(i).Position(1,j))^2+(endpgreen(1,2)-GreyWolves(i).Position(2,j))^2)<3.25
                            constraints(1,1) = 1;
                        end
                    elseif t==4
                        if sqrt((endpred(1,1)-GreyWolves(i).Position(1,j))^2+(endpred(1,2)-GreyWolves(i).Position(2,j))^2)<3.25
                            constraints(1,1) = 1;
                        end
                    else
                        for k=1:noofobs
                            if sqrt((obs(k,1)-GreyWolves(i).Position(1,j))^2+(obs(k,2)-GreyWolves(i).Position(2,j))^2) <2.35
                                constraints(1,1) = 1;
                            end
                        end      
                    end
                    t=t+1;
                end
%           Angle Constraint
                thetamax = 60;
                theta = 0;
                theta1 = 0;
                if j==2
                    theta = theta + abs(atan2(GreyWolves(i).Position(2,j)-GreyWolves(i).Position(2,j-1),GreyWolves(i).Position(1,j)-GreyWolves(i).Position(1,j-1))-atan2(GreyWolves(i).Position(2,j-1)-startp(1,2),GreyWolves(i).Position(1,j-1)-startp(1,1)))*180/pi;
                elseif j==1
                    theta =0;
                else
                   theta = theta + abs(atan2(GreyWolves(i).Position(2,j)-GreyWolves(i).Position(2,j-1),GreyWolves(i).Position(1,j)-GreyWolves(i).Position(1,j-1))-atan2(GreyWolves(i).Position(2,j-1)-GreyWolves(i).Position(2,j-2),GreyWolves(i).Position(1,j-1)-GreyWolves(i).Position(1,j-2)))*180/pi;
                   if j==nVar
                       theta1 = theta1 + abs(atan2(endp(1,2)-GreyWolves(i).Position(2,j),endp(1,1)-GreyWolves(i).Position(1,j))-atan2(GreyWolves(i).Position(2,j)-GreyWolves(i).Position(2,j-1),GreyWolves(i).Position(1,j)-GreyWolves(i).Position(1,j-1)))*180/pi;
                   end
                end
                if theta>thetamax
                    constraints(1,2) = 1;
                end
                if theta<thetamax       
                    if j==nVar
                        if theta1>thetamax
                            constraints(1,2) = 1; 
                        end
                    end
                end
                
%         Velocity Constraint
                [velxy, wavecurrent] = CurrentEffects(GreyWolves(i).Position(1,j),GreyWolves(i).Position(2,j),option,ti,headingangle,USVvelocity);
                USVcurrentvelocity = sqrt(velxy(1,1)^2+velxy(1,2)^2);
                Velocitynew = USVcurrentvelocity + wavecurrent;
                if Velocitynew<0
                    constraints(1,3) =1;
                end
            j=j+1;
            end
            GreyWolves(i).Position = GreyWolves(i).Position';

            sumofconstraints=0;
            for s=1:3
                sumofconstraints = sumofconstraints+constraints(1,s);
            end
            
            if sumofconstraints ==0
                GreyWolves(i).Cost=fobj(startp,endp,GreyWolves(i).Position,nVar,startpgreen,startpred,endpgreen,endpred,obs);
                GreyWolves(i).Best.Cost = GreyWolves(i).Cost;
                GreyWolves(i).Best.Position = GreyWolves(i).Position;
            else
                GreyWolves(i).Position = GreyWolves(i).Best.Position;
                GreyWolves(i).Cost = GreyWolves(i).Best.Cost;
            end
        end

        GreyWolves=DetermineDomination(GreyWolves);
        non_dominated_wolves=GetNonDominatedParticles(GreyWolves);

        Archive=[Archive
            non_dominated_wolves];

        Archive=DetermineDomination(Archive);
        Archive=GetNonDominatedParticles(Archive);

        for i=1:numel(Archive)
            [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
        end

        if numel(Archive)>Archive_size
            EXTRA=numel(Archive)-Archive_size;
            Archive=DeleteFromRep(Archive,EXTRA,gamma);

            Archive_costs=GetCosts(Archive);
            G=CreateHypercubes(Archive_costs,nGrid,alpha);    
        end

        disp(['In iteration ' num2str(ti) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
        save results

    %     Results

        costs=GetCosts(GreyWolves);
        Archive_costs=GetCosts(Archive);

        if drawing_flag==1
            hold off
            plot(costs(1,:),costs(2,:),'k.');
            hold on
            plot(Archive_costs(1,:),Archive_costs(2,:),'rd');
            legend('Grey wolves','Non-dominated solutions');
            drawnow
        end
    end

    [m,n] = size(Archive_costs);

    maximum = max(Archive_costs,[],2);
    minimum = min(Archive_costs,[],2);
    NormalisedArchiveCosts = zeros(m,n);

    for i=1:m
        for j=1:n
            if maximum(i,1)==minimum(i,1)
                NormalisedArchiveCosts(i,j)= 0;
            else
            NormalisedArchiveCosts(i,j) = 0+((Archive_costs(i,j)- minimum(i,1))*(1-0))/(maximum(i,1) - minimum(i,1));
            end
        end
    end

    NormalisedCosts = Pset*NormalisedArchiveCosts;
    BestNormalisedSolution = min(NormalisedCosts);

    no =1;
    k=1;
    while(no<n+1)
        if BestNormalisedSolution == NormalisedCosts(1,no)
            solutionnumber(k,1) = no;
            k=k+1;
        end
        no = no+1;
    end

    [a, b] = size(solutionnumber);
    SelectPath = randi([1 a],1,1);
    values = struct2cell(Archive);
    BestPath = values(1,SelectPath);
    BestCost = values(3,SelectPath);
    FinalBestCost = cell2mat(BestCost)
    FinalBestPath = cell2mat(BestPath)
    for j = 1:3
        disp(['Objective #' num2str(j) ':']);
        disp(['      Min = ' num2str(min(Archive_costs(j, :)))]);
        disp(['      Max = ' num2str(max(Archive_costs(j, :)))]);
        disp(['    Range = ' num2str(max(Archive_costs(j, :))-min(Archive_costs(j, :)))]);
        disp(['    St.D. = ' num2str(std(Archive_costs(j, :)))]);
        disp(['     Mean = ' num2str(mean(Archive_costs(j, :)))]);
        disp(' ');
    end
    toc
end
