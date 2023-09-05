function [FinalBestPath, FinalBestCost] = Traj_MOEAD(startp,endp,startpgreen,startpred,endpgreen,endpred,obs,USVvelocity,option,headingangle,nVar)
    tic
    drawing_flag = 0;
    Lusv = 24.25;
    Pset = [0.3,0.3,0.4];  
    %% Problem Definition
    
    TestProblem='F1F2F3';
    fobj = USVobj(TestProblem);
    VarSize = [2 nVar];   % Decision Variables Matrix Size

    nObj = 3;
    
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

    %% MOEA/D Settings

    Tmax = 50;  % Maximum Number of Iterations

    nPop = 30;    % Population Size (Number of Sub-Problems)

    nArchive = 30;

    T = max(ceil(0.15*nPop), 2);    % Number of Neighbors
    T = min(max(T, 2), 15);
    
    VarMin = [xmin, ymin];         % Decision Variables Lower Bound
    VarMax = [xmax, ymax];         % Decision Variables Upper Bound

    crossover_params.gamma = 0.5;
    %% Initialization

    % Create Sub-problems
    sp = CreateSubProblems(nObj, nPop, T);

    % Empty Individual
    empty_individual.Position = [];
    empty_individual.Cost = [];
    empty_individual.g = [];
    empty_individual.IsDominated = [];

    % Initialize Goal Point
    z = zeros(nObj, 1);
      
    noofobs = numel(obs)/2;
    tmp = noofobs+4;

    % Create Initial Population
    pop = repmat(empty_individual, nPop, 1);
    
    for i = 1:nPop
        pop(i).Position=zeros(2,nVar);       
        j=1;   
        while(j<nVar+1)
            check =0;
            t=1;
            if j==1
                pop(i).Position(1,j) = startp(1,1)+smalld;
                pop(i).Position(2,j) = unifrnd(startp(1,2)-Lusv/4,startp(1,2)+Lusv/4);
            else     
                pop(i).Position(1,j)= unifrnd(pop(i).Position(1,j-1),(pop(i).Position(1,j-1)+smalld));
                pop(i).Position(2,j)= unifrnd(ymin,ymax);
            end
    %% Safe Area Constraint
            while (t<tmp+1)
                if t==1
                    if sqrt((startpgreen(1,1)-pop(i).Position(1,j))^2+(startpgreen(1,2)-pop(i).Position(2,j))^2)<(Lusv+3.25)
                        check = 1;
                    end
                elseif t==2
                    if sqrt((startpred(1,1)-pop(i).Position(1,j))^2+(startpred(1,2)-pop(i).Position(2,j))^2)<(Lusv+3.25)
                        check = 2;
                    end
                elseif t==3
                    if sqrt((endpgreen(1,1)-pop(i).Position(1,j))^2+(endpgreen(1,2)-pop(i).Position(2,j))^2)<(Lusv+3.25)
                        check = 3;
                    end
                elseif t==4
                    if sqrt((endpred(1,1)-pop(i).Position(1,j))^2+(endpred(1,2)-pop(i).Position(2,j))^2)<(Lusv+3.25)
                        check = 4;
                    end
                else
                    for k=1:noofobs
                        if sqrt((obs(k,1)-pop(i).Position(1,j))^2+(obs(k,2)-pop(i).Position(2,j))^2) <(Lusv+2.35)
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
                theta = theta + abs(atan2(pop(i).Position(2,j)-pop(i).Position(2,j-1),pop(i).Position(1,j)-pop(i).Position(1,j-1))-atan2(pop(i).Position(2,j-1)-startp(1,2),pop(i).Position(1,j-1)-startp(1,1)))*180/pi;
            elseif j==1
                theta =0;
            else
                theta = theta + abs(atan2(pop(i).Position(2,j)-pop(i).Position(2,j-1),pop(i).Position(1,j)-pop(i).Position(1,j-1))-atan2(pop(i).Position(2,j-1)-pop(i).Position(2,j-2),pop(i).Position(1,j-1)-pop(i).Position(1,j-2)))*180/pi;
                if j==nVar
                    theta1 = theta1 + abs(atan2(endp(1,2)-pop(i).Position(2,j),endp(1,1)-pop(i).Position(1,j))-atan2(pop(i).Position(2,j)-pop(i).Position(2,j-1),pop(i).Position(1,j)-pop(i).Position(1,j-1)))*180/pi;
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
        pop(i).Position = sortrows(pop(i).Position');
        pop(i).Cost=fobj(startp,endp,pop(i).Position,nVar,startpgreen,startpred,endpgreen,endpred,obs)';
        z = min(z, pop(i).Cost);
    end

    for i = 1:nPop
        pop(i).g = DecomposedCost(pop(i), z, sp(i).lambda);
    end

    % Determine Population Domination Status
    pop = DetermineDomination2(pop);

    % Initialize Estimated Pareto Front
    RP = pop(~[pop.IsDominated]);

    %% MOEAD Main Loop

    for ti = 1:Tmax
        constraints = zeros(1,3);
        for i = 1:nPop

            % Reproduction (Crossover)
            K = randsample(T, 2);

            j1 = sp(i).Neighbors(K(1));
            p1 = pop(j1);

            j2 = sp(i).Neighbors(K(2));
            p2 = pop(j2);

            y = empty_individual;
            
            for t=1:2
                    crossover_params.VarMin = VarMin(t);
                    crossover_params.VarMax = VarMax(t);
                    for j =1:nVar
                        y.Position(j,t) = Crossover(p1.Position(j,t), p2.Position(j,t), crossover_params);
                    end
            end
            
%           Safe Area Constraint
            j=1;
            t=1;
            while(j<nVar+1)
                while (t<tmp+1)
                    if t==1
                        if sqrt((startpgreen(1,1)-y.Position(j,1))^2+(startpgreen(1,2)-y.Position(j,2))^2)<3.25
                            constraints(1,1) = 1;
                        end
                    elseif t==2
                        if sqrt((startpred(1,1)-y.Position(j,1))^2+(startpred(1,2)-y.Position(j,2))^2)<3.25
                            constraints(1,1) = 1;
                        end
                    elseif t==3
                        if sqrt((endpgreen(1,1)-y.Position(j,1))^2+(endpgreen(1,2)-y.Position(j,2))^2)<3.25
                            constraints(1,1) = 1;
                        end
                    elseif t==4
                        if sqrt((endpred(1,1)-y.Position(j,1))^2+(endpred(1,2)-y.Position(j,2))^2)<3.25
                            constraints(1,1) = 1;
                        end
                    else
                        for k=1:noofobs
                            if sqrt((obs(k,1)-y.Position(j,1))^2+(obs(k,2)-y.Position(j,2))^2) <2.35
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
                    theta = theta + abs(atan2(y.Position(j,2)-y.Position(j-1,2),y.Position(j,1)-y.Position(j-1,1))-atan2(y.Position(j-1,2)-startp(1,2),y.Position(j-1,1)-startp(1,1)))*180/pi;
                elseif j==1
                    theta =0;
                else
                   theta = theta + abs(atan2(y.Position(j,2)-y.Position(j-1,2),y.Position(j,1)-y.Position(j-1,1))-atan2(y.Position(j-1,2)-y.Position(j-2,2),y.Position(j-1,1)-y.Position(j-2,1)))*180/pi;
                   if j==nVar
                       theta1 = theta1 + abs(atan2(endp(1,2)-y.Position(j,2),endp(1,1)-y.Position(j,1))-atan2(y.Position(j,2)-y.Position(j-1,2),y.Position(j,1)-y.Position(j-1,1)))*180/pi;
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
                [velxy, wavecurrent] = CurrentEffects(y.Position(j,1),y.Position(j,2),option,ti,headingangle,USVvelocity);
                USVcurrentvelocity = sqrt(velxy(1,1)^2+velxy(1,2)^2);
                Velocitynew = USVcurrentvelocity + wavecurrent;
                if Velocitynew<0
                    constraints(1,3) =1;
                end
            j=j+1;
            end
            
            sumofconstraints=0;
            for s=1:3
                sumofconstraints = sumofconstraints+constraints(1,s);
            end

            y.Cost=fobj(startp,endp,y.Position,nVar,startpgreen,startpred,endpgreen,endpred,obs)';

            if sumofconstraints==0
                z = min(z, y.Cost);
            end

        end
        for j = 1:sp(i).Neighbors
            y.g = DecomposedCost(y, z, sp(j).lambda);
            if y.g <= pop(j).g
                pop(j) = y;
            end
       end

        % Determine Population Domination Status
        pop = DetermineDomination2(pop);

        ndpop = pop(~[pop.IsDominated]);

        RP = [RP
            ndpop]; %#ok

        RP = DetermineDomination2(RP);
        RP = RP(~[RP.IsDominated]);

        if numel(RP)>nArchive
            Extra = numel(RP)-nArchive;
            ToBeDeleted = randsample(numel(RP), Extra);
            RP(ToBeDeleted) = [];
        end
        
        if drawing_flag ==1
            % Plot RP
            figure(1);
            PlotCosts(RP);
            pause(0.01);
        end


        % Display Iteration Information
        disp(['Iteration ' num2str(ti) ': Number of Pareto Solutions = ' num2str(numel(RP))]);

    end

    %% Reults
 rep_costs = [RP.Cost];
 [m,n] = size(rep_costs);


    maximum = max(rep_costs,[],2);
    minimum = min(rep_costs,[],2);
    NormalisedRepCosts = zeros(m,n);

    for i=1:m
        for j=1:n
            if maximum(i,1)==minimum(i,1)
                NormalisedRepCosts(i,j)= 0;
            else
            NormalisedRepCosts(i,j) = 0+((rep_costs(i,j)- minimum(i,1))*(1-0))/(maximum(i,1) - minimum(i,1));
            end
        end
    end

    NormalisedCosts = Pset*NormalisedRepCosts;
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
    values = struct2cell(RP);
    BestPath = values(1,SelectPath);
    BestCost = values(2,SelectPath);
    FinalBestCost = (cell2mat(BestCost))'
    FinalBestPath = cell2mat(BestPath)
    for j = 1:3
        disp(['Objective #' num2str(j) ':']);
        disp(['      Min = ' num2str(min(rep_costs(j, :)))]);
        disp(['      Max = ' num2str(max(rep_costs(j, :)))]);
        disp(['    Range = ' num2str(max(rep_costs(j, :))-min(rep_costs(j, :)))]);
        disp(['    St.D. = ' num2str(std(rep_costs(j, :)))]);
        disp(['     Mean = ' num2str(mean(rep_costs(j, :)))]);
        disp(' ');
    end
    toc
end