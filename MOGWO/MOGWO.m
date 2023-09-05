clear
clc

nVar=3;
drawing_flag = 1;
nObj = 2;
Pset = [0.5,0.5];

lb = ones(1, nVar).*0; 
ub = ones(1, nVar).*1;

VarSize=[1 nVar];
GreyWolves_num=100;
MaxIt=500;  % Maximum Number of Iterations
Archive_size=100;   % Repository Size

alpha=0.1;  % Grid Inflation Parameter
nGrid=10;   % Number of Grids per each Dimension
beta=4;   % Leader Selection Pressure Parameter
gamma=2;    % Extra (to be deleted) Repository Member Selection Pressure

% Initialization

GreyWolves=CreateEmptyParticle(GreyWolves_num);

for i=1:GreyWolves_num
    GreyWolves(i).Velocity=0;
    GreyWolves(i).Position=zeros(1,nVar);
    GreyWolves(i).Position(1,:)=lb+(ub-lb).*rand(1,nVar);
    GreyWolves(i).Cost = zdt3(GreyWolves(i).Position);
    GreyWolves(i).Best.Position=GreyWolves(i).Position;
    GreyWolves(i).Best.Cost=GreyWolves(i).Cost;
end

GreyWolves=DetermineDomination(GreyWolves);

Archive=GetNonDominatedParticles(GreyWolves);

Archive_costs=GetCosts(Archive);
G=CreateHypercubes(Archive_costs,nGrid,alpha);

for i=1:numel(Archive)
    [Archive(i).GridIndex, Archive(i).GridSubIndex]=GetGridIndex(Archive(i),G);
end

% MOGWO main loop

for it=1:MaxIt
    a=2-it*((2)/MaxIt);
    for i=1:GreyWolves_num
        
        clear rep2
        clear rep3
        
        % Choose the alpha, beta, and delta grey wolves
        Delta=SelectLeader(Archive,beta);
        Beta=SelectLeader(Archive,beta);
        Alpha=SelectLeader(Archive,beta);
        
        % If there are less than three solutions in the least crowded
        % hypercube, the second least crowded hypercube is also found
        % to choose other leaders from.
        if size(Archive,1)>1
            counter=0;
            for newi=1:size(Archive,1)
                if sum(Delta.Position~=Archive(newi).Position)~=0
                    counter=counter+1;
                    rep2(counter,1)=Archive(newi);
                end
            end
            Beta=SelectLeader(rep2,beta);
        end
        
        % This scenario is the same if the second least crowded hypercube
        % has one solution, so the delta leader should be chosen from the
        % third least crowded hypercube.
        if size(Archive,1)>2
            counter=0;
            for newi=1:size(rep2,1)
                if sum(Beta.Position~=rep2(newi).Position)~=0
                    counter=counter+1;
                    rep3(counter,1)=rep2(newi);
                end
            end
            Alpha=SelectLeader(rep3,beta);
        end
        
        c=2.*rand(1, nVar);
        D=abs(c.*Delta.Position-GreyWolves(i).Position);
        A=2.*a.*rand(1, nVar)-a;
        X1=Delta.Position-A.*abs(D);
        
        c=2.*rand(1, nVar);
        D=abs(c.*Beta.Position-GreyWolves(i).Position);
        A=2.*a.*rand(1, nVar)-a;
        X2=Beta.Position-A.*abs(D);
        
        c=2.*rand(1, nVar);
        D=abs(c.*Alpha.Position-GreyWolves(i).Position);
        A=2.*a.*rand(1, nVar)-a;
        X3=Alpha.Position-A.*abs(D);
        
        GreyWolves(i).Position=(X1+X2+X3)./3;
        
        % Boundary checking
        GreyWolves(i).Position=min(max(GreyWolves(i).Position,lb),ub);
        GreyWolves(i).Cost = zdt3(GreyWolves(i).Position);
        
        costsum =  GreyWolves(i).Cost * Pset';
        bestcostsum =  GreyWolves(i).Best.Cost * Pset';
        if costsum<bestcostsum
            GreyWolves(i).Best.Position = GreyWolves(i).Position;
            GreyWolves(i).Best.Cost = GreyWolves(i).Cost;            
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
    
    disp(['In iteration ' num2str(it) ': Number of solutions in the archive = ' num2str(numel(Archive))]);
    save results
    
    % Results
    
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

True_Pareto=load('ZDT3.txt');
%% Plot data
figure(2);
plot(Archive_costs(1,:),Archive_costs(2,:),'o','LineWidth',2,...
    'MarkerEdgeColor','r','MarkerSize',2);
hold on
plot(True_Pareto(:,1),True_Pareto(:,2),'k'); 
title('Optimal Solution Pareto Set using MOGWO');
legend('MOGWO');
xlabel('F_1');
ylabel('F_2');

for j = 1:nObj
    disp(['Objective #' num2str(j) ':']);
    disp(['      Min = ' num2str(min(Archive_costs(j, :)))]);
    disp(['      Max = ' num2str(max(Archive_costs(j, :)))]);
    disp(['    Range = ' num2str(max(Archive_costs(j, :))-min(Archive_costs(j, :)))]);
    disp(['    St.D. = ' num2str(std(Archive_costs(j, :)))]);
    disp(['     Mean = ' num2str(mean(Archive_costs(j, :)))]);
    disp(' ');
end
%%  Metric Value
M_IGD=IGD(Archive_costs',True_Pareto);
M_GD=GD(Archive_costs',True_Pareto);
M_HV=HV(Archive_costs',True_Pareto);
M_Spacing=Spacing(Archive_costs',True_Pareto);
M_Spread=Spread(Archive_costs',True_Pareto);
M_DeltaP=DeltaP(Archive_costs',True_Pareto);
display(['The IGD Metric obtained by MOGWO is     : ', num2str(M_IGD)]);
display(['The GD Metric obtained by MOGWO is      : ', num2str(M_GD)]);
display(['The HV Metric obtained by MOGWO is      : ', num2str(M_HV)]);
display(['The Spacing Metric obtained by MOGWO is : ', num2str(M_Spacing)]);
display(['The Spread Metric obtained by MOGWO is  : ', num2str(M_Spread)]);
display(['The DeltaP Metric obtained by MOGWO is  : ', num2str(M_DeltaP)]);

