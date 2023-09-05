clc;
clear;
close all;

%% Problem Definition

nVar=3;         % Number of Decision Variables

VarSize = [1 nVar];   % Decision Variables Matrix Size

VarMin = 0;         % Decision Variables Lower Bound
VarMax = 1;         % Decision Variables Upper Bound

nObj = 2;

%% MOEA/D Settings

MaxIt = 500;  % Maximum Number of Iterations

nPop = 100;    % Population Size (Number of Sub-Problems)

nArchive = 100;

T = max(ceil(0.15*nPop), 2);    % Number of Neighbors
T = min(max(T, 2), 15);

crossover_params.gamma = 0.5;
crossover_params.VarMin = VarMin;
crossover_params.VarMax = VarMax;

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

% Create Initial Population
pop = repmat(empty_individual, nPop, 1);
for i = 1:nPop
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    pop(i).Cost = (zdt3(pop(i).Position))';
    z = min(z, pop(i).Cost);
end

for i = 1:nPop
    pop(i).g = DecomposedCost(pop(i), z, sp(i).lambda);
end

% Determine Population Domination Status
pop = DetermineDomination2(pop);

% Initialize Estimated Pareto Front
EP = pop(~[pop.IsDominated]);

%% Main Loop

for it = 1:MaxIt
    for i = 1:nPop
        
        % Reproduction (Crossover)
        K = randsample(T, 2);
        
        j1 = sp(i).Neighbors(K(1));
        p1 = pop(j1);
        
        j2 = sp(i).Neighbors(K(2));
        p2 = pop(j2);
        
        y = empty_individual;
        y.Position = Crossover(p1.Position, p2.Position, crossover_params);
        
        y.Cost = (zdt3(y.Position))';
        
        z = min(z, y.Cost);
        
        for j = sp(i).Neighbors
            y.g = DecomposedCost(y, z, sp(j).lambda);
            if y.g <= pop(j).g
                pop(j) = y;
            end
        end
        
    end
    
    % Determine Population Domination Status
	pop = DetermineDomination2(pop);
    
    ndpop = pop(~[pop.IsDominated]);
    
    EP = [EP
        ndpop]; %#ok
    
    EP = DetermineDomination2(EP);
    EP = EP(~[EP.IsDominated]);
    
    if numel(EP)>nArchive
        Extra = numel(EP)-nArchive;
        ToBeDeleted = randsample(numel(EP), Extra);
        EP(ToBeDeleted) = [];
    end
        
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Pareto Solutions = ' num2str(numel(EP))]);
    
end

%% Reults

EPC = [EP.Cost];
for j = 1:nObj
    
    disp(['Objective #' num2str(j) ':']);
    disp(['      Min = ' num2str(min(EPC(:,j)))]);
    disp(['      Max = ' num2str(max(EPC(:,j)))]);
    disp(['    Range = ' num2str(max(EPC(:,j))-min(EPC(:,j)))]);
    disp(['    St.D. = ' num2str(std(EPC(:,j)))]);
    disp(['     Mean = ' num2str(mean(EPC(:,j)))]);
    disp(' ');
    
end
True_Pareto=load('ZDT3.txt');
%% Plot data
figure(2);
plot(EPC(1,:),EPC(2,:),'o','LineWidth',2,...
    'MarkerEdgeColor','r','MarkerSize',2);
hold on
plot(True_Pareto(:,1),True_Pareto(:,2),'k'); 
title('Optimal Solution Pareto Set using MOEA/D');
legend('MOEA/D');
xlabel('F_1');
ylabel('F_2');
%%  Metric Value
M_IGD=IGD(EPC',True_Pareto);
M_GD=GD(EPC',True_Pareto);
M_HV=HV(EPC',True_Pareto);
M_Spacing=Spacing(EPC',True_Pareto);
M_Spread=Spread(EPC',True_Pareto);
M_DeltaP=DeltaP(EPC',True_Pareto);
display(['The IGD Metric obtained by MOEA/D is     : ', num2str(M_IGD)]);
display(['The GD Metric obtained by MOEA/D is      : ', num2str(M_GD)]);
display(['The HV Metric obtained by MOEA/D is      : ', num2str(M_HV)]);
display(['The Spacing Metric obtained by MOEA/D is : ', num2str(M_Spacing)]);
display(['The Spread Metric obtained by MOEA/D is  : ', num2str(M_Spread)]);
display(['The DeltaP Metric obtained by MOEA/D is  : ', num2str(M_DeltaP)]);
