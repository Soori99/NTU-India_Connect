clear
clc
drawing_map =1;
%% Number of Gates and Boundary Definitiion
mingates=2; %mininum number of gates 
maxgates=5; % maximum number of gates
nog = randi([maxgates maxgates],1,1); %number of gates
nogarr = zeros(maxgates,1);
for i=1:maxgates
    nogarr(i) = i;
end
xpositive = 800; %positive x axis boundary
xnegative = 0; %negative x axis boundary
ypositive = 600; %positive y axis boundary
ynegative = 0; %negative y axis boundary
centroidofgates = zeros(nog,2); 

angles = [-15; 0; 15]; %5 random angles between left green gate and right red gate
distbtwgates = 100; %fixed distance of 100 decimeter between green and red gate
%% Entrance and Exit Gates
entrancebluegate = [20,110]; %fixed entrance blue gate coordinates 
randangle1 = randi([3 3],1,1);

xtemp1 = entrancebluegate(1) + distbtwgates * sin(angles(randangle1)*pi/180);
ytemp1 = entrancebluegate(2) - distbtwgates * cos(angles(randangle1)*pi/180);
entranceredgate = [xtemp1 ytemp1]; %generation of corresponding entrance red gate varying the angle

xtemp1 = entranceredgate(1)+ entrancebluegate(1);
ytemp1 = entranceredgate(2)+ entrancebluegate(2);
sp = [xtemp1/2,ytemp1/2]; %centroid of entrance gates and starting point of USV
 
exitbluegate = [760,530]; %fixed exit blue gate coordinates 
randangle2 = randi([3 3],1,1);

xtemp2 = exitbluegate(1) + distbtwgates * sin(angles(randangle2)*pi/180);
ytemp2 = exitbluegate(2) - distbtwgates * cos(angles(randangle2)*pi/180);
exitredgate = [xtemp2 ytemp2]; %generation of corresponding exit red gate varying the angle

xtemp2 = exitredgate(1)+ exitbluegate(1);
ytemp2 = exitredgate(2)+ exitbluegate(2);
ep = [xtemp2/2,ytemp2/2]; %centroid of exit gates and ending point of USV
%% Boundary Specifications
x1 = linspace(0,800,20); % boundary totems
upperboundary = zeros(20,2); %upper boundary of course
lowerboundary = zeros(20,2); %lower boundary of course
for i=1:20
    upperboundary(i,1) = x1(i);
    lowerboundary(i,1) = x1(i);
    upperboundary(i,2) = 600;
    lowerboundary(i,2) = 0;
end
%% Gates 
greengate = [100,330;210,560;370,470;520,220;640,450]; %fixed coordinates of green gates
min1 = 1;
max1 = 5;
redgate = zeros(5,2);
for i=1:5
    randangle3 = randi([3 3],1,1);
    xtemp3 = greengate(i,1) + distbtwgates * sin(angles(randangle3)*pi/180);
    ytemp3 = greengate(i,2) - distbtwgates * cos(angles(randangle3)*pi/180); 
    redgate(i,1) = xtemp3;
    redgate(i,2) = ytemp3;
    %Generation of corresponding Red Gates
end
%% Gates Selection 
finalgatesselection = zeros(nog,4);
randomgates = zeros(nog,1);
s=1;
while (s~=nog+1)
    randomgates(s) = randi([min1 max1],1,1);
    if nogarr(randomgates(s))==0
        randomgates(s) = 0;
    else
        nogarr(randomgates(s)) = 0;
        finalgatesselection(s,1)= greengate(randomgates(s),1);
        finalgatesselection(s,2)= greengate(randomgates(s),2);
        finalgatesselection(s,3)= redgate(randomgates(s),1);
        finalgatesselection(s,4)= redgate(randomgates(s),2);
        s=s+1;
        %Final Gates Selection for USV course
    end
end
finalgatesselection = sortrows(finalgatesselection);
%% Centroid Generation
for i=1:nog
    xtemp3 = finalgatesselection(i,1)+ finalgatesselection(i,3);
    ytemp3 = finalgatesselection(i,2)+ finalgatesselection(i,4);
    for j=1:2
        if j==1
            centroidofgates(i,j) = xtemp3/2;
        else
            centroidofgates(i,j) = ytemp3/2;
        end
    end
end
centroidofgates = sortrows(centroidofgates); %Generation of centroid of each pair of gates
%% Obstacles
obs = [86,135;178,384;295,450;427,286;598,297;717,434];%Fixed coordinates of Obstacles
obs = sortrows(obs);
noofobs = [2;4;6]; 
randnoofobs = randi([3 3],1,1); %Random generation of obstacles every course
finalnoofobs = noofobs(randnoofobs);
finalobsselection = zeros(noofobs(randnoofobs),2);
s=1;
min2 = 1;
max2 = 6;
noofobsarr = zeros(6,1);
for i=1:6
    noofobsarr(i)=i;
end
randomobs = zeros(noofobs(randnoofobs),1);
while (s~=noofobs(randnoofobs)+1)
    randomobs(s) = randi([min2 max2],1,1);
    if noofobsarr(randomobs(s))==0
        randomobs(s) = 0;
    else
        noofobsarr(randomobs(s))=0;
        finalobsselection(s,1) = obs(randomobs(s),1);
        finalobsselection(s,2) = obs(randomobs(s),2);
        s=s+1;
        %Selection of final obstacles for the USV course
    end
end
finalobsselection = sortrows(finalobsselection);
% %% Iterative Path Generation
finalcentroidofgates = zeros(nog+2,2);
finalcentroidofgates(1,1) = sp(1,1);
finalcentroidofgates(1,2) = sp(1,2);
for i = 2:nog+1
    finalcentroidofgates(i,1) = centroidofgates(i-1,1);
    finalcentroidofgates(i,2) = centroidofgates(i-1,2);
end
finalcentroidofgates(nog+2,1) = ep(1,1);
finalcentroidofgates(nog+2,2) = ep(1,2);
leftgate = zeros(nog+2,2);
rightgate = zeros(nog+2,2);
leftgate(1,1) = entrancebluegate(1,1);
leftgate(1,2) = entrancebluegate(1,2);
rightgate(1,1) = entranceredgate(1,1);
rightgate(1,2) = entranceredgate(1,2);
for i = 2:nog+1
    leftgate(i,1) = finalgatesselection(i-1,1);
    leftgate(i,2) = finalgatesselection(i-1,2);
    rightgate(i,1) = finalgatesselection(i-1,3);
    rightgate(i,2) = finalgatesselection(i-1,4);
end
leftgate(nog+2,1) = exitbluegate(1,1);
leftgate(nog+2,2) = exitbluegate(1,2);
rightgate(nog+2,1) = exitredgate(1,1);
rightgate(nog+2,2) = exitredgate(1,2);

option = 1; %Water Current Scenario
headingangle = 22;
USVvelocity = 5;
points = [];
nVar =7;
p = ((nog+1)*nVar)+nog+2;
finalpoints = zeros(p,2);
finalcost = zeros(nog+1,3);
no = 1;
c=1;

for j=1:nog+1
    counter =1;
    for k=1:finalnoofobs
        if finalobsselection(k,1)<finalcentroidofgates(j+1,1) 
            if finalobsselection(k,1)>finalcentroidofgates(j,1)
                finalobs(counter,1) = finalobsselection(k,1);
                finalobs(counter,2) = finalobsselection(k,2);
                counter = counter +1;
            end
        end
    end
    [points, cost] = Traj_MOIGWO(finalcentroidofgates(j,1:2),finalcentroidofgates(j+1,1:2),leftgate(j,1:2),rightgate(j,1:2),leftgate(j+1,1:2),rightgate(j,1:2),finalobs,USVvelocity,option,headingangle,nVar);
    for i=1:nVar
        finalpoints(no,1) = points(i,1);
        finalpoints(no,2) = points(i,2);
        no=no+1;
    end
    for n=1:3
        finalcost(c,n) = cost(1,n);
    end
    c=c+1;
end

for j=1:nog+2
    finalpoints(no,1) = finalcentroidofgates(j,1);
    finalpoints(no,2) = finalcentroidofgates(j,2);
    no = no+1;
end
finalpoints = sortrows(finalpoints);
%% Plotting of Environment
if drawing_map ==1
    figure(1)
    hold on
    angles = linspace(0, 2*pi, 500);
    radius = 24.25;
    labels = {'SP';'WP1';'WP2';'WP3';'WP4';'WP5';'EP'};

    x1 = radius * cos(angles) + entrancebluegate(1);
    y1 = radius * sin(angles) + entrancebluegate(2);
    plot(x1, y1, 'b-', 'LineWidth', 2);
    c1 = plot(entrancebluegate(1),entrancebluegate(2),'o','MarkerSize',6.5,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1]);

    x2 = radius * cos(angles) + entranceredgate(1);
    y2 = radius * sin(angles) + entranceredgate(2);
    plot(x2, y2, 'r-', 'LineWidth', 2);
    c2 = plot(entranceredgate(1),entranceredgate(2),'o','MarkerSize',6.5,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0]);

    x3 = radius * cos(angles) + exitbluegate(1);
    y3 = radius * sin(angles) + exitbluegate(2);
    plot(x3, y3, 'b-', 'LineWidth', 2);
    plot(exitbluegate(1),exitbluegate(2),'o','MarkerSize',6.5,'MarkerEdgeColor',[0,0,1],'MarkerFaceColor',[0,0,1])

    x4 = radius * cos(angles) + exitredgate(1);
    y4 = radius * sin(angles) + exitredgate(2);
    plot(x4, y4, 'r-', 'LineWidth', 2);
    plot(exitredgate(1),exitredgate(2),'o','MarkerSize',6.5,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0])

    for k=1:nog
        c5 = plot(centroidofgates(k,1),centroidofgates(k,2),'*','MarkerSize',5,'MarkerEdgeColor',[0.25,0.25,0.25],'MarkerFaceColor',[0.25,0.25,0.25]);
    end

    for k=1:nog
        text(centroidofgates(k,1),centroidofgates(k,2),labels(k+1),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center','FontSize',10)
    end 

    for k=1:nog
        c8 = plot(finalgatesselection(k,1),finalgatesselection(k,2),'o','MarkerSize',6.5,'MarkerEdgeColor',[0,1,0],'MarkerFaceColor',[0,1,0]);
        x5 = radius * cos(angles) + finalgatesselection(k,1);
        y5 = radius * sin(angles) + finalgatesselection(k,2);
        plot(x5, y5, 'g-', 'LineWidth', 2);

        c9 = plot(finalgatesselection(k,3),finalgatesselection(k,4),'o','MarkerSize',6.5,'MarkerEdgeColor',[1,0,0],'MarkerFaceColor',[1,0,0]);
        x6 = radius * cos(angles) + finalgatesselection(k,3);
        y6 = radius * sin(angles) + finalgatesselection(k,4);
        plot(x6, y6, 'r-', 'LineWidth', 2);
    end

    c6 = plot(sp(1,1),sp(1,2),'*','MarkerSize',5,'MarkerEdgeColor',[0.25,0.25,0.25],'MarkerFaceColor',[0.25,0.25,0.25]);
    c7 = plot(ep(1,1),ep(1,2),'*','MarkerSize',5,'MarkerEdgeColor',[0.25,0.25,0.25],'MarkerFaceColor',[0.25,0.25,0.25]);

    text(sp(1,1),sp(1,2),labels(1),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center','FontSize',10)
    text(ep(1,1),ep(1,2),labels(7),'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center','FontSize',10)

    for k=1:noofobs(randnoofobs)
        c4 = plot(finalobsselection(k,1),finalobsselection(k,2),'o','MarkerSize',4.7,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[0,0,0]);
        x7 = radius * cos(angles) + finalobsselection(k,1);
        y7 = radius * sin(angles) + finalobsselection(k,2);
        plot(x7, y7, 'k-', 'LineWidth', 2);
    end

    c3 = plot(upperboundary(:,1),upperboundary(:,2),'-o','MarkerSize',4.7,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980]);
    plot(lowerboundary(:,1),lowerboundary(:,2),'-o','MarkerSize',4.7,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
    
    c10 = plot(finalpoints(:,1),finalpoints(:,2),'-','LineWidth', 1.5);
    
    axis([0 800 -10 610]);
    x = [0 0 800 800];
    y = [-10 610 610 -10];
    h = fill(x,y,'b');
    set(h,'FaceAlpha',0.1);

    xlabel('Length of course(x)','FontSize',14,'FontName','Calibri','FontWeight','bold')
    ylabel('Length of Vertical Boundary (y)','FontSize',14,'FontName','Calibri','FontWeight','bold')
    title('Task 2 - Maritime RobotX Challenge 2022','FontSize',14,'FontName','Calibri','FontWeight','bold')
    hold off

    subset = [c1;c2;c3;c4;c5;c6;c7;c8;c9;c10];
    legend(subset,'Entrance/Exit - Left Blue Gate','Entrance/Exit - Right Red Gate','Boundary','Obstacle','WP - Waypoint','SP - Starting Point','EP - End Point','Green Gate - Left','Red Gate - Right','Traj-MOIGWO')
    legend('Location','northwest','Orientation','vertical','FontSize',10,'FontName','Calibri','FontWeight','bold')
    legend boxoff
end

 
