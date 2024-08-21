        clc;
clear;
close all;

%Mohammadaadil Munvvarbhai Shaikh - 23282106 
%Mohammad Ameer Sohail - 23287773 
%Prajul Mullookkaran Pazhayapurayil - 23284633
%Athul Krishna Nalumakkal Sahul - 23233858 


[filename, pathname] = uigetfile;
filepath = [pathname, filename]; 
input = load(filepath);
conn = input.conn;
coord = input.coord;
force = input.force;
boundaryCond = input.boundaryCond;

x = coord(:,1);
y = coord(:,2);
nTruss = size(conn,1); 
nNode = size(x,1);
num_forces = size(force,1);
nBoundaryCond = size(boundaryCond,1);


f1 = figure;
hold on;

title(filename);
for i = 1:nTruss
    a = conn(i,1);
    b = conn(i, 2);
    x1 = coord(a,1);
    x2 = coord(b,1);
    y1 = coord(a,2);
    y2 = coord(b,2);
    plot([x1 x2], [y1 y2], '-ok');
    text((x1 + x2) / 2, (y1 + y2) / 2, sprintf('%d', i));
end

for i = 1: nBoundaryCond
    a = boundaryCond(i,1);
    b = boundaryCond(i,2);
    x = coord(a,1);
    y = coord(a,2);
    if b == 1
        p2 = plot(x - 0.1, y, 'r>'); %plotting lower 0.1 to the x
        %r> means red > triangle of this style
    else
        p2 = plot(x, y - 0.1, 'r^'); % plotting left 0.1 to the y
    end
end

for i = 1: num_forces
    a = force(i,1);
    x = coord(a,1);
    y = coord(a,2);
    %'MaxHeadSize', 1 / norm(F(i,2:3): sets the maximum size of the arrow head relative to the length of the force vector.
    p3 = quiver(x, y, force(i,2), force(i,3), 'MaxHeadSize', 1 / norm(force(i,2:3)), 'color', 'g');
end

axis equal;
hold off;

%Calculating vector EA (axial stiffness for each bar)
E = 3; % KN / mm2 (Young Modulus)
w = 5; % width of the bar lets assume
h = 5;
A = zeros(nTruss,1);
EA = zeros(nTruss,1);

for i = 1:nTruss
    a = conn(i,1);
    b = conn(i,2);
    x1 = coord(a,1);
    x2 = coord(b,1);
    y1 = coord(a,2);
    y2 = coord(b,2);
    len_bar = sqrt((x2 - x1)^2 + (y2 - y1)^2);
    A(i) = w * h;
    EA (i) = E * A(i);  
end

%Calling calcTrussStructure function 
% output = u--> displacement in x and y directions, S --> Internal bar force
[u,S]=calcTrussStructure(EA, nNode, nTruss, coord, conn, boundaryCond, force);

%plottig the deflected beam

% Separate the deflections in x and y directions
u_x = u(1:2:end);
u_y = u(2:2:end);

deformed_coord = coord + [u_x' u_y'];

f2 = figure;
hold on;

title('Deformed and undeformed structure')
for i = 1:nTruss
    a = conn(i,1);
    b = conn(i, 2);
    x1 = coord(a,1);
    x2 = coord(b,1);
    y1 = coord(a,2);
    y2 = coord(b,2);
    p1 = plot([x1 x2], [y1 y2], '-ok');
    text((x1 + x2) / 2, (y1 + y2) / 2, sprintf('%d', i));
end

for i = 1: nBoundaryCond
    a = boundaryCond(i,1);
    b = boundaryCond(i,2);
    x = coord(a,1);
    y = coord(a,2);
    if b == 1
        p2 = plot(x - 0.1, y, 'r>'); %plotting lower 0.1 to the x
        %r> means red > triangle of this style
    else
        p2 = plot(x, y - 0.1, 'r^'); % plotting left 0.1 to the y
    end
end


for i = 1: num_forces
    a = force(i,1);
    x = coord(a,1);
    y = coord(a,2);
    %'MaxHeadSize', 1 / norm(F(i,2:3): sets the maximum size of the arrow head relative to the length of the force vector.
    p3 = quiver(x, y, force(i,2), force(i,3), 'MaxHeadSize', 1 / norm(force(i,2:3)), 'color', 'g');
end

%plotting deformed beam

for i = 1:nTruss
    a = conn(i,1);
    b = conn(i,2);
    x1 = deformed_coord(a,1);
    x2 = deformed_coord(b,1);
    y1 = deformed_coord(a,2);
    y2 = deformed_coord(b,2);
    p4 = plot([x1 x2], [y1 y2], 'LineStyle','-', 'Color', 'blue');
end

% for i = 1: num_forces
%     a = force(i,1);
%     x = deformed_coord(a,1);
%     y = deformed_coord(a,2);
%     %'MaxHeadSize', 1 / norm(F(i,2:3): sets the maximum size of the arrow head relative to the length of the force vector.
%     p3 = quiver(x, y, force(i,2), force(i,3), 'MaxHeadSize', 1 / norm(force(i,2:3)), 'color', 'r');
% end

legend([p1,p2,p3,p4], 'Undeformed beam','Boundary conditions','forces','Deformed beam')
axis equal;
hold off;

%% Task 2

lowerBound = 0.01 * ones(nTruss, 1); % Lower bound vector
upperBound = 5 * ones(nTruss, 1);    % Upper bound vector


% initiate the search for the optimal scaling factors effectively, ensuring the process
% is grounded within the feasible region defined by the problem's constraints
startVector = rand(nTruss, 1) .* (upperBound - lowerBound) + lowerBound;

A =[];          %Nonlinear inequality constraints
b = [];         %linear inequality constraints
A_eqa = [];
b_eqa = [];


% function handle for an anonymous function
x_opt = fmincon(@(x_opt)objectiveFunction(x_opt,EA,nNode,nTruss,coord,conn,boundaryCond,force),startVector, A, b,A_eqa,b_eqa, lowerBound, upperBound);

EA_opt = x_opt .* EA;

[u_opt, s] = calcTrussStructure(EA_opt,nNode,nTruss,coord,conn,boundaryCond,force);

%seperating displacement in x and y direction

% Initialize empty vectors for x and y displacements
ux_opt = zeros(1, length(u_opt) / 2);
uy_opt = zeros(1, length(u_opt) / 2);

% using for loop
for i = 1:length(u_opt)
    if mod(i, 2) == 1
        % Odd indices correspond to x-direction displacements
        ux_opt((i + 1) / 2) = u_opt(i);
    else
        % Even indices correspond to y-direction displacements
        uy_opt(i / 2) = u_opt(i);
    end
end

%plotting the graph
f3 = figure;
hold on
title('Optimized Structure')

%Undeformed beam
for i = 1:nTruss
    a = conn(i,1);
    b = conn(i, 2);
    x1 = coord(a,1);
    x2 = coord(b,1);
    y1 = coord(a,2);
    y2 = coord(b,2);
    p1 = plot([x1 x2], [y1 y2], '-ok');
    text((x1 + x2) / 2, (y1 + y2) / 2, sprintf('%d', i));
end

for i = 1: nBoundaryCond
    a = boundaryCond(i,1);
    b = boundaryCond(i,2);
    x = coord(a,1);
    y = coord(a,2);
    if b == 1
        p2 = plot(x - 0.1, y, 'r>'); %plotting lower 0.1 to the x
        %r> means red > triangle of this style
    else
        p2 = plot(x, y - 0.1, 'r^'); % plotting left 0.1 to the y
    end
end


for i = 1: num_forces
    a = force(i,1);
    x = coord(a,1);
    y = coord(a,2);
    %'MaxHeadSize', 1 / norm(F(i,2:3): sets the maximum size of the arrow head relative to the length of the force vector.
    p3 = quiver(x, y, force(i,2), force(i,3), 'MaxHeadSize', 1 / norm(force(i,2:3)), 'color', 'g');
end

%Deformed beam with optimized displacements
for i = 1:nTruss
    a = conn(i,1);
    b = conn(i,2);
    x1 = deformed_coord(a,1);
    x2 = deformed_coord(b,1);
    y1 = deformed_coord(a,2);
    y2 = deformed_coord(b,2);
    p4 = plot([x1 x2], [y1 y2], 'LineStyle','-', 'Color', 'blue');
end

%plotting  optimizied Truss
deformed_coord_opt = coord + [ux_opt' uy_opt'];

for i = 1:nTruss
    a = conn(i,1);
    b = conn(i,2);
    x1 = deformed_coord_opt(a,1);
    x2 = deformed_coord_opt(b,1);
    y1 = deformed_coord_opt(a,2);
    y2 = deformed_coord_opt(b,2);
    p5 = plot([x1 x2], [y1 y2], 'LineStyle','-', 'Color', 'r', 'LineWidth',x_opt(i));
end


legend([p1, p2, p3, p4,p5] , 'Undeformed Beam', 'Boundary Conditions', 'forces', 'Deformed Truss', 'optimized structure')
axis equal;
hold off;

%% Task 3 
%Minimization of the y-displacement at the force application point with equality (mass) constraints


A =[];          %linear inequality constraints  Ax =b
b = [];         %linear inequality constraints
A_eqa = EA';    %linear equality constraints   A_eqa * x = b_eqa
b_eqa = 825;


x_opt_3 = fmincon(@(x_opt)objectiveFunction(x_opt,EA,nNode,nTruss,coord,conn,boundaryCond,force),startVector, A, b,A_eqa,b_eqa, lowerBound, upperBound);

EA_opt_3 = x_opt_3 .* EA;

[u_opt_3, s] = calcTrussStructure(EA_opt_3,nNode,nTruss,coord,conn,boundaryCond,force);

%seperating displacement in x and y direction

% Initialize empty vectors for x and y displacements
ux_opt_3 = zeros(1, length(u_opt) / 2);
uy_opt_3 = zeros(1, length(u_opt) / 2);

% using for loop
for i = 1:length(u_opt_3)
    if mod(i, 2) == 1
        % Odd indices correspond to x-direction displacements
        ux_opt_3((i + 1) / 2) = u_opt_3(i);
    else
        % Even indices correspond to y-direction displacements
        uy_opt_3(i / 2) = u_opt_3(i);
    end
end

%plotting the graph
f3 = figure;
hold on
title('Optimized Structure with equality (mass) constraints')

%Undeformed beam
for i = 1:nTruss
    a = conn(i,1);
    b = conn(i, 2);
    x1 = coord(a,1);
    x2 = coord(b,1);
    y1 = coord(a,2);
    y2 = coord(b,2);
    p1 = plot([x1 x2], [y1 y2], '-ok');
    text((x1 + x2) / 2, (y1 + y2) / 2, sprintf('%d', i));
end

for i = 1: nBoundaryCond
    a = boundaryCond(i,1);
    b = boundaryCond(i,2);
    x = coord(a,1);
    y = coord(a,2);
    if b == 1
        p2 = plot(x - 0.1, y, 'r>'); %plotting lower 0.1 to the x
        %r> means red > triangle of this style
    else
        p2 = plot(x, y - 0.1, 'r^'); % plotting left 0.1 to the y
    end
end


for i = 1: num_forces
    a = force(i,1);
    x = coord(a,1);
    y = coord(a,2);
    %'MaxHeadSize', 1 / norm(F(i,2:3): sets the maximum size of the arrow head relative to the length of the force vector.
    p3 = quiver(x, y, force(i,2), force(i,3), 'MaxHeadSize', 1 / norm(force(i,2:3)), 'color', 'g');
end

%Deformed beam with optimized displacements
for i = 1:nTruss
    a = conn(i,1);
    b = conn(i,2);
    x1 = deformed_coord(a,1);
    x2 = deformed_coord(b,1);
    y1 = deformed_coord(a,2);
    y2 = deformed_coord(b,2);
    p4 = plot([x1 x2], [y1 y2], 'LineStyle','-', 'Color', 'blue');
end


%plotting optimized Truss
deformed_coord_opt_3 = coord + [ux_opt_3' uy_opt_3'];

for i = 1:nTruss
    a = conn(i,1);
    b = conn(i,2);
    x1 = deformed_coord_opt_3(a,1);
    x2 = deformed_coord_opt_3(b,1);
    y1 = deformed_coord_opt_3(a,2);
    y2 = deformed_coord_opt_3(b,2);
    p5 = plot([x1 x2], [y1 y2], 'LineStyle','-', 'Color', 'r', 'LineWidth',x_opt_3(i));
end


legend([p1, p2, p3, p4,p5] , 'Undeformed Beam', 'Boundary Conditions', 'forces','Deformed Truss' ,'optimized structure')
axis equal;
hold off;
%% Task 4
% Minimization of the y-displacement at the force application point with equality (mass) and inequality(stress) constraints
A_eqa = [];
b_eqa = [];
sum_EA = 500;       %reduced total stiffness
A = w * h;
 

[x_opt4] = fmincon(@(x_opt_3)objectiveFunction(x_opt_3,EA,nNode,nTruss,coord,conn,boundaryCond,force), startVector, [], [], A_eqa, b_eqa,lowerBound,upperBound,@(x_opt_3)constraintFunction(x_opt_3,EA,A,s,nTruss));
 
%scaling axial stiffness 
EA_opt_4 = x_opt4.*EA;
[u_opt4,s]=calcTrussStructure(EA_opt_4,nNode,nTruss,coord,conn,boundaryCond,force);

%seperating displacement in x and y direction

% Initialize empty vectors for x and y displacements
ux_opt_4 = zeros(1, length(u_opt4) / 2);
uy_opt_4 = zeros(1, length(u_opt4) / 2);

% using for loop
for i = 1:length(u_opt4)
    if mod(i, 2) == 1
        % Odd indices correspond to x-direction displacements
        ux_opt_4((i + 1) / 2) = u_opt4(i);
    else
        % Even indices correspond to y-direction displacements
        uy_opt_4(i / 2) = u_opt4(i);
    end
end

%plotting the graph
f3 = figure;
hold on
title('Optimized Structure with equality (mass) and inequality constraints')

%Undeformed beam
for i = 1:nTruss
    a = conn(i,1);
    b = conn(i, 2);
    x1 = coord(a,1);
    x2 = coord(b,1);
    y1 = coord(a,2);
    y2 = coord(b,2);
    p1 = plot([x1 x2], [y1 y2], '-ok');
    text((x1 + x2) / 2, (y1 + y2) / 2, sprintf('%d', i));
end

for i = 1: nBoundaryCond
    a = boundaryCond(i,1);
    b = boundaryCond(i,2);
    x = coord(a,1);
    y = coord(a,2);
    if b == 1
        p2 = plot(x - 0.1, y, 'r>'); %plotting lower 0.1 to the x
        %r> means red > triangle of this style
    else
        p2 = plot(x, y - 0.1, 'r^'); % plotting left 0.1 to the y
    end
end


for i = 1: num_forces
    a = force(i,1);
    x = coord(a,1);
    y = coord(a,2);
    %'MaxHeadSize', 1 / norm(F(i,2:3): sets the maximum size of the arrow head relative to the length of the force vector.
    p3 = quiver(x, y, force(i,2), force(i,3), 'MaxHeadSize', 1 / norm(force(i,2:3)), 'color', 'g');
end

%Plotting Deformed Truss
for i = 1:nTruss
    a = conn(i,1);
    b = conn(i,2);
    x1 = deformed_coord(a,1);
    x2 = deformed_coord(b,1);
    y1 = deformed_coord(a,2);
    y2 = deformed_coord(b,2);
    p4 = plot([x1 x2], [y1 y2], 'LineStyle','-', 'Color', 'blue');
end

%plotting optimized Truss
deformed_coord_opt_4 = coord + [ux_opt_4' uy_opt_4'];

for i = 1:nTruss
    a = conn(i,1);
    b = conn(i,2);
    x1 = deformed_coord_opt_4(a,1);
    x2 = deformed_coord_opt_4(b,1);
    y1 = deformed_coord_opt_4(a,2);
    y2 = deformed_coord_opt_4(b,2);
    p5 = plot([x1 x2], [y1 y2], 'LineStyle','-', 'Color', 'r', 'LineWidth',x_opt4(i));
end


legend([p1, p2, p3, p4, p5] , 'Undeformed Beam', 'Boundary Conditions', 'forces', 'Deformed Truss','optimized structure')
axis equal;
hold off;

%% Task 5
%  Evaluation of the efficiency of different optimization algorithms
%MaxIter = Maximum number of iteration for slected sqp Algorithm

A_eqa = [];
b_eqa = [];
sum_EA = 500;       %reduced total stiffness
A = w * h;
 
algorithm = optimoptions(@fmincon,'Display','iter','Algorithm','sqp','MaxIter',20);
[x_opt5] = fmincon(@(x_opt_3)objectiveFunction(x_opt_3,EA,nNode,nTruss,coord,conn,boundaryCond,force), startVector, [], [], A_eqa, b_eqa,lowerBound,upperBound,@(x_opt_3)constraintFunction(x_opt_3,EA,A,s,nTruss),algorithm);
 
%scaling axial stiffness 
EA_opt_5 = x_opt4.*EA;
[u_opt5,s]=calcTrussStructure(EA_opt_4,nNode,nTruss,coord,conn,boundaryCond,force);

%seperating displacement in x and y direction

% Initialize empty vectors for x and y displacements
ux_opt_5 = zeros(1, length(u_opt5) / 2);
uy_opt_5 = zeros(1, length(u_opt5) / 2);

% using for loop
for i = 1:length(u_opt5)
    if mod(i, 2) == 1
        % Odd indices correspond to x-direction displacements
        ux_opt_5((i + 1) / 2) = u_opt5(i);
    else
        % Even indices correspond to y-direction displacements
        uy_opt_5(i / 2) = u_opt5(i);
    end
end

%plotting the graph
f3 = figure;
hold on
title('Optimized Structure by SQP(Sequential Quadratic Programming)')

%Undeformed beam
for i = 1:nTruss
    a = conn(i,1);
    b = conn(i, 2);
    x1 = coord(a,1);
    x2 = coord(b,1);
    y1 = coord(a,2);
    y2 = coord(b,2);
    p1 = plot([x1 x2], [y1 y2], '-ok');
    text((x1 + x2) / 2, (y1 + y2) / 2, sprintf('%d', i));
end

for i = 1: nBoundaryCond
    a = boundaryCond(i,1);
    b = boundaryCond(i,2);
    x = coord(a,1);
    y = coord(a,2);
    if b == 1
        p2 = plot(x - 0.1, y, 'r>'); %plotting lower 0.1 to the x
        %r> means red > triangle of this style
    else
        p2 = plot(x, y - 0.1, 'r^'); % plotting left 0.1 to the y
    end
end


for i = 1: num_forces
    a = force(i,1);
    x = coord(a,1);
    y = coord(a,2);
    %'MaxHeadSize', 1 / norm(F(i,2:3): sets the maximum size of the arrow head relative to the length of the force vector.
    p3 = quiver(x, y, force(i,2), force(i,3), 'MaxHeadSize', 1 / norm(force(i,2:3)), 'color', 'g');
end

%Plotting Deformed Truss
for i = 1:nTruss
    a = conn(i,1);
    b = conn(i,2);
    x1 = deformed_coord(a,1);
    x2 = deformed_coord(b,1);
    y1 = deformed_coord(a,2);
    y2 = deformed_coord(b,2);
    p4 = plot([x1 x2], [y1 y2], 'LineStyle','-', 'Color', 'blue');
end

%plotting optimized Truss
deformed_coord_opt_5 = coord + [ux_opt_5' uy_opt_5'];

for i = 1:nTruss
    a = conn(i,1);
    b = conn(i,2);
    x1 = deformed_coord_opt_5(a,1);
    x2 = deformed_coord_opt_5(b,1);
    y1 = deformed_coord_opt_5(a,2);
    y2 = deformed_coord_opt_5(b,2);
    p5 = plot([x1 x2], [y1 y2], 'LineStyle','-', 'Color', 'r', 'LineWidth',x_opt5(i));
end


legend([p1, p2, p3, p4, p5] , 'Undeformed Beam', 'Boundary Conditions', 'forces', 'Deformed Truss','optimized structure')
axis equal;
hold off;

%% Task 6 Saving Results

save Results_Group37_V4_Task7_crane.mat

%% Task 7 
%optimize a crane’struss structure

%Select V4_Input_crane
%Change the equality constaints value to 2175 in constraintFunction by
%commenting out another value



