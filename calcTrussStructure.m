function [u,S]=calcTrussStructure(EA,nNode,nTruss,coord,conn,boundaryCond,force)

% Definition of the degrees of freedom (dofs)
alldofs=1:2*nNode;
fixeddofs=[];
for i=1: size(boundaryCond, 1)
    if boundaryCond(i,2)==1
        dofs=boundaryCond(i,1)*2-1;
    else
        dofs=boundaryCond(i,1)*2;
    end
    fixeddofs=union(fixeddofs,dofs);
end
freedofs=setdiff(alldofs,fixeddofs);

% Setup of the force vector
f=zeros(2*nNode,1);
for i=1: size(force, 1)
    f(force(i,1)*2-1)=force(i,2);
    f(force(i,1)*2)=force(i,3);
end

% Calculating truss lengths LL
LL=((coord(conn(1:nTruss,2),1)-coord(conn(1:nTruss,1),1)).^2+(coord(conn(1:nTruss,2),2)...
    -coord(conn(1:nTruss,1),2)).^2).^0.5;

% Calculating orientation/rotation of each truss
for i=1:nTruss
    y(i+1)=coord(conn(i,2),2);
    y(i)=coord(conn(i,1),2);
    x(i+1)=coord(conn(i,2),1);
    x(i)=coord(conn(i,1),1);
    
    if (y(i+1)>y(i)) && (x(i+1)>x(i))
        v(i)=atan((y(i)-y(i+1))/(x(i)-x(i+1)));
    end
    if (y(i+1)>y(i)) && (x(i+1)<x(i))
        v(i)=pi+atan((y(i)-y(i+1))/(x(i)-x(i+1)));
    end
    if (y(i+1)<y(i)) && (x(i+1)<x(i))
        v(i)=pi+atan((y(i)-y(i+1))/(x(i)-x(i+1)));
    end
    if (y(i+1)<y(i)) && (x(i+1)>x(i))
        v(i)=atan((y(i)-y(i+1))/(x(i)-x(i+1)));
    end
    
    if (y(i+1)==y(i)) && (x(i+1)>x(i))
        v(i)=0;
    end
    if (y(i+1)==y(i)) && (x(i+1)<x(i))
        v(i)=pi;
    end
    if (y(i+1)<y(i)) && (x(i+1)==x(i))
        v(i)=pi+pi/2;
    end
    if (y(i+1)>y(i)) && (x(i+1)==x(i))
        v(i)=pi/2;
    end
end

% Assembling of the stiffness matrix
K=zeros(length(alldofs),length(alldofs));
for i=1:nTruss
    
    % Allocation matrix
    A=zeros(2,2*nNode);
    for j=1:2
        A((2*(j-1)+1),(2*(conn(i,j))-1))=1;
        A((2*(j-1)+2),(2*(conn(i,j))))=1;
    end
    cellA{i,1}=A;
    
    % Transformation matrix
    D=zeros(2,4);
    D=[cos(v(i)) sin(v(i)) 0 0;...
        0 0 cos(v(i)) sin(v(i))];
    cellD{i,1}=D;
    
    % Element stiffness matrix
    ks=EA(i)/LL(i)*[1 -1;...
        -1 1];
    
    % Summing the stiffness matrices
    K=K+A'*D'*ks*D*A;
end

% Calculating the deformations
u(fixeddofs)=0;
u(freedofs)=K(freedofs,freedofs)\f(freedofs);

% Calculating the truss forces
for i=1:nTruss
    S(i)=EA(i)/LL(i)*[-1 1]*cellD{i,1}*cellA{i,1}*u';
end
