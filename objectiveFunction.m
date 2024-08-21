%Mohammadaadil Munvvarbhai Shaikh - 23282106 
%Mohammad Ameer Sohail - 23287773 
%Prajul Mullookkaran Pazhayapurayil - 23284633
%Athul Krishna Nalumakkal Sahul - 23233858 



function Y = objectiveFunction(x_opt,EA,nNode,nTruss,coord,conn,boundaryCond,force)
    %To find a Objective value -Displacement at the point of Application
    % Scale the initial axial stiffness EA with the scaling values
    scaledEA = EA .* x_opt;
    
    % Calculate the displacements using calcTrussStructure 
    [u, s] = calcTrussStructure(scaledEA,nNode,nTruss,coord,conn,boundaryCond,force);
    
    % How he find ??
    % The objective value is the negative magnitude of the displacement
    Y = abs(u(1,2* force(1,1)));
end