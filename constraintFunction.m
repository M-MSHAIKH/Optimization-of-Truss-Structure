%Mohammadaadil Munvvarbhai Shaikh - 23282106 
%Mohammad Ameer Sohail - 23287773 
%Prajul Mullookkaran Pazhayapurayil - 23284633
%Athul Krishna Nalumakkal Sahul - 23233858 



function [c, c_eqa] = constraintFunction(x_opt,EA,A,s,nTruss)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
EA_scale = x_opt .* EA;
A_scale = x_opt .* A;

% For normal truss 
%c_eqa = sum(EA_scale)-500;
% For larger crane system
c_eqa = sum(EA_scale)-2175;

% internal force/scaled area - (normal stress)
c = abs((s'./(x_opt.*A)))-0.080.*ones(nTruss,1);
end