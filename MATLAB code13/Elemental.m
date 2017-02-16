function k = Elemental(EI,L)
%
%EI for each element is taken as input and Elemental stiffness matrix is
%returned


%   k matrix assembly
k = [24*EI/(L^3) -24*EI/(L^3); -24*EI/(L^3) 24*EI/(L^3)];
