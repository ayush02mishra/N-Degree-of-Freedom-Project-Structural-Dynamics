function [K, M] = Assemble(N,dof,L,m,mType,EI,id)
%
%       Input parameters
%       N = number of elements
%       dof = degree of freedon, L = length of all elements, 
%       m = mass, mType = mass distribution, EI = Matrix of all elemental
%       stiffness, t = time Array, id = nodes for each element,
%       F1, F2 and f0 = external force parameters, omega = external force
%       frequency.
%


%
K = zeros(dof*N,dof*N);
switch mType
    case 1
        M = m*eye(N);
    case 2
        M = diag(m:m/N:m/N);
end

k = Elemental(EI(1),L);
K (1,1) = k(1,1);
k = Elemental(EI(N),L);
K (N,N) = k(2,2);
for i = 2:N;
    K(id(i-1,:),id(i-1,:)) = K(id(i-1,:),id(i-1,:)) + Elemental(EI(i),L);
end