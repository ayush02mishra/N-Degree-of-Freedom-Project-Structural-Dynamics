function [S,yield] = ElastoplasticStress(u,k,S0)
S = zeros(numel(u),1);
yield = zeros(numel(u),1);
S(1) = k(1)*u(1);
UPlastic = 0;
for i = 2:numel(u)
    S(i) = k(i)*((u(i)-UPlastic));
    if abs(S(i))>S0(i)
        yield(i) = 1;
        S(i) = (sign(S(i)))*S0(i);
        UPlastic = u(i)-(S0(i)/k(i))*sign(S(i));
    end
end