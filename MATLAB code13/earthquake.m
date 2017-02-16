% Shp = [b,n,m,c]'
%   function quake = earthquake(EQtype,Shp,N,Peak,Amp,Frq,Phs,tD,t)

function quake = earthquake(EQtype,Shp,N,Peak,Amp,Frq,Phs,tD,t)
Amp = Amp./norm(Amp);
if EQtype == 1
    tao = Shp(1)*t/(Shp(2)*tD);
    shape = Peak*(tao.^Shp(2)).*exp(Shp(2).*(1-tao));
else
    if EQtype == 2
    tao = t/tD;
    a = 2*Shp(3)*pi/Shp(1);
    tao0 = pi/(2*a);
    shape = Peak*((1+1/(4*a*a)).*((sin(a.*tao)).^2).*exp(tao0-tao) + Shp(4).*tao.*(1-tao));
    end
end
quake1 = zeros(numel(t),1);
% plot(t,shape);
figure;
for time = 1:numel(t);
    for i = 1:N;
        quake1(time) = quake1(time) + shape(time)*Amp(i)*sin(Frq(i)*time + Phs(i));
    end
end
quake = quake1;
% plot(t,quake);
end