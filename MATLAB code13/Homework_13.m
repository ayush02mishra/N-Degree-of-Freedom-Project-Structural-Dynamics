%   Ayush Mishra
%   Homework 11
%
%   Time Parameters
tD = 50;
dt = .2;
t = 1:dt:tD;
timeSteps = numel(t);

RayleighDamping = 1;        %   Else General Orthogonal Damping

%       Design Parameters
L = 5;
LStorey = 1;
m = 50;
mType = 1;      %   Type of mass distribution: 1 is constant, 2 linearly decreasing with minimum at top.
EI = 50*[1 1 1 1 1 1 1 1 1 1]';

%       Problem Parameters
N = 10;
dof = 1;


%   Earthquake Parameters
Shp = [2 3 3 1]';
Amp = [2 3 5 1 3 4]';
Frq = [1 5 8 30 4 6]';
Phs = [0 0 0 0 0 0]';
Peak = .2;
NQuake = 3;
EQtype = 2;
eqInterval = numel(t)/2;
t0 = t(1:eqInterval);

%       Initial Conditions
u0 = [0 0 0 0 0 0 0 0 0 0]';
v0 = [0 0 0 0 0 0 0 0 0 0]';
%   Assembly
id = [(1:N-1)', (2:N)'];
[K,M] = Assemble(N,dof,L,m,mType,EI,id);

F = zeros(timeSteps,1);
Quake = m*earthquake(EQtype,Shp,NQuake,Peak,Amp,Frq,Phs,eqInterval*timeStep,t0);
F(1:eqInterval) = Quake;
%   Yield Stress
S0 = 300*[1 1 1 1 1 1 1 1 1 1];

%   Modal freq and mode shapes
[V,D]= eig(K,M);
disp('Modal Frequencies');
disp(diag(D));

%       Damping Parameters
if RayleighDamping ==0
DampingModes = [3 6 8];
ModalDampingValues = [100 100 100];
CModal = zeros(N,1);
CModal(DampingModes) = ModalDampingValues;
CModal = diag(CModal);
C = M*V*CModal*V'*M;
else
    mu1 = .01;
    mu2 = .01;
    C = mu1*M + mu2*K;
end

%       Newmarks's Parameters
h = dt;
B = 1/4;
Y = 1/2;
T = h*h*(1/2-B);
KEff = M + T*K + h*(1-Y)*C;


u = zeros(timeSteps,N);
v = zeros(timeSteps,N);
a = zeros(timeSteps,N);
u(1,:) = u0;
v(1,:) = v0;
a(1,:) = 0;


%_______________________________________________________________________
R = zeros(timeSteps,N);
%       Newton's iterator parameters
tol = 10^-4;
maxIteration = 1000;
an= zeros(maxIteration,N);
un= zeros(maxIteration,N);
vn= zeros(maxIteration,N);

for i = 1:timeSteps
    if i == 1
        u(i,:) = u0;
        v(i,:) = v0;
        R(i,:) = ElastoplasticStress(u(i,:),EI,S0);
        CForce = ((v(i,:) + h*Y*a(i,:))*C);
        a(i,:) = (F(i)*ones(N,1)' - R(i,:)-CForce)/M;
    else
        %       --------Newton's iteration
        %initial assumption
        u(i,:) = u(i-1,:);
        R(i,:) = ElastoplasticStress(u(i,:),EI,S0);
        CForce = (v(i-1,:) + h*Y*a(i-1,:))*C;
        a(i,:) = (F(i)*ones(N,1)' - R(i,:)-CForce)/M;
        
        j = 2;
        an(j,:) = a(i,:);
        un(j,:) = u(i,:);
        vn(j,:) = v(i,:);
        while ((j<maxIteration && abs(an(j)-an(j-1))>tol) || j<3)
            j= j+1;
            CForce = (vn(j,:) + h*Y*an(j,:))*C;
            R(i,:) = ElastoplasticStress(un(j,:),EI,S0);
            Cn = vn(j-1) + h*Y*an(j-1);
            Bn = un(j-1) + h*vn(j-1) + h*h*B*an(j-1);
            an(j,:) = an(j-1) - (F(i-1)*ones(N,1)' - R(i,:) - CForce)/(M);
            vn(j,:) = Cn + h*(1-Y)*an(j);
            un(j,:) = Bn + h*h*(1/2-B)*an(j);
        end
        a(i,:) = an(j,:);
        u(i,:) = un(j,:);
        v(i,:) = vn(j,:);
    end
end



%. GRAPHICAL OUTPUT OF RESULTS
xmin = 2*min(min(u)); xmax = 2*max(max(u));
ymin = 0; ymax = N*LStorey;

%. Plot 2D motion in space as a movie
fig5=figure(5);
kframes = 3; kforget = 1;
numframes = floor(timeSteps/kframes);
Mov(:,numframes) = getframe(fig5);

for k=1:numframes
    grid on; axis ([xmin xmax ymin ymax]); axis square; hold on;
    xlabel('x'); ylabel('y'); title('Trajectory');
    kk = k - kforget; kk = (kk + abs(kk))/2;
    istart = 1 + kframes*kk; iend = kframes*k;
    xao = 0; yao = 0;
    for i=1:N
        %....... Plot the masses
        xa = u(iend,i); ya = i*LStorey;
        lnx = [xao,xa]; lny = [yao,ya];
        plot(lnx,lny,'-','Color','b','Linewidth',2);
        plot(xa,ya,'--rs','Color','b','Linewidth',2,...
            'MarkerFaceColor','k',...
            'MarkerEdgeColor','k',...
            'MarkerSize',9);
        xao = xa; yao = ya;
    end % for i
    
    %..... Store frame for movie
    Mov(k) = getframe(fig5);
    clf;
end % for k

for i = 1:N
    plot (t,u(:,i));
    hold on;
end
