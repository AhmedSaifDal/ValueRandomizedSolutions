%% This code does the following:
% 1. Generates random instances of the Distributionally-robust capacitated facility location problem (DRCFLP)
% 2. Solves the deterministic strategy DRCFLP, and its linear relaxation, by calling the function DRUFLP
% 3. Solves the randomized strategy DRCFLP by calling the function DRUFLPr
% 4. Tests the solutions obrained on out-of-sample data
% 5. Writes the results to an Excel file

%% Set the problem parameters
I=10; % Number of customers
J=I; % Number of potential facility locations
N=10; % Number of data points
for rat= [3 5] % For these capacity ratios
ff = 10; % Setup costs of facilities
f = ff*ones(1,J);
C = [-speye(I);speye(I)]; % Constraint matrix for the support set
epsn = [0 20 40 60 80 100 120 140 160 180 200 220 240 260 280 300];  % epsilon values to be tested
fname = sprintf('CFLP_%d_%d_%d_%d.xlsx', I, rat, ff, N); % Name of the output file (MS Excel)
for nnn = 1:length(epsn) % For all values of epsilon
eps = epsn(nnn);
for ins = 1:10 % Test ten randomly generated instances
rng(ins) % Set the seed for random numbers equal to the instance number (for reproduceability of results)
mean_demand = 10*ones(I,1)+10*rand(I,1); % Mean demand
max_dev=0.5+0.5*rand(I,1); % Maximum (relative) deviation from the mean
d = [-mean_demand.*(1-max_dev);mean_demand.*(1+max_dev)]; % RHS of the support set constraints
L = length(d);
lat = rand(max(I,J),1); % Latitude of the points
lon = rand(max(I,J),1); % Longitude of the poiunts
c=zeros(I,J);
for i =1:I
for j = 1:J
c(i,j)=sqrt((lat(i)-lat(j))^2+(lon(i)-lon(j))^2); % Euclidean distance between points (shipping cost) 
end
end
xi = repmat(mean_demand,1,N)+repmat(max_dev,1,N).*repmat(mean_demand,1,N).*(2*rand(I,N)-1); % Data sample of size N, drawn uniformly and independently from the support set at random
V = ceil(rat*sum(mean_demand)/J)*ones(J,1); % Capacity of facilities

%% Solve the deterministic  strategy problem and its LP-relaxation
[xd, vd, cpud, iterd, xLR, vLR, cpuLR, iterLR] = DRCFLP(xi,eps,f,c,V,C,d);

%% Solve the randomized strategy problem
[xr, pr, vr, cpur, cpum1, cpus1, cpum2, cpus2, itert1, itert2, iterr] = DRCFLPr(xi,eps,f,c,V,C,d);

%% Generate out-of-sample test points
NN= 100* N; % Number of data points for out-of-sample testing
xitest = repmat(mean_demand,1,NN)+repmat(max_dev,1,NN).*repmat(mean_demand,1,NN).*(2*rand(I,NN)-1); % Generate out-of-sample data
%% Test the deterministic and randomized solutions on out-of-sample points
OSdd=zeros(1, NN);
OSrr=zeros(1, NN);
% Create a Gurobi model object for the out-of-sample test problem
params.outputflag = 0;
clear modeld;
modelo.vtype = repmat('C',1,I*J);
modelo.sense = [repmat('<',1,J), repmat('=',1,I)];
model.modelsense = 'min';
model.lb = zeros(1,I*J);
model.ub = ones(1,I*J);
for nn = 1:NN
% Out-of-sample tests of the deterministic strategy solution
modelo.obj = reshape(c,1,I*J).*repmat(xitest(:,nn)',1,J);
modelo.A = [kron(speye(J),xitest(:,nn)');repmat(speye(I),1,J)] ;
modelo.rhs = [V.*xd; ones(I,1)];
gurobi_write(modelo, 'lp.lp');
% Solve the out-of-sample test problem based on the deterministic strategy
% solution
result = gurobi(modelo, params);
OSVd = result.objval;  
OSdd(nn) = f*xd+OSVd; % Compute the expected OS cost for the deterministic strategy solution
% Out-of-sample tests of the randomized strategy solution
ind=find(rand <= cumsum(pr),1,'first'); % Pick the solution randomly according to the probabilities in the optimal randomized strategy
xrr = xr(:,ind);
modelo.rhs = [V.*xrr; ones(I,1)];
gurobi_write(modelo, 'lp.lp');
% Solve the out-of-sample test problem based on the randomized strategy
% solution
result = gurobi(modelo, params);
OSVr = result.objval;
OSrr(nn) = f*xrr+OSVr; % Compute the expected OS cost for the randomized strategy solution
end
OSd = mean(OSdd); % Average OS cost - deterministic
OSr = mean(OSrr); % Average OS cost - randomized
Result = [ins,vd,cpud,iterd,vr,length(pr),cpur,cpum1, cpus1, cpum2, cpus2, itert1, itert2, iterr,vLR,cpuLR,OSd,OSr]; % Results, to be stored in the excel file
row = sprintf('A%d', ins+2); % Row number
xlswrite(fname,Result,nnn,row); % Write to excel
 end
xlswrite(fname,eps,nnn,'A1'); % Write the epsilon value
end
end

