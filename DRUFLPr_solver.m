%% This code does the following:
% 1. Generates random instances of the Distributionally-robust uncapacitated facility location problem (DRUFLP)
% 2. Solves the deterministic and randomized strategy DRUFLP by calling the
% function DRUFLPr
% 3. Tests the solutions obrained on out-of-sample data
% 4. Writes the results to an Excel file

%% Generate random DRUFLP instances
I=100; % Number of customers
J=I; % Number of potential facility locations
N=10; % Number of data points
ff = 10; % Setup costs of facilities
f = ff*ones(1,J);
C = [-speye(I);speye(I)]; % Constraint matrix for the support set
epsn = [0 200 400 600 800 1000 1200 1400 1600 1800 2000]; % epsilon values to be tested
fname = sprintf('DRUFLP_%d_%d_%d.xlsx', I, ff, N); % Name of the output file (MS Excel)
for nnn = 1:length(epsn)  % For all values of epsilon
eps = epsn(nnn);
Result = [];
for ins = 1:10 % Test ten randomly generated instances
rng(ins) % Set the seed for random numbers equal to the instance number (for reproduceability of results)
mean_demand = 10*ones(I,1)+10*rand(I,1); % Mean demand
max_dev=0.5+0.5*rand(I,1); % Maximum (relative) deviation from the mea
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

%% Solve the DRUFLP with and without randomization
[xr, yr, pr, vr, cpur, iter, xd, yd, vd, cpud, xLR, vLR, cpuLR, cpum, cpus] = DRUFLPr(xi,eps,f,c,C,d,1);

%% Test on out-of-sample data
NN = 1000*N; % Number of data points for out-of-sample testing
xitest = repmat(mean_demand,1,NN)+repmat(max_dev,1,NN).*repmat(mean_demand,1,NN).*(2*rand(I,NN)-1); % Generate out-of-sample data
% Deterministic strategy
OSd = f*xd+mean(reshape((c.*yd),1,I*J)*repmat(xitest,J,1)); % Compute the expected OS cost
% Randomized strategy
DROr = zeros(1,NN);
for nn=1:NN
ind=find(rand <= cumsum(pr),1,'first'); % Pick the solution randomly according to the probabilities in the optimal randomized strategy
DROr(nn) = f*xr(:,ind)+reshape(c,1,I*J).*yr(:,ind)'*repmat(xitest(:,nn),J,1); % Compute the expected OS costs
end
OSr = mean(DROr); % Average OS cost

%% Write results to an excel file
Result = [ins,vLR,vd,cpud,iter, vr,length(pr),cpur, cpum, cpus, OSd, OSr]; % Results, to be stored in the excel file
row = sprintf('A%d', ins+2); % Row number
xlswrite(fname,Result,nnn,row); % Write to excel
end
xlswrite(fname,eps,nnn,'A1'); % Write the epsilon value
end
