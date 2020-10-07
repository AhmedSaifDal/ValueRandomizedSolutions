% Distributionally-robust assignment problem with
% Wasserstein ambiguity set of radius epsilon and with randomization
clc
clear
I=10; % Number of customers
J=I; % Number of servers (equal to the number of customers in balanced problems)
N=10; % Number of data points
C = [-speye(I*J);speye(I*J)]; % Constraint matrix for the support set
epsn = [0 200 400 600 800 1000 1200 1600 2000 2400 2800 3200 3600 4000 6000 10000 15000 25000]; % epsilon values to be tested
fname = sprintf('DRAP%d%d.xlsx',I,N); % Name of the output file (MS Excel)
for nnn = 1:length(epsn) % For all values of epsilon
eps = epsn(nnn)
for ins = 1:10 % Test ten randomly generated instances
rng(ins) % Set the seed for random numbers equal to the instance number (for reproduceability of results)
mean_cost = 10*ones(I,J)+10*rand(I,J); % Mean assignment cost
max_dev= 0.5+0.5*rand(I,J); % Maximum (relative) deviation from the mean
d = [-mean_cost.*(1-max_dev),mean_cost.*(1+max_dev)]; % RHS of the support set constraints
xi = repmat(mean_cost.*(1-max_dev),1,1,N)+2*repmat(mean_cost.*max_dev,1,1,N).*rand(I,J,N); % Data sample of size N, drawn uniformly and independently from the support set at random
%% Solve the DRO problem with and without randomization
[xr, pr, vr, cpur, iter, xd, vd, cpud] = DRAPr(xi,eps,C,d,1); % Call the function DRAPr to solve the deterministic and the randomized strategy problems
%% Test on out-of-sample data
NN = 100*N; % Number of data points for out-of-sample testing
xitest = repmat(mean_cost.*(1-max_dev),1,1,NN)+2*repmat(mean_cost.*max_dev,1,1,NN).*rand(I,J,NN); % Generate out-of-sample data
% % Deterministic strategy
DROd = permute(sum(sum(repmat(reshape(xd,I,J),1,1,NN).*xitest),2),[3 2 1]); % Compute the OS assignment costs
OSd = mean(permute(sum(sum(repmat(reshape(xd,I,J),1,1,NN).*xitest),2),[3 2 1])); % Average OS cost
% % Randomized strategy
for nn=1:NN
ind=find(rand <= cumsum(pr),1,'first'); % Pick the assignment randomly according to the probabilities in the optimal randomized strategy
DROr(nn) = sum(sum(reshape(xr(:,ind),I,J).*xitest(:,:,nn)),2); % Compute the expected OS assignment costs
end
OSr = mean(DROr); % Average OS cost
Result = [ins, vd, cpud, vr, length(pr),cpur,OSd,OSr]; % Results, to be stored in the excel file
row = sprintf('A%d', ins+2); % Row number
xlswrite(fname,Result,nnn,row); % Write to excel
end
xlswrite(fname,eps,nnn,'A1'); % Write the epsilon value in excel
end
