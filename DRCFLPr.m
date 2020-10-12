%% DRCFLPr solves the RANDOMIZED distributionally-robust capacitated facility location problem with a Wasserstein distributional ambiguity set and a polyhedral support set
% I is the number of demand points
% J is the number of potential facility locations
% N is the number of data points
% L is the number of constraints in the support matrix
% xi is the data points matrix of size I*N, representing the demands
% eps is the radius of the Wasserstein ball, a positive scaler
% f is the facility setup cost column vector of length J 
% c is the unit shipping cost matrix of size I*J
% C is the support  matrix of size L*I
% d is the RHS of the support set, a column vector of length L
% x is an optional initial solution, a column vector of size J, if no solution is provided, initial solution is to open all facilities 
% xr is a matrix that contains the solutions that have positive probability, each is a column vector of length J
% pr is the optimal randomized strategy, i.e., probabilities of the solutions
% vr is the optimal value of the randomized problem
% cpur is the computational time of the randomized problem is seconds
% cpum1 is the cumulative computational time for solving the primal master problems
% cpus1 is the cumulative computational time for solving the primal subproblems
% cpum2 is the cumulative computational time for solving the dual master problems
% cpus2 is the cumulative computational time for solving the dual subproblems
% iter1 is the number of primal problem iterations
% iter2 is the number of dual problem iterations
% iterr is the number of algorithm iterations

function [xr, pr, vr, cpur, cpum1, cpus1, cpum2, cpus2, itert1, itert2, iterr] = DRCFLPr(xi,eps,f,c,V,C,d,x)
% Extract problem parameters from the data set
J = length(f);
[I, N] = size(xi);
L = length(d);
if min(size(c)==[I, J])==0 || min(size(V)==[J, 1])==0 || min(size(C)==[L I])==0 % Use the remaining inputs to check the dimensional consistency of the inputs
    error('Incorrect input size')
end
if (~exist('x', 'var')) % If no initial solution is provided, open all facilities
        x = ones(J,1);
end
M = 200; % Big M
Tol = 0.02; % Tolerance

% Set Gurobi parameters
params.outputflag = 0;
 
% Some matrices to be used in the models
dd = repmat(d,1,N)-C*xi;
if min(min(dd),[],2)<0
error('Data points outside the support set')
end
cc = reshape(c,I*J,1);
%% Solve the DRO problem with randomization
%%
S = zeros(1,N); % initial S-bar
G = cell(1,N); % N sets of evaluated costs
K = 1; % Number of feasible solutions
H = zeros(1,N); 
Ver = cell(1,N); % N sets of vertices, one for each sample points
lam = 1; % Initial lambda value
p = 1; % Intial probability distribution of the randomized strategy, (just the deterministic solution)
UB = 1e10*ones; % Initial upper bound
LB = -1e10*ones; % Initial lower bound

% Create a Gurobi model object for the second-stage problem
clear modelev
modelev.A = [kron(speye(J),ones(1,I));repmat(speye(I),1,J)];
modelev.sense = [repmat('<',1,J),repmat('=',1,I)];
modelev.vtype = repmat('C',1,I*J);
modelev.modelsense = 'min';
modelev.obj = cc';
modelev.lb = zeros(1,I*J);
   
iterr = 0; % Initialize iteration count
cpum1=0; % Initialize primal master problem computational time
cpus1=0; % Initialize primal subproblem computatuonal time
cpum2=0; % Initialize dual master problem computational time
cpus2=0; % Initialize dual subproblem computational time
itert1=0; % Initialize iteration count for the primal problem
itert2=0; % Initialize iteration count for the dual problem
tic % Start the clock
fprintf('\n      Randomized Strategy \n\n')
while UB - LB >Tol % While the gap is greater than the optimality tolerance 
iterr = iterr+1; % Increment the iterations
UB1 = UB; % Set UB1 equal to UB
LB1 = -1e10*ones; % LB1
iter1 = 0; % Initialize primal iteration count
%% Solve the primal subproblem to generate new vertices
%%
% Order of variables: xi_i, gamma, delta_i, nu_ik, mu_jk, alpha_l, beta, phi^+_i,
% phi^-_i, B1_l, B2, B3_i, B4_i, B5, B6_i
% Create a Gurobi model object for the primal subproblem 
clear modelspp;
modelspp.A = [[sparse(I*J*K,I+1+I),kron(speye(K),repmat(speye(I),J,1)),-kron(speye(J*K),ones(I,1)),sparse(I*J*K,L+1+I+I+L+1+I+I+1+I)];... % 95
           [C,sparse(L,1+I+I*K+J*K+L+1+I+I+L+1+I+I+1+I)];... % 96
           [-C, sparse(L,1+I+I*K+J*K+L+1+I+I),M*speye(L),sparse(L,1+I+I+1+I)];... % 96
           [sparse(L,I+1+I+I*K+J*K),speye(L),sparse(L,1+I+I),-M*speye(L),sparse(L,1+I+I+1+I)];... % 97       
           [sparse(1,I),-1,ones(1,I),sparse(1,I*K+J*K+L+1+I+I+L+1+I+I+1+I)];... % 98
           [sparse(1,I),1,-ones(1,I),sparse(1,I*K+J*K+L+1+I+I+L),M,sparse(1,I+I+1+I)];... % 98
           [sparse(1,I+1+I+I*K+J*K+L),1,sparse(1,I+I+L),-M,sparse(1,I+I+1+I)];... % 99
           [sparse(I,I+1+I+I*K+J*K+L+1),speye(I),sparse(I,I+L+1),-M*speye(I),sparse(I,I+1+I)];... %101
           [sparse(I,I+1+I+I*K+J*K+L+1+I),speye(I),sparse(I,L+1+I),-M*speye(I,I),sparse(I,1+I)];... % 103
           [sparse(1,I+1+I+I*K+J*K+L),1,sparse(1,I+I+L+1+I+I+1+I)];... % 105
           [sparse(1,I+1+I+I*K+J*K+L),-1,sparse(1,I+I+L+1+I+I),M,sparse(1,I)];... % 105 
           [sparse(1,I),1,sparse(1,I+I*K+J*K+L+1+I+I+L+1+I+I),-M,sparse(1,I)];... % 106
           [sparse(I,I+1+I+I*K+J*K+L),-ones(I,1), speye(I), speye(I),sparse(I,L+1+I+I+1+I)];... % 107
           [sparse(I,I+1+I+I*K+J*K+L), ones(I,1),-speye(I),-speye(I),sparse(I,L+1+I+I+1),M*speye(I)];... % 107
           [sparse(I,I+1),speye(I),sparse(I,I*K+J*K+L+1+I+I+L+1+I+I+1),-M*speye(I)];...% 108
           [speye(I),sparse(I,1),-speye(I),sparse(I,I*K+J*K+L+1+I+I+L+1+I+I+1+I)];... % 100
           [-speye(I),sparse(I,1),speye(I),sparse(I,I*K+J*K+L+1+I+I+L+1),M*speye(I),sparse(I,I+1+I)];... % 100 
           [-speye(I),sparse(I,1),-speye(I),sparse(I,I*K+J*K+L+1+I+I+L+1+I+I+1+I)];... % 102
           [speye(I),sparse(I,1),speye(I),sparse(I,I*K+J*K+L+1+I+I+L+1+I),M*speye(I),sparse(I,1+I)];... % 102
           [sparse(I,I+1+I),-repmat(speye(I),1,K),sparse(I,J*K),C',sparse(I,1),speye(I),-speye(I),sparse(I,L+1+I+I+1+I)]]; %104
 modelspp.vtype = [repmat('C',1,I+1+I+I*K+J*K+L+1+I+I),repmat('B',1,L+1+I+I+1+I)];
 modelspp.sense = [repmat('<',1,I*J*K+3*L+9*I+6),repmat('=',1,I)];
 modelspp.modelsense = 'max';
 modelspp.lb = [zeros(1,I+1+I),-1e10*ones(1,I*K),zeros(1,J*K+L+1+I+I+L+1+I+I+1+I)];
 fprintf('Primal%d\n',iterr)
while UB1-LB1>Tol/2
iter1 = iter1+1; % Increment the primal iterations
cpubs1 = toc; % beginning of SP1 time
 for n=1:N
modelspp.rhs = [repmat(cc,K,1).*kron(p',ones(I*J,1));d;M*ones(L,1)-d;zeros(L,1);0;M;0;zeros(I,1);zeros(I,1);lam;M-lam;0;zeros(I,1);M*ones(I,1);zeros(I,1);xi(:,n);M*ones(I,1)-xi(:,n);-xi(:,n);M*ones(I,1)+xi(:,n);zeros(I,1)];
modelspp.obj = [zeros(1,I),0,zeros(1,I+I*K),-repmat(V',1,K).*reshape(x,1,J*K),d',0,xi(:,n)',-xi(:,n)',zeros(1,L+1+I+I+1+I)];
gurobi_write(modelspp, 'mip1.lp');

% Solve the primal subprobem for each n
result = gurobi(modelspp, params);
sol = result.x;
S(n) = result.objval;
ver = sol(1:I+1); % Extract the vertix
if iterr+iter1 > 2 % If this is not the first iteration
if min(max(abs(Ver{n}-repmat(ver,1,H(n)))))>0.001 % Check if the vertix is new
Ver{n} = [Ver{n},ver]; % Add the new vertix to the set of vertices
H(n) = H(n)+1; % Update the number of vertices
% Evaluate the 2nd stage objective value for the new vertices
for k = 1:K
modelev.rhs = [V.*x(:,k);ver(1:I)];
gurobi_write(modelev, 'eval.lp');
% Solve the evaluation problem
result = gurobi(modelev, params);
G{n}(H(n),k) = result.objval; % Evaluated value   
end
end
else % If this is the first iteration
Ver{n} = ver; % Add the new vertix
H(n) = 1; % The number of vertices is one
% Evaluate the 2nd stage objective value for the new vertices
for k = 1:K    
modelev.rhs = [V.*x(:,k);ver(1:I)];
gurobi_write(modelev, 'eval.lp');
% Solve the evaluation problem
result = gurobi(modelev, params);
G{n}(H(n),k) = result.objval; % Evaluated value
end     
end
end
cpues1 = toc; %end of SP1 time
cpus1 = cpus1 + cpues1 - cpubs1; % Update the computational time for the primal subproblem
UB1 = min(UB1,sum(S)/N+p*(f*x)'+eps*lam); % Update UB1
HH = sum(H); % Total number of vertices

%% Solve the restricted primal master problem to update the value of lambda and p
%%
% Order of variables: p_k, lambda, s_n
% Create a Gurobi model object for the primal master problem 
clear modelmpp
modelmpp.obj = [f*x,eps,ones(1,N)/N];
Aineqm = [];
for n= 1:N
   Aineqm = [Aineqm;G{n},-Ver{n}(end,:)',-repmat(sparse(1,n,1,1,N),H(n),1)]; % Create the inequality constraint matrix
end
modelmpp.A = [Aineqm;ones(1,K),sparse(1,1+N)];
modelmpp.rhs = [zeros(HH,1);1];
modelmpp.lb = zeros(1,K+1+N);
modelmpp.ub = [ones(1,K),M,1e10*ones(1,N)];
modelmpp.modelsense = 'min';
modelmpp.sense = [repmat('<',1,HH),'='];
modelmpp.vtype = repmat('C',1,K+1+N);
gurobi_write(modelmpp, 'mip2.lp');

% Solve the primal master problem
cpubm1 = toc;
result = gurobi(modelmpp, params); 
cpuem1 = toc;
cpum1 = cpum1 + cpuem1 - cpubm1; % Update the computational time for the primal master problem
solm = result.x;
LB1 = result.objval; % Update LB1
dual = -result.pi; % Extract the dual variables
p = solm(1:K)'; % Probabilities
lam = solm(K+1); % Update lambda
cpu1 = toc; % Computational time so far
fprintf('iter1=%2g, LB1=%2.3f, UB1=%2.3f, H=%2g, cpu=%2.3f\n', iter1, LB1, UB1, HH, cpu1) % Print outputs
end
itert1 = itert1 + iter1; % Update the number of primal problem iterations
UB = UB1; % Update overall UB
if UB - LB <Tol % Check if the gap is less than the tolerance, if so, print outputs and terminate
    fprintf('iter1=%2g, LB1=%2.3f, UB1=%2.3f, H=%2g, cpu=%2.3f\n', iter1, LB1, UB1, HH, cpu1)
    break
end

%% Solve the dual subproblem to generate new feasible solutions
%%
fprintf('Dual%d\n',iterr)
UB2 = 1e6*ones; % UB2
LB2 = LB; % Set LB2 equal to LB
%q = dual.ineqlin;
q = dual(1:HH); 
iter2 = 0; % Initialize dual iteration count
while UB2-LB2>Tol/2
iter2=iter2+1; % Increment the primal iterations

% Solve the subproblem to generate a new x
% Order of variables: x_j, z_ijh, w_l
ver1 = cell2mat(Ver);
% Create a Gurobi model object for the dual subproblem
 clear modelspd
 modelspd.A = [sparse(I*HH,J),-kron(speye(HH),repmat(speye(I),1,J)),sparse(I*HH,L);repmat(-diag(V),HH,1),kron(speye(J*HH),ones(1,I)), sparse(J*HH,L);-V',sparse(1,I*J*HH),d';sparse(I,J+I*J*HH),C'];
 modelspd.rhs = [-reshape(ver1(1:I,:),I*HH,1);zeros(J*HH,1);0;ones(I,1)];
 modelspd.vtype = [repmat('B',1,J),repmat('C',1,I*J*HH+L)];
 modelspd.modelsense = 'min';
 modelspd.sense = [repmat('<',1,I*HH+J*HH+1),repmat('=',1,I)];
 modelspd.obj = [f,repmat(cc',1,HH).*kron(q',ones(1,I*J)), zeros(1,L)];
 gurobi_write(modelspd, 'mip3.lp');
% Solve the dual subproblem
 cpubs2 = toc;
 result = gurobi(modelspd, params);
 cpues2 = toc;
 cpus2 = cpus2 + cpues2 - cpubs2; % Update the computational time for the dual subproblem
 xs = result.x;
 vs = result.objval;
xnew = xs(1:J); % New feasible solution
if min(max(abs(x-repmat(xnew,1,K))),[],2)>0.001 % Vheck if it is really a new solution
     x = [x,xnew]; % Add it to the list of solutions
     K = K+1; % Update the number of solutions
 % Evaluate the 2nd stage objective value for the new solution
 for n = 1:N
 for h = 1:H(n)
modelev.rhs = [V.*xnew;Ver{n}(1:I,h)];
gurobi_write(modelev, 'eval.lp');
result = gurobi(modelev, params); 
% Solve the evaluation problem
     G{n}(h,K) = result.objval;
 end    
 end
end
LB2 = max(vs,LB2); % Update LP2
 
%% Solve the restricted dual master problem to update the dual variables q
%%
% Order of variables: w, q_nh, psi+_in, psi-_in
% Create a Gurobi model object for the dual master problem
clear modelmpd
modelmpd.A = [sparse(N,1),-Aineqm(:,end-N+1:end)';ones(K,1),-cell2mat(G')';0,ver1(I+1,:)];
 modelmpd.rhs = [ones(N,1)/N;(f*x)';eps];
 modelmpd.vtype = repmat('C',1,1+HH);
 modelmpd.modelsense = 'min';
 modelmpd.sense = [repmat('=',1,N),repmat('<',1,K+1)];
 modelmpd.obj = [-1,zeros(1,HH)];
 gurobi_write(modelmpd, 'mip3.lp');
% Solve the dual master problem
cpubm2 = toc;
result = gurobi(modelmpd, params);
cpuem2 = toc;
cpum2 = cpum2 + cpuem2 - cpubm2; % Update the computational time for the dual master problem
sold = result.x; 
duald = -result.pi; % Extract the dual varialbles
UB2 = sold(1); % Update UB2
q = sold(2:HH+1);
cpu2 = toc; % Computational time so far
fprintf('iter2=%2g, LB2=%2.3f, UB2=%2.3f, K=%2g,  cpu=%2.3f\n', iter2,  LB2, UB2, K, cpu2) % Print outputs
end
itert2 = itert2 + iter2; % Update the number of primal problem iterations
 p = duald(N+1:N+K)'; % Probabilities
 LB = UB2; % Update overall UB
end
xr = x(:,p~=0); % Filter the solutions with nonzero probability only
pr = p(p~=0); % Probabilities of the filtered solutions
vr = LB; % Optimal value of the randomized problems
cpur = toc; % CPU time for the randomized problem
end