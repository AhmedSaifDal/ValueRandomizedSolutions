%% DRUFLPr solves the RANDOMIZED distributionally-robust uncapacitated facility location problem with a Wasserstein distributional ambiguity set and a polyhedral support set
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
% display is an optional boolean input, if 'true' displays the iterations
% xr, yr are matrices that contain the solutions with positive
% probabilities, each is a column vector of length J, and I*J, respectively
% pr is the optimal randomized strategy, i.e., probabilities of the solutions
% vr is the optimal value of the randomized problem
% cpur is the computational time of the randomized problem is seconds
% iter is the number of iterations required to solve the randomized problem
% xd, yd are the solution of the deterministic problem, column vectors of
% length J and I*J, respectively
% vd is the optimal value of the deterministic problem
% cpud is the computational time of the deterministic problem in seconds
% xLR is the solution of the LP-relaxation deterministic problem
% vLR is the optimal value of the LP-relxed deterministic problem
% cpuLR is the computational time of the LP-relaxed deterministic problem
% cpum is the cumulative computational time for solving the master problems
% cpus is the cumulative computational time for solving the sub-problems

function [xr, yr, pr, vr, cpur, iter, xd, yd, vd, cpud, xLR, vLR, cpuLR, cpum, cpus] = DRUFLPr(xi,eps,f,c,C,d,display)
% Extract problem parameters from the data set
J = length(f); 
[I, N] = size(xi);
L = length(d);
if min(size(c)==[I, J])==0 || min(size(C)==[L I])==0 % Use the support set constraint matrix C to check the dimensional consistency of the inputs
    error('Incorrect input size')
end
if (~exist('display', 'var')) % By default, display is off
        display = false;
end

% Some matrices to be used in the models
d1 = repmat(d,1,N)-(C*xi);
dd = sparse(N,N*L);
for i = 1:N
   dd(i,(i-1)*L+1:i*L)=d1(:,i)'; 
end
cc = zeros(I,I*J);
for j=1:J
   cc(:,(j-1)*I+1:j*I)=diag(c(:,j)); 
end

% Set Gurobi parameters
params.outputflag = 0;

%% Without Randomization
%%
% Create a Gurobi model object for the deterministic problem
% %% Mathematical model: order of variables:  x_j, y_ij, lam_ln, s_n, lam
clear modeld
modeld.A = [sparse(I,J),repmat(speye(I),1,J),sparse(I,N+N*L+1);-kron(speye(J),ones(I,1)),speye(I*J),sparse(I*J,N+N*L+1);sparse(N,J),repmat(reshape(c,1,I*J),N,1).*repmat(xi',1,J),dd,-speye(N),sparse(N,1);sparse(I*N,J), -repmat(cc,N,1),kron(speye(N),C'),sparse(I*N,N),-ones(I*N,1) ; sparse(I*N,J), repmat(cc,N,1),-kron(speye(N),C'),sparse(I*N,N),-ones(I*N,1)];
modeld.rhs = [ones(I,1);zeros(I*J+N+2*N*I,1)];    
modeld.obj = [f,zeros(1,I*J+L*N),ones(1,N)/N,eps];
modeld.sense = [repmat('=',1,I),repmat('<',1,I*J+N+2*N*I)];
modeld.modelsense = 'min';
modeld.vtype = [repmat('B',1,J),repmat('C',1,I*J+N*L+N+1)];
gurobi_write(modeld, 'mip.lp');

% Solve the deterministic problem
tic
result = gurobi(modeld, params);
cpud = toc; % CPU time for the deterministic problem
sold = result.x;
vd = result.objval; % Deterministic problem optimal value
xd = sold(1:J);
yd = reshape(sold(J+1:J+I*J),I,J);
ysol = sold(J+1:J+I*J);
xsol = xd;
xsols = xd;
ysols = ysol;

%% Solve the LP-relaxed Deterministic Problem
modeld.vtype = repmat('C',1,J+I*J+N*L+N+1); % Use the same model of the deterministic problem, just make all variables continuous
gurobi_write(modeld, 'mip.lp');
tic
result = gurobi(modeld, params); % Solve the LP-relaxed deterministic problem
solLR = result.x;
vLR = result.objval; % Optimal value of the LP-relaxed deterministic problem
cpuLR = toc;
xLR = solLR(1:J);

%% With Randomization
% Use the deterministic solution as the starting point
% Order of variables: lambda, s_n, lambda_ln, p_k
% Initialize the master problem matrices and vectors
OFr = [eps, ones(1,N)/N,zeros(1,L*N)]; % The objective function
Aineq1r = [sparse(N,1),-speye(N),dd]; % Inequality constraint 1
Aineq2r = [-ones(I*N,1),sparse(I*N,N),kron(speye(N),C');-ones(I*N,1),sparse(I*N,N),-kron(speye(N),C')]; % Inequality constraint 2
Aeqr = sparse(1,1+N+N*L); % Equality constraint
llr = zeros(1,1+N+N*L); % Lower bounds of variables
UB = 100000; % Initialize the upper bound, an arbitrary large number
LB = 0; % Initialize the lower bound
iter = 0; % Start counting the iterations
cpum = 0; % Start counting the master problems cpu time
cpus = 0; % Start counting the subproblem cpu time
tic % Start the clock
while UB-LB > 0.02 % While the gap is greater than the optimality tolerance
iter = iter+1; % Increment the iterations
% Update the master problem matrices and vectors
OFr = [OFr, sum(f.*xsol')]; % Update the objective function vector
Aineq1r = [Aineq1r, xi'*sum(c.*reshape(ysol,I,J),2)]; % Inequality constraint 1
Aineq2r = [Aineq2r, [-repmat(sum(c.*reshape(ysol,I,J),2),N,1);repmat(sum(c.*reshape(ysol,I,J),2),N,1)]]; % Inequality constraint 2
Aeqr = [Aeqr,1]; % Equality constraint
llr = [llr,0]; % Lower bounds of variables

% Create a Gurobi model object for the master problem
clear modelm
modelm.obj = OFr;
modelm.A = [Aineq1r;Aineq2r;Aeqr];
modelm.rhs = [zeros(N+2*I*N,1);1];
modelm.lb = llr;
modelm.modelsense = 'min';
modelm.sense = [repmat('<',1,N+2*I*N),'='];
modelm.vtype = repmat('C',1,length(llr));
gurobi_write(modelm, 'mip2.lp');

% Solve the maste problem
cpubm = toc;
resultm = gurobi(modelm, params);
cpuem = toc;
cpum = cpum + cpuem - cpubm;
solm = resultm.x;
UB = resultm.objval; % Update the upper bound
lam = -resultm.pi; % Extract the dual variables
if UB-LB<0.001 % Check if the gap is small enough, if so, stop
    break
end

%% Solve a subproblem to generate a new solution
psi1 = reshape(lam(N+1:N+I*N),I,N);
psi2 = reshape(lam(N+N*I+1:N+2*I*N),I,N);

% Create a Gurobi model object for the subproblem    
clear models
models.obj = [f,repmat(sum(xi/N-psi1+psi2,2)',1,J).*reshape(c,1,I*J)];
models.A = [-kron(speye(J),ones(I,1)),speye(I*J);sparse(I,J),repmat(speye(I),1,J)];
models.rhs = [zeros(I*J,1);ones(I,1)];
models.lb = zeros(1,J+I*J);
models.ub = ones(1,J+I*J);
models.modelsense = 'min';
models.sense = [repmat('<',1,I*J),repmat('=',1,I)];
models.vtype = [repmat('I',1,J),repmat('C',1,I*J)];
gurobi_write(models, 'mip2.lp');

% Solve the subproblem
cpubs = toc;
results = gurobi(models, params);
cpues = toc;
cpus = cpus + cpues - cpubs;
sols = results.x;
vs = results.objval;
LB = max(LB,vs); % Update the lower bound
xsol = sols(1:J); % Extract the new solution (x & y)
ysol = sols(J+1:J+I*J);
xsols = [xsols,xsol]; % Add the new solution to the solutions matrix
ysols = [ysols,sols(J+1:J+I*J)];
if display == true
fprintf('iter=%2g, LB=%2.4f, UB=%2.4f \n', iter, LB, UB); % Dispaly outputs
end
end
cpur = toc; % CPU time for the randomized problem
vr = UB; % Optimal value of the randomized problems 
p = solm(1+N+N*L+1:end); % All probabilities
xr = xsols(:,p~=0); % Filter the solutions with nonzero probability only
yr = ysols(:,p~=0);
pr = p(p~=0); % Probabilities of the filtered solutions
end