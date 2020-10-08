%% DRAPr solves the RANDOMIZED distributionally-robust assignment problem with a Wasserstein distributional ambiguity set and a polyhedral support set
% I is the number of customers
% J is the number of servers
% N is the number of data points
% L is the number of constraints in the support set matrix
% xi is the data points matrix of size (I*J)*N, representing the assignment
% cost
% eps is the radius of the Wasserstein ball, a positive scaler
% C is the support  matrix of size L*I
% d is the RHS of the support set, a column vector of length L
% display is an optional boolean input, if 'true' displays the iterations
% xr is a matrix that contains the solutions that have positive probability, each is a column vector of length I*J
% pr is the optimal solution, i.e., probabilities of the solutions
% comprising the optimal randomized strategy
% vr is the optimal value of the randomized problem
% cpur is the computational time of the randomized problem is seconds
% iter is the number of iterations required to solve the randomized problem
% xd is the solution of the deterministic problem, a column vector of length I*J
% vd is the optimal value of the deterministic problem
% cpud is the computational time of the deterministic problem is seconds

function [xr, pr, vr, cpur, iter, xd, vd, cpud] = DRAPr(xi,eps,C,d,display)
[I, J, N] = size(xi); % Extract problem parameters from the data set
L = numel(d); % Extract the parameter L from the support set definition
if  min(size(C)==[L I*J])==0 % Use the support set constraint matrix C to check the dimensional consistency of the inputs
    error('Incorrect input size') 
end
if (~exist('display', 'var')) % By default, display is off
        display = false;
end

% Some matrices to be used in the models
d1 = repmat(reshape(d,2*I*J,1),1,N)-(C*reshape(xi,I*J,N));
dd = zeros(N,N*L);
for n = 1:N
   dd(n,(n-1)*L+1:n*L)=d1(:,n)'; 
end
xxi = reshape (xi,I*J,N);

% Set Gurobi parameters
params.outputflag = 0;

%% Without Randomization
%%
% Mathematical model: order of variables:  x_ij, gamma_ln, s_n, lam

% Create a Gurobi model object for the determistic problem
clear modeld
modeld.A = [repmat(speye(I),1,J),sparse(I,N*L+N+1);kron(speye(J),ones(1,I)),sparse(J,N*L+N+1);xxi',dd,-speye(N),sparse(N,1);-repmat(speye(I*J),N,1),kron(speye(N),C'),sparse(I*J*N,N),-ones(I*J*N,1) ; repmat(speye(I*J),N,1),-kron(speye(N),C'),sparse(I*J*N,N),-ones(I*J*N,1)];
modeld.rhs = [ones(I+J,1);zeros(N+2*N*I*J,1)];
modeld.obj = [zeros(1,I*J+L*N),ones(1,N)/N,eps];
modeld.sense = [repmat('=',1,I+J),repmat('<',1,N+2*N*I*J)];
modeld.modelsense = 'min';
modeld.vtype = [repmat('I',1,I*J),repmat('C',1,N*L+N+1)];
gurobi_write(modeld, 'mip.lp');

% Solve the deterministic problem
tic
result = gurobi(modeld, params);
cpud = toc; % CPU time for the deterministic problem
sold = result.x; 
vd = result.objval; % Deterministic problem optimal value
xd = sold(1:I*J); % Deteministic problem optimal solution


%% With Randomization
%% Solve the LP-relaxed deterministic problem
% Create a Gurobi model object for the LP-relaxed determistic problem
modeld.vtype = repmat('C',1,I*J+1+N+N*L); % Use the same model of the deterministic problem, just make all variables continuous
gurobi_write(modeld, 'mip.lp');
result = gurobi(modeld, params); % Solve the LP-relaxed deterministic problem
solr = result.x;
vr = result.objval; % Optimal value of the LP-relaxed deterministic problem, also the optimal value of the randomized problem
xrel = solr(1:I*J); % Solution of the LP-relaxed problem

% Initialize parameters for the randomized problem
xnew = xd; % Use the deterministic problem solution to initiate the randomized problem
diff = 100000; % Some arbitrary large number
f = ones(1,I*J); % Initial objective function of the master problem
Aeq = sparse(1,I*J); % Initial equality constraints of the master problem
Aineq = [-speye(I*J);-speye(I*J)]; % Initial inequality constraints of the master problem
ll = zeros(1,I*J); % Lower bound for variables in the master problem
xsols1 = xd; % Initialize the solutions matrix using the deterministic solution
tic
iter = 0; % Start counting the iterations
while diff > 0.02 % If the difference is greater than the tolerance
    iter = iter+1; % Increment the iterations
%% Solve the master problem to find the optimal probabilities of solutions
% order of variables: \theta_ij, p_k
f = [f,0]; % Update the objective function
Aeq = [Aeq,1]; % Update the equality constraints
Aineq = [Aineq,[-xnew;xnew]]; % Update the inequality constraints
ll = [ll,0]; % Update the lower bounds

% Create a Gurobi model object for the master problem
clear modelm
modelm.A = [Aineq;Aeq];
modelm.rhs = [-xrel;xrel;1];
modelm.obj = f;
modelm.sense = [repmat('<',1,2*I*J),'='];
modelm.modelsense = 'min';
modelm.vtype = repmat('C',1,length(ll));
modelm.lb = ll;
gurobi_write(modelm, 'mip.lp');

% Solve the maste problem
resultm = gurobi(modelm, params);
solm = resultm.x;
diff = resultm.objval; % The difference
dual = -resultm.pi; % Extract the dual variables

if display == true
    fprintf('iter=%2g, diff=%2.4f\n', iter, diff); % Display outputs
end
if diff <= 0.1 % Check if the difference is small enough, if so, stop
    break
end

%% Solve a subproblem to generate a new solution
psi1 = dual(1:I*J);
psi2 = dual(I*J+1:2*I*J);

% Create a Gurobi model object for the subproblem
clear models
models.A = [repmat(speye(I),1,J);kron(speye(J),ones(1,I))];
models.rhs = ones(I+J,1);
models.obj = (-psi1+psi2)';
models.sense = repmat('<',1,I+J);
models.modelsense = 'min';
models.vtype = repmat('C',1,I*J);
gurobi_write(models, 'mip.lp');

% Solve the subproblem
results = gurobi(models, params);
xnew = results.x; % Extract the new solution
xsols1 = [xsols1,xnew]; % Add the new solution to the solutions matrix
end
cpur = toc; % CPU time for the randomized problem
p1 = solm(I*J+1:end)'; % Probabilties of the solutions in the optimal strategy
xr = xsols1(:,p1~=0); % Filter the solutions with nonzero probability only
pr = p1(p1~=0)'; % Probabilities of the filtered solutions
end