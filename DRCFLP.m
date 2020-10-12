%% DRCFLP solves the two-stage distributionally-robust capacitated facility location problem with a Wasserstein distributional ambiguity set and a polyhedral support set
% I is the number of demand points
% J is the number of potential facility locations
% N is the number of data points
% L is the number of constraints in the support matrix
% xi is the data points matrix of size I*N, representing the demands
% eps is the radius of the Wasserstein ball, a positive scaler
% f is the facility setup cost column vector of length J 
% c is the unit shipping cost matrix of size I*J
% V is the capacity of facilities, a column vector of length J
% C is the support  matrix of size L*I
% d is the RHS of the support set, a column vector of length L
% x is an optional initial solution, a column vector of size J, if no solution is provided, initial solution is to open all facilities 
% xd is optimal solution of the deterministic strategy problem, a column vector of length J
% vd is the optimal value of the deterministic strategy problem
% cpud is the computational time for the deterministic strategy problem is seconds
% iterd is the number of iterations for the deterministic strategy problem
% xLR is optimal solution of the LP-relaxed deterministic strategy problem, a column vector of length J
% vLR is the optimal value of the LP-relaxed deterministic strategy problem
% cpuLR is the computational time for the LP-relaxed deterministic strategy problem is seconds
% iterLR is the number of iterations for the LP-relaxed deterministic strategy problem 

function [xd, vd, cpud, iterd, xLR, vLR, cpuLR, iterLR] = DRCFLP(xi,eps,f,c,V,C,d,x)
% Extract problem parameters from the data set
J = length(f); % Number of facilities
[I, N] = size(xi); % Number of demand points and sample points
L = length(d); % Number of constraints defining the support set
if min(size(c)==[I, J])==0 || min(size(V)==[J, 1])==0 || min(size(C)==[L I])==0 % Check if the input dimensions are correct
error('Incorrect input size')
end
M = 200; % Big M
Tol = 0.02; % Tolerance
if (~exist('x', 'var')) % If no initial solution is provided, open all facilities
x = ones(J,1);
end

% Set Gurobi parameters
params.outputflag = 0;

% Some matrices to be used in the models
cc = reshape(c,I*J,1);

%% Solve the deterministic strategy DRO problem
%%
S = zeros(1,N); % initial S-bar
H = zeros(1,N); % Number of vertices for each point
Ver = cell(1,N); % Set of vertices for each point
lam = 1; % Initial lambda value
UB = 1e10*ones; % Initial upper bound
LB = -1e10*ones; % Initial lower bound
iterd = 0; %Initialize primal iteration count

fprintf('\n      Deterministic Strategy \n\n')
%% Solve the subproblem to generate new vertices
%%
% Order of variables: \xi_i, \gamma, delta_i, z_ij, \nu_i, \mu_j, B1_j, B2_ij
% Create a Gurobi model object for the subproblem 
clear modelsp;
modelsp.A = [[sparse(I*J,I+1+I),repmat(speye(I),J,1),-kron(speye(J),ones(I,1)),sparse(I*J,L+1+I+I+L+1+I+I+1+I)];... % 95
           [C,sparse(L,1+I+I+J+L+1+I+I+L+1+I+I+1+I)];... 
           [-C, sparse(L,1+I+I+J+L+1+I+I),M*speye(L),sparse(L,1+I+I+1+I)];...
           [sparse(L,I+1+I+I+J),speye(L),sparse(L,1+I+I),-M*speye(L),sparse(L,1+I+I+1+I)];...       
           [sparse(1,I),-1,ones(1,I),sparse(1,I+J+L+1+I+I+L+1+I+I+1+I)];... 
           [sparse(1,I),1,-ones(1,I),sparse(1,I+J+L+1+I+I+L),M,sparse(1,I+I+1+I)];... 
           [sparse(1,I+1+I+I+J+L),1,sparse(1,I+I+L),-M,sparse(1,I+I+1+I)];... 
           [sparse(I,I+1+I+I+J+L+1),speye(I),sparse(I,I+L+1),-M*speye(I),sparse(I,I+1+I)];... 
           [sparse(I,I+1+I+I+J+L+1+I),speye(I),sparse(I,L+1+I),-M*speye(I,I),sparse(I,1+I)];... 
           [sparse(1,I+1+I+I+J+L),1,sparse(1,I+I+L+1+I+I+1+I)];... 
           [sparse(1,I+1+I+I+J+L),-1,sparse(1,I+I+L+1+I+I),M,sparse(1,I)];... 
           [sparse(1,I),1,sparse(1,I+I+J+L+1+I+I+L+1+I+I),-M,sparse(1,I)];...
           [sparse(I,I+1+I+I+J+L),-ones(I,1), speye(I), speye(I),sparse(I,L+1+I+I+1+I)];...
           [sparse(I,I+1+I+I+J+L), ones(I,1),-speye(I),-speye(I),sparse(I,L+1+I+I+1),M*speye(I)];...
           [sparse(I,I+1),speye(I),sparse(I,I+J+L+1+I+I+L+1+I+I+1),-M*speye(I)];...
           [speye(I),sparse(I,1),-speye(I),sparse(I,I+J+L+1+I+I+L+1+I+I+1+I)];... 
           [-speye(I),sparse(I,1),speye(I),sparse(I,I+J+L+1+I+I+L+1),M*speye(I),sparse(I,I+1+I)];...
           [-speye(I),sparse(I,1),-speye(I),sparse(I,I+J+L+1+I+I+L+1+I+I+1+I)];... 
           [speye(I),sparse(I,1),speye(I),sparse(I,I+J+L+1+I+I+L+1+I),M*speye(I),sparse(I,1+I)];...
           [sparse(I,I+1+I),-speye(I),sparse(I,J),C',sparse(I,1),speye(I),-speye(I),sparse(I,L+1+I+I+1+I)]];
modelsp.vtype = [repmat('C',1,I+1+I+I+J+L+1+I+I),repmat('B',1,L+1+I+I+1+I)];
modelsp.sense = [repmat('<',1,I*J+3*L+9*I+6),repmat('=',1,I)];
modelsp.modelsense = 'min';

tic % Start the clock
while UB - LB >Tol % While the gap is greater than the optimality tolerance
iterd = iterd+1; % Increment the iterations
 for n=1:N
modelsp.rhs = [cc;d;M*ones(L,1)-d;zeros(L,1);0;M;0;zeros(I,1);zeros(I,1);lam;M-lam;0;zeros(I,1);M*ones(I,1);zeros(I,1);xi(:,n);M*ones(I,1)-xi(:,n);-xi(:,n);M*ones(I,1)+xi(:,n);zeros(I,1)];
modelsp.obj = [zeros(1,I),0,zeros(1,I+I),(V.*x)',-d',0,-xi(:,n)',xi(:,n)',zeros(1,L+1+I+I+1+I)];
gurobi_write(modelsp, 'mip1.lp');

% Solve the subproblem
result = gurobi(modelsp, params);
sol = result.x;
S(n) = -result.objval;
ver = sol(1:I+1); % Extract the vertix
if iterd > 1 % If this is not the first iteration
    if min(max(abs(Ver{n}-repmat(ver,1,H(n)))))>0.001 % Check if the vertix is new
    Ver{n} = [Ver{n},ver]; % Add the new vertix to the set of vertices
    H(n) = H(n)+1; % Update the number of vertices
    end
else % If this is the first iteration
    Ver{n} = ver; % Add the new vertix
    H(n) = 1; % The number of vertices is one     
end
 end
UB = min(UB,sum(S)/N+f*x+eps*lam); % Update the upper bound
if UB - LB <Tol % Check if the gap is less than the tolerance, if so, print outputs and terminate
    cpud = toc;  % CPU time for the deterministic problem
    fprintf('iter=%2g, LB=%2.3f, UB=%2.3f, H=%2g, cpu=%2.3f \n', iterd, LB, UB, HH, cpud);
    break
end

%% Solve the restricted primal problem to update the value of lambda and x
%%
% Order of variables: x_j, z_ijh, lambda, s_n, w_l
HH = sum(H); % Total number of vertices
ver1 = cell2mat(Ver);
str= zeros(HH,N);
count = 1;
for n= 1:N
  str(count:count+H(n)-1,n)=-1;
  count = count+H(n);
end
% Create a Gurobi model object for the master problem
clear modelm
modelm.A = [sparse(HH,J),kron(speye(HH),cc'),-ver1(I+1,:)',str,sparse(HH,L);repmat(-diag(V),HH,1),kron(speye(J*HH),ones(1,I)),sparse(J*HH,1+N),sparse(J*HH,L);-V',sparse(1,I*J*HH+1+N),d';sparse(I*HH,J),kron(speye(HH),repmat(speye(I),1,J)),sparse(I*HH,1+N+L);sparse(I,J+I*J*HH+1+N),C'];
modelm.obj = [f,zeros(1,I*J*HH),eps,ones(1,N)/N, zeros(1,L)];
modelm.rhs = [zeros(HH+HH*J+1,1);reshape(ver1(1:I,:),I*HH,1);ones(I,1)];
modelm.vtype = [repmat('B',1,J),repmat('C',1,I*J*HH+1+N+L)];
modelm.sense = [repmat('<',1,HH+J*HH+1),repmat('=',1,I*HH+I)];
modelm.modelsense = 'min';
gurobi_write(modelm, 'mip1.lp');

% Solve the master problem
resultm = gurobi(modelm, params);
solm = resultm.x;
LB = resultm.objval; % Update the lower bound
x = solm(1:J); % Update the x value
lam = solm(J+I*J*HH+1); % Update the lambda valie
cpud = toc;  % CPU time for the deterministic problem
fprintf('iter=%2g, LB=%2.3f, UB=%2.3f, H=%2g, cpu=%2.3f \n', iterd, LB, UB, HH, cpud); % Print outputs
end
xd=x; % The deterministic strategy solution
vd=LB; % The deterministic strategy optimal value


%% Solve the LP-relaxed deterministic strategy DRO problem 
%%
S = zeros(1,N); % initial S-bar
H = zeros(1,N); % Number of vertices for each point
Ver = cell(1,N); % Set of vertices for each point
lam = 1; % Initial lambda value
UB = 1e10*ones; % Initial upper bound
LB = -1e10*ones; % Initial lower bound
iterLR = 0; %Initialize primal iteration count

fprintf('\n      LP-relaxed Deterministic Strategy \n\n')
%% Solve the subproblem to generate new vertices
%%
tic % Start the clock
while UB - LB >Tol % While the gap is greater than the optimality tolerance
% Create a Gurobi model object for the subproblem
iterLR = iterLR+1; % Increment the iterations
 for n=1:N
modelsp.rhs = [cc;d;M*ones(L,1)-d;zeros(L,1);0;M;0;zeros(I,1);zeros(I,1);lam;M-lam;0;zeros(I,1);M*ones(I,1);zeros(I,1);xi(:,n);M*ones(I,1)-xi(:,n);-xi(:,n);M*ones(I,1)+xi(:,n);zeros(I,1)];
modelsp.obj = [zeros(1,I),0,zeros(1,I+I),(V.*x)',-d',0,-xi(:,n)',xi(:,n)',zeros(1,L+1+I+I+1+I)];
gurobi_write(modelsp, 'mip1.lp');

% Solve the subproblem
result = gurobi(modelsp, params);
sol = result.x;
S(n) = -result.objval;
ver = sol(1:I+1); % Extract the vertix
if iterLR > 1 % If this is not the first iteration
    if min(max(abs(Ver{n}-repmat(ver,1,H(n)))))>0.001 % Check if the vertix is new
    Ver{n} = [Ver{n},ver]; % Add the new vertix to the set of vertices
    H(n) = H(n)+1; % Update the number of vertices
    end
else % If this is the first iteration
    Ver{n} = ver; % Add the new vertix
    H(n) = 1; % The number of vertices is one     
end
 end
UB = min(UB,sum(S)/N+f*x+eps*lam); % Update the upper bound
if UB - LB <Tol % Check if the gap is less than the tolerance, if so, print outputs and terminate
    cpuLR = toc;  % CPU time for the deterministic problem
    fprintf('iter=%2g, LB=%2.3f, UB=%2.3f, H=%2g, cpu=%2.3f \n', iterLR, LB, UB, HH, cpuLR);
    break
end

%% Solve the restricted primal problem to update the value of lambda and x
%%
% Order of variables: x_j, z_ijh, lambda, s_n, w_l
HH = sum(H); % Total number of vertices
ver1 = cell2mat(Ver);
str= zeros(HH,N);
count = 1;
for n= 1:N
  str(count:count+H(n)-1,n)=-1;
  count = count+H(n);
end
% Create a Gurobi model object for the master problem
clear modelm
modelm.A = [sparse(HH,J),kron(speye(HH),cc'),-ver1(I+1,:)',str,sparse(HH,L);repmat(-diag(V),HH,1),kron(speye(J*HH),ones(1,I)),sparse(J*HH,1+N),sparse(J*HH,L);-V',sparse(1,I*J*HH+1+N),d';sparse(I*HH,J),kron(speye(HH),repmat(speye(I),1,J)),sparse(I*HH,1+N+L);sparse(I,J+I*J*HH+1+N),C'];
modelm.obj = [f,zeros(1,I*J*HH),eps,ones(1,N)/N, zeros(1,L)];
modelm.rhs = [zeros(HH+HH*J+1,1);reshape(ver1(1:I,:),I*HH,1);ones(I,1)];
modelm.vtype = repmat('C',1,J+I*J*HH+1+N+L);
modelm.sense = [repmat('<',1,HH+J*HH+1),repmat('=',1,I*HH+I)];
modelm.modelsense = 'min';
gurobi_write(modelm, 'mip1.lp');

% Solve the master problem
resultm = gurobi(modelm, params);
solm = resultm.x;
LB = resultm.objval; % Update the lower bound
x = solm(1:J); % Update the x value
lam = solm(J+I*J*HH+1); % Update the lambda valie
cpuLR = toc;  % CPU time for the LP-relaxed deterministic problem
fprintf('iter=%2g, LB=%2.3f, UB=%2.3f, H=%2g, cpu=%2.3f \n', iterLR, LB, UB, HH, cpuLR); % Print outputs
end
xLR=x; % The deterministic strategy solution
vLR=LB; % The LP-relaxed deterministic strategy optimal value
end