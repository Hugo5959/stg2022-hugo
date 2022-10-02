%% Control law of the Model Predictive Control (MPC) Problem to Tracking Target Sets
function kappa = kappa(x0,Omega,T,N)
% Kappa is the control inputs
% x0 is the current state
% Omega is the target Set
% T discretization time
% N is the horizon of prediction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL OF THE DYNAMICAL OF THE VESSEL
b = 1.83/2; m11 = 180.05; m22 = 179; m33 = 248.2;
X_u = 2; Y_v = -0.5; N_r = -0.65; % parameters of the model
%b = 1.83/2; m11 = 180.05; m22 = 179; m33 = 248.2;
%X_u = 51.3; Y_v = 40; N_r = 400; % parameters of the model
f  = @(x,u)[x(4).*cos(x(3))-x(5).*sin(x(3));
        x(4).*sin(x(3))+x(5).*cos(x(3));
        x(6);
        (u(1)+u(2))/m11+m22/m11*x(5).*x(6)+X_u/m11*x(4);
        -m11/m22*x(4).*x(6)+Y_v/m22.*x(5);
        b*(u(1)-u(2))/m33+(m22-m11)/m33*x(4).*x(5)+N_r/m33.*x(6)]; 
% x(1)=x,x(2)=y,x(3)=psi,x(4)=u,x(5)=v,x(6)=r States of the system of the USV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE Omega: Omega = {x\in\R^2: Aox <= b0}
Ao = Omega.A;
bo = Omega.b;        % x in Omega iff Ao*x<= bo

% CONTROL PARAMETERS
F1_max = 250; F2_max = 250; % Fi is the force of motor i,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OPTIMIZATION
opti = casadi.Opti(); % Optimization problem in CASADI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECISION VARIABLES
X = opti.variable(6,N+1); % predicted states
x = X(1,:); y = X(2,:); psi = X(3,:); u = X(4,:); v = X(5,:); r = X(6,:); % states
Xa = opti.variable(2,N+1); % artificial variables for auxiliar calculus (distance from a point to a set)
U = opti.variable(2,N);   % control variables F1 and F2
F1 = U(1,:); F2 = U(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COST FUNCTION
V = 0; % cost function
Q = zeros(2,2); Q(1,1) = 30; Q(2,2)= 30; % weight parameters for state
R = zeros(2,2); R(1,1) = 0.005; R(2,2) = 0.005; % weight parameters for control
for k = 1:N
    V = V + ([x(k)-Xa(1,k);y(k)-Xa(2,k)])'*Q*([x(k)-Xa(1,k);y(k)-Xa(2,k)])+ U(:,k)'*R*U(:,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRAINTS
opti.subject_to(0<=F1<=F1_max);           % u in U
opti.subject_to(0<=F2<=F2_max);           % u in U
opti.subject_to(X(:,1)==x0);   % x0 = x  Initial condition
%opti.subject_to(Ao*[x(N+1);y(N+1)]<=bo); % x_N+1 in Omega
for k=1:N
    opti.subject_to(Ao*[Xa(1,k);Xa(2,k)]<=bo); % The artificial state belongs to Omega   
end
% x(j+1) = f(x(j),u(j)) Dynamic of the system restriction:
for k=1:N 
   % Runge-Kutta 4 integration
   k1 = f(X(:,k),         U(:,k));
   k2 = f(X(:,k)+T/2*k1, U(:,k));
   k3 = f(X(:,k)+T/2*k2, U(:,k));
   k4 = f(X(:,k)+T*k3,   U(:,k));
   x_next = X(:,k) + T/6*(k1+2*k2+2*k3+k4);
   opti.subject_to(X(:,k+1)==x_next);      % x_j+1 = f(x_j,u_j)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOLVER
opti.minimize(V); 
opti.solver('ipopt'); % set numerical 
sol = opti.solve();   % actual solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOLUTION
v1=sol.value(F1);
v2=sol.value(F2);
kappa = [v1(1) v2(1)]';