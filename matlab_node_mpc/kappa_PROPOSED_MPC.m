function kappa = kappa_PROPOSED_MPC(x0,Omega1,Omega2,Omega3,T,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL OF THE DYNAMICAL OF THE VESSEL
b = 1.83/2; m11 = 180.05; m22 = 179; m33 = 248.2;
X_u = 2; Y_v = -0.5; N_r = -0.65; % parameters of the model
f  = @(x,u)[x(4).*cos(x(3))-x(5).*sin(x(3));
        x(4).*sin(x(3))+x(5).*cos(x(3));
        x(6);
        (u(1)+u(2))/m11+m22/m11*x(5).*x(6)+X_u/m11*x(4);
        -m11/m22*x(4).*x(6)+Y_v/m22.*x(5);
        b*(u(1)-u(2))/m33+(m22-m11)/m33*x(4).*x(5)+N_r/m33.*x(6)]; % model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets
% Step is the lenght of set Omega1
Step = 1;
Poly=Polyhedron('lb',[-Step/6, -Step/6],'ub',[Step/6, Step/6]);
%
OmegaIn = Omega1-Poly;
AoIn = OmegaIn.A;
boIn = OmegaIn.b;        % x in Omega iff Ao*x<= bo
% DEFINE Omegax1:
Ao1 = Omega1.A;
bo1 = Omega1.b;        % x in Omega iff Ao*x<= bo
% DEFINE Omegax2:
Ao2 = Omega2.A;
bo2 = Omega2.b;        % x in Omega iff Ao*x<= bo
% DEFINE Omegax3:
Ao3 = Omega3.A;
bo3 = Omega3.b;        % x in Omega iff Ao*x<= bo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prediction horizon for Omega1, Omega2, Omega3
N1 = floor(0.7*N);
N2 = floor(0.3*N);
N3 = N-N1-N2;


% CONTROL PARAMETERS
F1_max = 250; F2_max = 250; 
% Bounds of inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OPTIMIZATION

opti = casadi.Opti(); % Optimization problem in CASADI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECISION VARIABLES
X = opti.variable(6,N+1); % predicted states
x = X(1,:); y = X(2,:); psi = X(3,:); u = X(4,:); v = X(5,:); r = X(6,:); 
Xa = opti.variable(2,N+1); % artificial variables 
U = opti.variable(2,N);   % control variables
F1 = U(1,:); F2 = U(2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COST FUNCTION

V = 0; % cost function
Q = zeros(2,2); Q(1,1) = 30; Q(2,2) = 30; 
R = zeros(2,2); R(1,1) = .005; R(2,2) = .005;
for k = 1:N
    V = V + ([x(k)-Xa(1,k);y(k)-Xa(2,k)])'*Q*([x(k)-Xa(1,k);y(k)-Xa(2,k)])+ U(:,k)'*R*U(:,k);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRAINTS
opti.subject_to(0<=F1<=F1_max);           % u in U
opti.subject_to(0<=F2<=F2_max);           % u in U
%opti.subject_to(0<=psi<=2*pi);           % pi in [0,2pi]
%for k=1:N1
%    opti.subject_to(AoIn*[Xa(1,k);Xa(2,k)]<=boIn); % xak in OmegaIn   
%end
for k=1:N1
    opti.subject_to(Ao1*[Xa(1,k);Xa(2,k)]<=bo1); % xak in Omega1   
end
for k=N1+1:N1+N2
    opti.subject_to(Ao2*[Xa(1,k);Xa(2,k)]<=bo2); % xak in Omega2   
end
for k=N1+N2+1:N1+N2+N3+1
    opti.subject_to(Ao3*[Xa(1,k);Xa(2,k)]<=bo3); % xak in Omega3   
end
opti.subject_to(X(:,1)==x0);   % x0 = x 
%opti.subject_to(Ao1*[x(N1+1);y(N1+1)]<=bo1); % x_N+1 in Omega
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