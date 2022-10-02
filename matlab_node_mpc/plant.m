function x = plant(x0,u0)


% Parameters 
b = 1.83/2; m11 = 180.05; m22 = 179; m33 = 248.2;
X_u = 5; Y_v = -0.5; N_r = -0.65; % parameters of the model
%b = 1.83/2; m11 = 180.05; m22 = 179; m33 = 248.2;
%X_u = 51.3; Y_v = 40; N_r = 400; % parameters of the model
F1_max = 100; F2_max = 100; % Bounds on inputs
T = 1; % sample time (s)

% dynamic model
f  = @(x,u)[x(4).*cos(x(3))-x(5).*sin(x(3));
        x(4).*sin(x(3))+x(5).*cos(x(3));
        x(6);
        (u(1)+u(2))/m11+m22/m11*x(5).*x(6)+X_u/m11*x(4);
        -m11/m22*x(4).*x(6)+Y_v/m22.*x(5);
        b*(u(1)-u(2))/m33+(m22-m11)/m33*x(4).*x(5)+N_r/m33.*x(6)];


% Runge-Kutta 4 integration
   k1 = f(x0,        u0);
   k2 = f(x0+T/2*k1, u0);
   k3 = f(x0+T/2*k2, u0);
   k4 = f(x0+T*k3,   u0);
   x = x0 + T/6*(k1+2*k2+2*k3+k4)+0.01*rand(1,6)';
   