%% Model Predictive Control (MPC) Problem to Tracking Target Sets
% This code presents an MPC problem to tracking targets sets on the surface with an USV
clear all;clc;
addpath('/home/bot/WorkspaceMatlab/matlab_node_mpc/casadi-linux-matlabR2014b-v3.5.5');
import casadi.*;
set(0,'defaultFigureRenderer','painters');
% Omega is the target sets to reach (target set en 6 dimensions)
P= Polyhedron('lb',[-.5 -.5],'ub',[.5 .5]); % Square centered on the origin
Omega = [40 40]'+P; % Target Set: Square centered on the point (3,3)

%% MPC SIMULATION - 

% INITIALIZATION: INITIAL STATE x0, CONTROL PARAMETER (N), Discretization time (T)
T = 1;                                    % sample time (s)
N = 5;                                   % Prediction Horizon
x0 = [0; 0; 0; 0; 0; 0];               % initial state
x = x0;                                   % current state 
xx = [x0];                                % Sequence of controlled states
uu = [];                                  % Sequence of states
Steps = 0;                                 % Step of simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECEDING HORIZON STRATEGY
while isInside(Omega,[x(1);x(2)]) == 0   % Check if the state x reach the set Omega
       kappa_MPC = kappa(x,Omega,T,N);    % Compute the control law of the MPC
       fprintf("(%f)",kappa_MPC)
       Steps = Steps + 1;
       X = ['time step: ', num2str(Steps)];
       disp(X)
       uu = [uu kappa_MPC];                    % Sequence of control
       x = plant(x,kappa_MPC); 
       fprintf("state : (%f,%f)",x(4),x(5));% Actualize the current state
       xx = [xx x];                            % Sequence of state    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PLOTS of the surface (x,y)
vecr=[0.7;0;0];vecv=[0.7;0.5;0];vecb=[0.7;0.5;0.5];vecy=[0.4;0.3;0]; % Color vector
rob_diam = 0.4;                           % Diameter of vessel
figure(1);hold on
plot(Omega,'Color', vecr, 'Alpha', 0.1,'edgecolor',vecr,'LineWidth',0.1)
plot(xx(1,:)',xx(2,:)'','o--');
quiver(xx(1,:)',xx(2,:)',cos(xx(3,:)'),sin(xx(3,:)'),rob_diam,'r')
title('Controlled Trajectory')
xlabel('x')
ylabel('y')
axis([0 4 0 4])
% CONTROLS AND POSITION
tt=[0:1:length(uu(1,:))-1];t=[0:1:length(xx(1,:))-1];
figure(2)
subplot(411)
stairs(tt,uu(1,:),'k','linewidth',1.5)
ylabel('F1')
subplot(412)
stairs(tt,uu(2,:),'k','linewidth',1.5)
ylabel('F2')
subplot(413)
plot(t,xx(1,:),'o--');
ylabel('x')
subplot(414)
plot(t,xx(2,:),'o--');
ylabel('y')
% VELOCITY
figure(3)
subplot(411)
plot(t,xx(3,:),'o--')
ylabel('\psi')
subplot(412)
plot(t,xx(4,:),'o--');
ylabel('u')
subplot(413)
plot(t,xx(5,:),'o--');
ylabel('v')
subplot(414)
plot(t,xx(6,:),'o--');
ylabel('r')

