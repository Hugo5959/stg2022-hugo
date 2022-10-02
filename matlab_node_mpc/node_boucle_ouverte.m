%%Parametres initiaux
% Omega is the target sets to reach
global T;
global P;
global N;
global x;
global input_gauche ;                         
global input_droit;
x = [0; 0; pi/2; 0; 0; 0];
P= Polyhedron('lb',[-.5 -.5],'ub',[.5 .5]); % Square centered on the origin
T = 1;                                    % sample time (s)
N = 5;                                     % Prediction Horizon
input_droit = 5;
input_gauche = 1;
%pyenv('Version','/usr/bin/python3.9') 


rosshutdown; %suppression du dernier node specifique matlab
rosinit; %création du node spécifique matlab, adresse ip du master
%%Creation du node MPC a partir de l'adresse ip du master

node1 = ros.Node('/MPCh','10.89.2.241'); %creation du node qui se connecte au ros master communication avec un ROS qui tourne sur la même machine ici machine


% Create ROS subscribers and publishers
sub = ros.Subscriber(node1,'/vessel_state','geometry_msgs/Twist',...
    @callback,...
    'DataFormat','struct'); %suscriber en callback  
sub1 = ros.Subscriber(node1,'/data_sensors','geometry_msgs/Twist',...
    @callback1,...
    'DataFormat','struct'); %suscriber en callback  

pub_left = ros.Publisher(node1,'/wamv/thrusters/left_thrust_cmd','std_msgs/Float32'); %publisher de la commande
pub_right = ros.Publisher(node1,'/wamv/thrusters/right_thrust_cmd','std_msgs/Float32');
pub_sensors = ros.Publisher(node1, '/input_vessel','geometry_msgs/Point');

msg_left = rosmessage(pub_left);
msg_right = rosmessage(pub_right);
input = rosmessage(pub_sensors);

r = rosrate(1); %frequence du pub 20Hz
reset(r);
while (25)
    fprintf('Node is alive..  \n');
    msg_left.Data = input_gauche;
    msg_right.Data = input_droit;
    input.X = input_gauche;
    input.Y = input_droit;
    send(pub_left,msg_left);
    send(pub_right,msg_right);
    send(pub_sensors,input);
    %if isInside(Omega,[x(1);x(2)]) == 1
    %    rosshutdown;
    waitfor(r); 
    
end

%Suscriber callback function


function callback(~,state)
P= Polyhedron('lb',[-.5 -.5],'ub',[.5 .5]);
Omega = [20 20]'+P; % Target Set: Square centered on the point (3,3)
global T;
global x;
N = 5;                                   % Prediction Horizon
global input_gauche;
global input_droit;

fprintf("x : (%f)\n",state.Linear.X);

x =[state.Linear.X , state.Linear.Y , state.Linear.Z , state.Angular.X , state.Angular.Y , state.Angular.Z];
x = x.';
kappa_MPC = kappa(x,Omega,T,N);  %reprend le code de simulation
input_gauche = kappa_MPC(1);
input_droit = kappa_MPC(2);
end


function callback1(~,state)


fprintf("x : (%f,%f,%f,%f,%f,%f)\n",state.Linear.X,state.Linear.Y,state.Linear.Z,state.Angular.X,state.Angular.Y,state.Angular.Z)


end
