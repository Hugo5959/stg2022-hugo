%% Importations

load('Zone.mat');
load('X_Zone_Rebuilt.mat');

addpath('/home/bot/WorkspaceMatlab/matlab_node_mpc/casadi-linux-matlabR2014b-v3.5.5');
import casadi.*;


%% Paramètres globaux

global T;                                 % Sample time pour le MPC (s)
global N;                                 % Prediction Horizon MPC
global x_1;                               % Vecteur d'etat du robot (x,y,yaw,vx,vy,vyaw)
global input_gauche_1 ;                   % Commande moteur gauche
global input_droit_1;                     % Commande moteur droit
%global map;                              % Liste des targets sets (inutile dans le cas où on utilise la carte du lac)
global cnt;                               % Compteur pour la liste des targets sets
global Ox_1;                              % Taux d'oxygène
global t;                                 % Temps (s)
global donx_1;                            % Liste abscisse x
global dony_1;                            % Liste abscisse y
global donox_1;                           % Liste donnée oxygène
global dont;                              % Liste des timestamp  
global Oxi_1;                             % Taux d'oxygène sans modification
global donoxori_1;                        % Liste données oxygène origine
%% Initialisation                                                         

cnt = 1;
x_1 = [0; 0; pi/2; 0; 0; 0];               
t = 0;
T = 1;                                    
N = 5;                                   
Ox_1 = 0;  
Oxi_1 = 0;
input_droit_1 = 0;
input_gauche_1 = 0;
donx_1 = [];
dony_1 = [];
donox_1 = [];
dont = [];
donoxori_1 = [];
%map = grille_parcours1(60,80,50,100,5);  % Choix entre grille_parcours1 (zigzag) ou grille_parcours2 (spirale);

%% ROS Architecture

rosshutdown;                             % Suppression du dernier node specifique matlab
rosinit;                                 % Création du node spécifique matlab


node2 = ros.Node('/MPC2','10.1.160.141');  % Création du node MPC qui se connecte au ros master communication avec un ROS qui tourne sur la même machine ici machine
                                           % Il faut entrer l'adresse IP du master en paramètre
                                        

% Create ROS subscribers and publishers

sub_1 = ros.Subscriber(node2,'/vessel_state_1','geometry_msgs/Twist',...
    @callback,...
    'DataFormat','struct'); % Suscriber en callback de l'état du robot

pub_left_1 = ros.Publisher(node2,'/wamv1/thrusters/left_thrust_cmd','std_msgs/Float32');       % Publisher de la commande gauche
pub_right_1 = ros.Publisher(node2,'/wamv1/thrusters/right_thrust_cmd','std_msgs/Float32');     % Publisher de la commande droite
pub_ox_1 = ros.Publisher(node2,'/taux_ox_1','std_msgs/Float32');                               % Publisher du taux d'oxygène

% Message ROS associée aux publishers

msg_left_1 = rosmessage(pub_left_1);
msg_right_1 = rosmessage(pub_right_1);
msg_ox_1 = rosmessage(pub_ox_1);

r = rosrate(0.7);                     % Fréquence des publishers 1 Hz (fréquence lente conseillée pour éviter que ça plante)
reset(r);

while (1)

    fprintf('Node is alive..  \n');
    
    % Envoi des messages

    msg_left_1.Data = input_gauche_1;
    msg_right_1.Data = input_droit_1;

    send(pub_left_1,msg_left_1);
    send(pub_right_1,msg_right_1);

    for i=1:length(Zone) 

        if isInside(Zone(i),[x_1(2);x_1(1)]) 
            Ox_1 = X_Zone_Rebuilt(i).Oxygene ;
            Oxi_1 = Ox_1 + ((-abs( rem(t,86200)-43200 ))* (2.8/43200))+ 1.4;  % Evolution du taux d'oxygène selon l'heure de la journée (cycle de 24h, soit 86200s)
        end
    end
    donoxori_1(end+1) = Ox_1;
    donox_1(end + 1) = Oxi_1;
    dont(end +1) = t;
    fprintf("ox : (%f)",Oxi_1);
    msg_ox_1.Data = Oxi_1;
    send(pub_ox_1, msg_ox_1);
    t = t + (1/0.7);                      % Incrémentation du temps
    
    waitfor(r); 
    
end

%% Suscriber callback function

function callback(~,state_1)

global T;
global x_1;
global N ;                                   
global input_gauche_1;
global input_droit_1;
global cnt;
global donx_1;
global dony_1;
%global map;
load('Position_sorted.mat');

% Création des targets sets avec Polyhédron

P = Polyhedron('lb',[-0.2 -0.2],'ub',[0.2 0.2]);          % Polyhedron pour le MPC
Q = Polyhedron('lb',[-2.5 -2.5],'ub',[2.5 2.5]);          % Polyhedron pour la détection de zone atteinte

%Omega1 = [map(cnt,1) map(cnt,2)]'+P;
%Omega2 = [map(cnt+1,1) map(cnt+1,2)]'+P;
%Omega3 = [map(cnt+2,1) map(cnt+2,2)]'+P;
%Omegak = [map(cnt,1) map(cnt,2)]'+Q;

Omega1 = [Position_sorted(453-cnt,2) Position_sorted(453-cnt,1)]'+P;
%Omega2 = [Position_sorted(cnt+1,2) Position_sorted(cnt+1,1)]'+P;
%Omega3 = [Position_sorted(cnt+2,2) Position_sorted(cnt+2,1)]'+P;
Omegak = [Position_sorted(453-cnt,2)   Position_sorted(453-cnt,1)]'+Q;

x_1 = [state_1.Linear.X , state_1.Linear.Y , state_1.Linear.Z , state_1.Angular.X , state_1.Angular.Y , state_1.Angular.Z];
donx_1(end + 1) = x_1(1);
dony_1(end + 1) = x_1(2);

if isInside(Omegak,[x_1(1);x_1(2)])                 % Vérification si la target set est atteinte et on passe à la suivante (avec omega k ou 1)
    cnt = cnt + 1;
end

x_1 = x_1.';
%kappa_MPC = kappa_PROPOSED_MPC(x,Omega1,Omega2,Omega3,T,N);         % Reprend le code de l'algo MPC multi -target
kappa_MPC = kappa(x_1,Omega1,T,N);                                   % MPC avec une seule target
input_gauche_1 = kappa_MPC(1);
input_droit_1 = kappa_MPC(2);

if input_droit_1+input_gauche_1 < 0.00001             % Permet d'éviter l'arrêt du robot car il considère une zone inatteignable
    input_gauche_1 = 250;
    input_droit_1 = 0;    
end

fprintf("input : (%f,%f)\n",input_gauche_1,input_droit_1);
fprintf("map = (%f,%f)\n",Position_sorted(453-cnt,1),Position_sorted(453-cnt,2));
fprintf("cnt = (%f)\n", cnt)
fprintf("yaw : (%f)\n",state_1.Linear.Z);
fprintf("x : (%f)\n",state_1.Linear.X);
fprintf("y : (%f)\n",state_1.Linear.Y);

end


