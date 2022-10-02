%clear all
%close all
%clc
%load Lac_Heron_2021.mat
Step = 5; 
Shape = 1;


% %Trajectoire - Zone Delimitation
% figure (111)
%  hold on
%  grid on
% plot(Coord_E(N_data1+Sample_Between+40:N_data-451),Coord_N(N_data1+Sample_Between+40:N_data-451),'r');
% plot(Coord_E(N_data1+Sample_Between),Coord_N(N_data1+Sample_Between),'b*');
% title('Trajectory')
% xlabel('E')
% ylabel('N')
% %legend({'','','','',''},'Location','southeast')

% Determination of Points of Boundaries (Rectangle)
Min_Coord_E=min(Coord_E(N_data1+Sample_Between+40:N_data-451));  % min from East
Max_Coord_E=max(Coord_E(N_data1+Sample_Between+40:N_data-451));  % max from East
Min_Coord_N=min(Coord_N(N_data1+Sample_Between+40:N_data-451));  % min from North
Max_Coord_N=max(Coord_N(N_data1+Sample_Between+40:N_data-451));  % max from North

% Determination of Points around the explored Zone
% Coordonnate of N according to min E 
Coord_N_Min_E=Coord_N(N_data1+Sample_Between+40+find(Coord_E(N_data1+Sample_Between+40:N_data-451)==Min_Coord_E));
Coord_N_Min_E=Coord_N_Min_E(1);
% Coordonnate of N according to max E
Coord_N_Max_E=Coord_N(N_data1+Sample_Between+40+find(Coord_E(N_data1+Sample_Between+40:N_data-451)==Max_Coord_E));
Coord_N_Max_E=Coord_N_Max_E(1);
% Coordonnate of E according to min N 
Coord_E_Min_N=Coord_E(N_data1+Sample_Between+40+find(Coord_N(N_data1+Sample_Between+40:N_data-451)==Min_Coord_N));
Coord_E_Min_N=Coord_E_Min_N(1);
% Coordonnate of E according to max N
Coord_E_Max_N=Coord_E(N_data1+Sample_Between+40+find(Coord_N(N_data1+Sample_Between+40:N_data-451)==Max_Coord_N));
Coord_E_Max_N=Coord_E_Max_N(end);

% Plot of Extreme points - drawing a Rectangle
% plot(Min_Coord_E,Min_Coord_N,'k*');
% plot(Min_Coord_E,Max_Coord_N,'k*');
% plot(Max_Coord_E,Min_Coord_N,'k*');
% plot(Max_Coord_E,Max_Coord_N,'k*');

E_Rectangle=[Min_Coord_E Min_Coord_E Max_Coord_E Max_Coord_E Min_Coord_E];
N_Rectangle=[Min_Coord_N Max_Coord_N Max_Coord_N Min_Coord_N Min_Coord_N];
% plot(E_Rectangle,N_Rectangle,'k');

% Plot of Explored Zone - drawing a Quadrilateral
% plot(Min_Coord_E,Coord_N_Min_E,'m*');
% plot(Max_Coord_E,Coord_N_Max_E,'m*');
% plot(Coord_E_Min_N,Min_Coord_N,'m*');
% plot(Coord_E_Max_N,Max_Coord_N,'m*');

E_Quad=[Min_Coord_E Coord_E_Max_N Max_Coord_E Coord_E_Min_N Min_Coord_E];
N_Quad=[Coord_N_Min_E Max_Coord_N Coord_N_Max_E Min_Coord_N Coord_N_Min_E];

% plot(E_Quad,N_Quad,'m');

%% ************************************************************************
% ************************** Rectangle ZONE *******************************
% *********************** Distance Computation [m] ************************

% Boundary Rectangle 
% Ref point 1 - Origine (Left Bottom)
B1=Min_Coord_E;
C1=Min_Coord_N;
% Ref point 2 - Left Top
B2=Min_Coord_E;
C2=Max_Coord_N;
% Ref point 3 - Right Top
B3=Max_Coord_E;
C3=Max_Coord_N;
% Ref point 4- Right Bottom
B4=Max_Coord_E;
C4=Min_Coord_N;

% Computation of distances in meter
Lenght_left=acos(sin(deg2rad(B1))*sin(deg2rad(B2))+cos(deg2rad(B1))*cos(deg2rad(B2))*cos(deg2rad(C1-C2)))*6371*1000;
Lenght_top=acos(sin(deg2rad(B2))*sin(deg2rad(B3))+cos(deg2rad(B2))*cos(deg2rad(B3))*cos(deg2rad(C2-C3)))*6371*1000;
Lenght_right=acos(sin(deg2rad(B3))*sin(deg2rad(B4))+cos(deg2rad(B3))*cos(deg2rad(B4))*cos(deg2rad(C3-C4)))*6371*1000;
Lenght_bottom=acos(sin(deg2rad(B4))*sin(deg2rad(B1))+cos(deg2rad(B4))*cos(deg2rad(B1))*cos(deg2rad(C4-C1)))*6371*1000;

%% ************************************************************************
% ************************** Explored ZONE ********************************
% ********************** Distance Computation [m] *************************

% Boundary Rectangle 
% Ref point 1 - Origine (Left Bottom)
Bz1=Min_Coord_E;
Cz1=Coord_N_Min_E;
% Ref point 2 - Left Top
Bz2=Coord_E_Max_N;
Cz2=Max_Coord_N;
% Ref point 3 - Right Top
Bz3=Max_Coord_E;
Cz3=Coord_N_Max_E;
% Ref point 4- Right Bottom
Bz4=Coord_E_Min_N;
Cz4=Min_Coord_N;

% Computation of distances in meter
Lenght_z_left=acos(sin(deg2rad(Bz1))*sin(deg2rad(Bz2))+cos(deg2rad(Bz1))*cos(deg2rad(Bz2))*cos(deg2rad(Cz1-Cz2)))*6371*1000;
Lenght_z_top=acos(sin(deg2rad(Bz2))*sin(deg2rad(Bz3))+cos(deg2rad(Bz2))*cos(deg2rad(Bz3))*cos(deg2rad(Cz2-Cz3)))*6371*1000;
Lenght_z_right=acos(sin(deg2rad(Bz3))*sin(deg2rad(Bz4))+cos(deg2rad(Bz3))*cos(deg2rad(Bz4))*cos(deg2rad(Cz3-Cz4)))*6371*1000;
Lenght_z_bottom=acos(sin(deg2rad(Bz4))*sin(deg2rad(Bz1))+cos(deg2rad(Bz4))*cos(deg2rad(Bz1))*cos(deg2rad(Cz4-Cz1)))*6371*1000;


%% ************************************************************************
% ************************** Rectangle ZONE *******************************
% ****************** Modification of Coordonate --> m *********************

Coord_E_Zone2=Coord_E(N_data1+Sample_Between+40:N_data-451);
Coord_N_Zone2=Coord_N(N_data1+Sample_Between+40:N_data-451);

% Determiantion of Coordonate x (E) and y (N) in [m]
Coord_E_Zone2_m=acos(sin(deg2rad(B1))*sin(deg2rad(Coord_E_Zone2))+cos(deg2rad(B1))*cos(deg2rad(Coord_E_Zone2))*cos(deg2rad(C1-C1)))*6371*1000;
Coord_N_Zone2_m=acos(sin(deg2rad(B1))*sin(deg2rad(B1))+cos(deg2rad(B1))*cos(deg2rad(B1))*cos(deg2rad(C1-Coord_N_Zone2)))*6371*1000;

%% *****************  DETERMINATION ZONE AUTOMATIC ************************

% Determination Boundaries
Index_Bound1 = boundary(Coord_E_Zone2_m,Coord_N_Zone2_m);  % Determination of the Boundary of the Zone
Index_Bound = boundary(Coord_E_Zone2_m,Coord_N_Zone2_m,0.05); % Determination of the Boundary of the Zone less detailled - 0.05
% plot(Coord_E_Zone2_m(Index_Bound1),Coord_N_Zone2_m(Index_Bound1),'m','LineWidth',2);
% plot(Coord_E_Zone2_m(Index_Bound),Coord_N_Zone2_m(Index_Bound),'k','LineWidth',2);

% Generation of the Polyhedron
Poly_Boundaries=[Coord_E_Zone2_m(Index_Bound) Coord_N_Zone2_m(Index_Bound)];
Poly_Creation=Polyhedron(Poly_Boundaries);

%% ************************** MESHING *************************************
% Total coverage Zone - Square
Quadrillage_E=0:round(Lenght_top/Step); 
Quadrillage_N=0:round(Lenght_left/Step);

% Determination of the Rectangle limits
Xmin=Quadrillage_E(1)*Step;
Xmax=Quadrillage_E(end)*Step;
Ymin=Quadrillage_N(1)*Step;
Ymax=Quadrillage_N(end)*Step;

XYmax=max(Xmax,Ymax);
 
% Creation Polyhedron
Poly_Rectangle=Polyhedron('lb',[Xmin, Ymin],'ub',[XYmax, XYmax]);

% Determination of the meshing step
Y_mesh=Ymin:Step*(1+cosd(60))/2:XYmax;
X1_mesh=Xmin:Step*cosd(30):XYmax;
X2_mesh=Step*cosd(30)/2:Step*cosd(30):XYmax;

%% Creation of Meshing centers
PT_Centers=[];
for i=1:length(Y_mesh)
    if mod(i,2)==1
        y_temp=ones(1,length(X1_mesh))*Y_mesh(i);
        PT_Centers=[PT_Centers; X1_mesh' y_temp'];
    else
        y_temp=ones(1,length(X2_mesh))*Y_mesh(i);
        PT_Centers=[PT_Centers; X2_mesh' y_temp'];
    end
end

%% ****************** Checking for Point inside the Zone ******************
% Creation of Vector of Position
X0_H=PT_Centers;
% Building of an Empty Vector
X_H_inside=0*X0_H;

% Detection if the Points are inside the Polyhedron
for i=1:length(X0_H)
    if isnan(X0_H(i,:))==0
         X_H_inside(i,:)=isInside(Poly_Creation,X0_H(i,:)');
    end
end

% Keep the Coordonnates of the points inside the Polyhedron
Position_H=X0_H.*X_H_inside;

% Detection of the index of the Points that are INside the Polyhedron
index_H=[find(Position_H(:,1)~=0) find(Position_H(:,1)~=0)];

% Keep the Points that are Inside the Polyhedron
j=1;
for i=1:length(index_H)
    Position2_H(j,:)=Position_H(index_H(i),:);
    j=j+1;
end

% Keep the Point that are not NaN
j=1;
for i=1:length(Position2_H)
    if isnan(Position2_H(i,1))==1   

    else
        Position3_H(j,:)=Position2_H(i,:);
        j=j+1;
    end
end

Position_H=Position3_H;
clear Position2_H Position3_H

%% ****************** Determination of the Local Zone *********************
% Creation of Hexagonal Zone with relative coordonnates
%Step=Step*1.25 % ***** For OVERLAPPING *****

Center=[0 0];
Pt_x=[-(Step*cosd(30))/2 0      (Step*cosd(30))/2 (Step*cosd(30))/2 0       -(Step*cosd(30))/2]; % x coord
Pt_y=[(Step*cosd(60))/2  Step/2 (Step*cosd(60))/2 -(Step*cosd(60))/2  -Step/2 -(Step*cosd(60))/2]; % y coord
Pt=[Pt_x' Pt_y']; % Determination of the points at the boundary

% Creation of the Hexagone Zone
Poly_Hexagone=Polyhedron(Pt);

% for i=1:length(Position_H)
%     pos_H=[Position_H(i,1),Position_H(i,2)]';
%     Zone(i)=pos_H+Poly_Hexagone;
% end

% figure(1);
% subplot(224)
% plot(Poly_Creation)
% hold on; grid on;
% plot(Zone)
% plot(Position_H(:,1),Position_H(:,2),'k*')
% plot(Coord_E_Zone2_m,Coord_N_Zone2_m,'k','LineWidth',2);
% plot(Coord_E_Zone2_m(Index_Bound),Coord_N_Zone2_m(Index_Bound),'g','LineWidth',2);
% title('Trajectory')
% xlabel('E')
% ylabel('N')

% size(Zone)

%% ****************** Checking for Point inside the Zone ******************

clear Position
clear Label

Position=[];
Label=[];

coordi = 1;
coordj = 1;

last_y = Position_H(1,2);

for ind_positionH = 1:size(Position_H,1)
    
    if inpolygon(Position_H(ind_positionH,1),Position_H(ind_positionH,2),Poly_Boundaries(:,1),Poly_Boundaries(:,2))==1
        
        Position=[Position; [Position_H(ind_positionH,1) Position_H(ind_positionH,2)]];
        Label=[Label; coordi,coordj];
           
    end  
    
    coordj = coordj + 1;

    if Position_H(ind_positionH,2) > last_y

        coordi = coordi + 1;
        coordj = 1;

    end

    last_y = Position_H(ind_positionH,2);
           
end

%% ****************** Route creation *****************

route=[];

Position_sorted=[];

horizontalroute = 1; % 1 Ruta horizontal
                     % 0 Ruta vertical
                    
% Options considering the route selected was horizontal
                    
route_type_a = 1; % route_type = 1 -> right  left  right....
                  % route_type = 2 -> right  left  right....

route_type_b = 1; % route_type = 1 -> down to up
                  % route_type = 2 -> up to down

if horizontalroute == 1

    if route_type_b == 1

        rango_fila = min(Label(:,1)):max(Label(:,1));

    elseif route_type_b == 2

        rango_fila = max(Label(:,1)):-1:min(Label(:,1));    

    end

    if route_type_a == 1

        for fila=rango_fila

            celdas_fila=(Label(:,1)==fila).*Label;
            celdas_fila(~any(celdas_fila,2),:)=[];

            new_Position_sorted=(Label(:,1)==fila).*Position;
            new_Position_sorted(~any(new_Position_sorted,2),:)=[];

            if fila == 1
            new_Position_sorted = [0 0; new_Position_sorted];
            end

            if mod(fila,2)==1

                [~, nuevas_celdas_index] = sort(celdas_fila(:, 2),'ascend');

            elseif mod(fila,2)==0

                [~, nuevas_celdas_index] = sort(celdas_fila(:, 2),'descend');

            end

            nuevas_celdas=celdas_fila(nuevas_celdas_index,:);

            new_Position_sorted=new_Position_sorted(nuevas_celdas_index,:);

            route=[route; nuevas_celdas];
            Position_sorted=[Position_sorted; new_Position_sorted];

        end

    elseif route_type_a == 2

        for fila=rango_fila

            celdas_fila=(Label(:,1)==fila).*Label;
            celdas_fila(~any(celdas_fila,2),:)=[];

            new_Position_sorted=(Label(:,1)==fila).*Position;
            new_Position_sorted(~any(new_Position_sorted,2),:)=[];

            if fila == 1
            new_Position_sorted = [0 0; new_Position_sorted];
            end

            if mod(fila,2)==1

                [~, nuevas_celdas_index] = sort(celdas_fila(:, 2),'descend');

            elseif mod(fila,2)==0

                [~, nuevas_celdas_index] = sort(celdas_fila(:, 2),'ascend');

            end

            nuevas_celdas=celdas_fila(nuevas_celdas_index,:);

            new_Position_sorted=new_Position_sorted(nuevas_celdas_index,:);

            route=[route; nuevas_celdas];
            Position_sorted=[Position_sorted; new_Position_sorted];

        end        

    end
    
end

% Vertical route

if horizontalroute == 0

    if route_type_b == 1

        rango_column = min(Label(:,2)):max(Label(:,2));

    elseif route_type_b == 2

        rango_column = max(Label(:,2)):-1:min(Label(:,2));    

    end
          
    if route_type_a == 1
    
        for column=rango_column

            celdas_column=(Label(:,2)==column).*Label;
            celdas_column(~any(celdas_column,2),:)=[];

            new_Position_sorted=(Label(:,2)==column).*Position;
            new_Position_sorted(~any(new_Position_sorted,2),:)=[];

            if column == 1
            new_Position_sorted = [0 0; new_Position_sorted];
            end

            if mod(column,2)==1

                [~, nuevas_celdas_index] = sort(celdas_column(:, 1),'ascend');

            elseif mod(column,2)==0

                [~, nuevas_celdas_index] = sort(celdas_column(:, 1),'descend');

            end

            nuevas_celdas=celdas_column(nuevas_celdas_index,:);

            new_Position_sorted=new_Position_sorted(nuevas_celdas_index,:);

            route=[route; nuevas_celdas];
            Position_sorted=[Position_sorted; new_Position_sorted];

        end    
    
    elseif route_type_a == 2
        
        for column=rango_column

            celdas_column=(Label(:,2)==column).*Label;
            celdas_column(~any(celdas_column,2),:)=[];

            new_Position_sorted=(Label(:,2)==column).*Position;
            new_Position_sorted(~any(new_Position_sorted,2),:)=[];

            if column == 1
            new_Position_sorted = [0 0; new_Position_sorted];
            end

            if mod(column,2)==1

                [~, nuevas_celdas_index] = sort(celdas_column(:, 1),'descend');

            elseif mod(column,2)==0

                [~, nuevas_celdas_index] = sort(celdas_column(:, 1),'ascend');

            end

            nuevas_celdas=celdas_column(nuevas_celdas_index,:);

            new_Position_sorted=new_Position_sorted(nuevas_celdas_index,:);

            route=[route; nuevas_celdas];
            Position_sorted=[Position_sorted; new_Position_sorted];

        end    
        
    end
    
end


%% Zone Creation

%clear Zone
P = Polyhedron('lb',[-0.5 -0.5],'ub',[0.5 0.5]);
for i=1:length(Position_sorted)
    pos_H=[Position_sorted(i,1),Position_sorted(i,2)]';
    Zone_parcours(i)=pos_H+P;
end

