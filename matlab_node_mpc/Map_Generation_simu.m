% Map_Generation(Step,Shape)
% Step is the size of the meshing element - ex : 2m
% Shape = 0 is Square element; 
% Shape = 1 is Hexagonal element; 


% NEW CODE ##### ZONE 2 #####
function Map_Generation_simu(Step,Shape)

%Step=5
%Shape=1
load Lac_Heron_2021.mat
load donox.mat
load dont.mat
load donx.mat
load dony.mat
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

%% **************** Sample per experiment **********************
% N_start=1994;%00;  % Start Sample - 1 min
% N_data1=4846; % End of the 1st experiment
% N_data2=4955; % End of the 2nd experiment
% Sample_Between=N_data2-N_data1+160;
% N_data=9031;  % End Sample 

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
Coord_E_Zone2_m=dony';
Coord_N_Zone2_m=donx';

%% *****************  DETERMINATION ZONE AUTOMATIC ************************
% figure(112)
%  hold on
%  grid on
% plot(Coord_E_Zone2_m,Coord_N_Zone2_m,'k');
% %plot(Coord_E(N_data1+Sample_Between+40:N_data-451),Coord_N(N_data1+Sample_Between+40:N_data-451),'r');
% %plot(Coord_E(N_data1+Sample_Between),Coord_N(N_data1+Sample_Between),'b*');
% title('Trajectory')
% xlabel('E')
% ylabel('N')

% Determination Boundaries
Index_Bound1 = boundary(Coord_E_Zone2_m,Coord_N_Zone2_m);  % Determination of the Boundary of the Zone
Index_Bound = boundary(Coord_E_Zone2_m,Coord_N_Zone2_m,0.05); % Determination of the Boundary of the Zone less detailled - 0.05
% plot(Coord_E_Zone2_m(Index_Bound1),Coord_N_Zone2_m(Index_Bound1),'m','LineWidth',2);
% plot(Coord_E_Zone2_m(Index_Bound),Coord_N_Zone2_m(Index_Bound),'k','LineWidth',2);

% Generation of the Polyhedron
Poly_Boundaries=[Coord_E_Zone2_m(Index_Bound) Coord_N_Zone2_m(Index_Bound)];
Poly_Creation=Polyhedron(Poly_Boundaries);

figure(113);
plot(Poly_Creation);
 hold on
 grid on
plot(Coord_E_Zone2_m,Coord_N_Zone2_m,'b');
plot(Coord_E_Zone2_m(Index_Bound1),Coord_N_Zone2_m(Index_Bound1),'g','LineWidth',2);
plot(Coord_E_Zone2_m(Index_Bound),Coord_N_Zone2_m(Index_Bound),'c','LineWidth',2);
title('Trajectory')
xlabel('E')
ylabel('N')


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

switch Shape
    case 0  % Rectangle Meshing
        %% Creation of Meshing centers
        % Determination of the meshing step
        Mesh_Step=XYmax/Step+1;
        % Meshing
        [X, Y]= meshGrid(Poly_Rectangle,Mesh_Step);

        %% ****************** Checking for Point inside the Zone ******************
        % Creation of Vector of Position
        X0=[X(:) Y(:)];
        % Building of an Empty Vector
        X_inside=0*X0;

        % Detection if the Points are inside the Polyhedron
        for i=1:length(X0)
            if isnan(X0(i,:))==0
                 X_inside(i,:)=isInside(Poly_Creation,X0(i,:)');
            end
        end

        % Keep the Coordonnates of the points inside the Polyhedron
        Position=X0.*X_inside;

        % Detection of the index of the Points that are INside the Polyhedron
        index=[find(Position(:,1)~=0) find(Position(:,1)~=0)];

        % Keep the Points that are Inside the Polyhedron
        j=1;
        for i=1:length(index)
            Position2(j,:)=Position(index(i),:);
            j=j+1;
        end

        % Keep the Point that are not NaN
        j=1;
        for i=1:length(Position2)
            if isnan(Position2(i,1))==1   

            else
                Position3(j,:)=Position2(i,:);
                j=j+1;
            end
        end

        Position=Position3;
        clear Position2 Position3

        %% ****************** Determination of the Local Zone *********************

        Mini_Poly=Polyhedron('lb',[-Step/2, -Step/2],'ub',[Step/2, Step/2]);

        for i=1:length(Position)
            pos=[Position(i,1),Position(i,2)]';
            Zone(i)=pos+Mini_Poly;
        end

        figure(114);
        plot(Poly_Creation)
        hold on; grid on;
%         plot(Zone)
%         plot(Position(:,1),Position(:,2),'k*')
%         plot(Coord_E_Zone2_m,Coord_N_Zone2_m,'k','LineWidth',2);
%         plot(Coord_E_Zone2_m(Index_Bound),Coord_N_Zone2_m(Index_Bound),'g','LineWidth',2);
        title('Trajectory')
        xlabel('E')
        ylabel('N')

%         figure(1114);
%         subplot(221)
%         plot(Poly_Creation)
%         hold on; grid on;
%         plot(Zone)
%          plot(Position(:,1),Position(:,2),'k*')
%         plot(Coord_E_Zone2_m,Coord_N_Zone2_m,'k','LineWidth',2);
%         plot(Coord_E_Zone2_m(Index_Bound),Coord_N_Zone2_m(Index_Bound),'g','LineWidth',2);
%         title('Trajectory')
%         xlabel('E')
%         ylabel('N')

        



    case 1 % Hexagonal Meshing
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

        for i=1:length(Position_H)
            pos_H=[Position_H(i,1),Position_H(i,2)]';
            %fprintf("(%f),\n",pos_H);
            Zone(i)=pos_H+Poly_Hexagone;
        end

        figure(1114);
        %subplot(224)
        plot(Poly_Creation)
        hold on; grid on;
        plot(Zone)
        plot(Position_H(:,1),Position_H(:,2),'k*')
        plot(Coord_E_Zone2_m,Coord_N_Zone2_m,'k','LineWidth',2);
        plot(Coord_E_Zone2_m(Index_Bound),Coord_N_Zone2_m(Index_Bound),'g','LineWidth',2);
        title('Trajectory')
        xlabel('E')
        ylabel('N')

        size(Zone)


    otherwise % None
        a=pi
end



%% ****************** Association Measures to each Zone *******************
Time_Sample=[Time(N_data1+Sample_Between+40:N_data-451)]; % Instant of the sample
Temp1_Dyn=[Temp1(N_data1+Sample_Between+40:N_data-451)]; % Temperature for the Zone
Oxy_Dyn=donox'; % Oxygene for the Zone
pH_Dyn=[pH(N_data1+Sample_Between+40:N_data-451)]; % pH for the Zone 
Redox_Dyn=[Redox(N_data1+Sample_Between+40:N_data-451)]; % Redox for the Zone
Conduc_Dyn=[Conduc(N_data1+Sample_Between+40:N_data-451)]; % Conduc for the Zone
Turbi_Filter_B_Dyn=[Turbi_Filter_B(N_data1+Sample_Between+40:N_data-451)]; % Turbi_Filter_B for the Zone

x=Coord_E_Zone2_m';
y=Coord_N_Zone2_m';

% Coordonnates of Points belong the Explored Zone
coord_xy=[x(:) y(:)];

% Building of a Empty Vectors and detection of the Points inside each Zone
for k=1:length(Zone) 
    X_Zone(k).X_Zone_inside=0*coord_xy;
    for i=1:length(coord_xy)
         X_Zone(k).X_Zone_inside(i,:)=isInside(Zone(k),coord_xy(i,:)');
    end
    % Determination of mean of Parameter inside each Zone
   X_Zone(k).Time=Time_Sample(find(X_Zone(k).X_Zone_inside(:,1)==1));  % Instants of the sample per zone
   X_Zone(k).Temperature=mean(Temp1_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
    X_Zone(k).Oxygene=mean(Oxy_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
%    X_Zone(k).pH=mean(pH_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
%    X_Zone(k).Redox=mean(Redox_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
%    X_Zone(k).Conduc=mean(Conduc_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
%    X_Zone(k).Turbidity=mean(Turbi_Filter_B_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
   %fprintf("(%f)",mean(Oxy_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1))));
   X_Zone(k).Temperature_Var=var(Temp1_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
   X_Zone(k).Oxygene_Var=var(Oxy_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
%    X_Zone(k).pH_Var=var(pH_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
%    X_Zone(k).Redox_Var=var(Redox_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
%    X_Zone(k).Conduc_Var=var(Conduc_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
%    X_Zone(k).Turbidity_Var=var(Turbi_Filter_B_Dyn(find(X_Zone(k).X_Zone_inside(:,1)==1)));
   
end


%% ************************************************************************
% **************************  PLOT GRADIENT   *****************************
%% ************************************************************************

% %% PLOT OF ZONES WITH GRADIENT OF TEMPERATURE
% for i=1:length(Zone) 
% Temperature_Zone(i)=X_Zone(i).Temperature;
% end
% 
% min_Temp=min(Temperature_Zone); % Min of the Temperature
% max_Temp=max(Temperature_Zone); % Max of the Temperature
% 
% Index_min=0; % Min Color
% Index_max=255;
% Index=round((Index_max/(max_Temp-min_Temp))*Temperature_Zone-(Index_max/(max_Temp-min_Temp))*min_Temp,0); % Computation of the Color depending of the Temperature at time k
% 
% index_color=(max(Index)-min(Index))/(5*255); % Number of intervalles
% Ticks_color=min(Index)/255:index_color:max(Index)/255; % Number of Ticks
% %index_color=(max_Temp-min_Temp)/5; % Number of intervalles
% %Ticks_color=min_Temp:index_color:max_Temp
% Ticks_Value=round(Ticks_color*255*((max_Temp-min_Temp)/max(Index))+min_Temp,2); % Values of Ticks
% 
% Color_Range=colormap;
% 
% for i=1:length(Index) 
%     if isnan(Index(i))==1  % White Zone for Nan
%         lineColored(i,1)=1;
%         lineColored(i,2)=1;
%         lineColored(i,3)=1;
%     else
%         lineColored(i,:)=Color_Range(Index(i)+1,:);
%     end
% end
% 
%  figure(120);
%  for i=1:length(Zone) 
%   hold on
%   grid on
% plot(Zone(i),'color', lineColored(i,1:3),...
%     'FaceColor', 'no',...
% 	'EdgeColor', 'm',...
% 	'LineWidth', 1);
%  end
%  title('Temperature - Zone 2')
%  xlabel('E')
%  ylabel('N')
%     c_bar=colorbar('Ticks',Ticks_color,...
%          'TickLabels',Ticks_Value);
% c_bar.Label.String = 'Temperature';

 %% PLOT OF ZONES WITH GRADIENT OF OXYGENE
for i=1:length(Zone) 
Oxygene_Zone(i)=X_Zone(i).Oxygene;
end

min_Temp=min(Oxygene_Zone); % Min of the Temperature
max_Temp=max(Oxygene_Zone); % Max of the Temperature

Index_min=0; % Min Color
Index_max=255;
Index=round((Index_max/(max_Temp-min_Temp))*Oxygene_Zone-(Index_max/(max_Temp-min_Temp))*min_Temp,0); % Computation of the Color depending of the Temperature at time k

index_color=(max(Index)-min(Index))/5/255; % Number of intervalles
Ticks_color=min(Index)/255:index_color:max(Index)/255; % Number of Ticks
Ticks_Value=round(Ticks_color*255*((max_Temp-min_Temp)/max(Index))+min_Temp,2); % Values of Ticks


Color_Range=colormap;

for i=1:length(Index) 
    if isnan(Index(i))==1  % White Zone for Nan
        lineColored(i,1)=1;
        lineColored(i,2)=1;
        lineColored(i,3)=1;
    else
        lineColored(i,:)=Color_Range(Index(i)+1,:);
    end
end

 figure(121);
 for i=1:length(Zone) 
  hold on
  grid on
plot(Zone(i),'color', lineColored(i,1:3),...
    'FaceColor', 'no',...
	'EdgeColor', 'm',...
	'LineWidth', 1);
 end
 title('Oxygene - Zone 2')
 xlabel('E')
 ylabel('N')
    c_bar=colorbar('Ticks',Ticks_color,...
         'TickLabels',Ticks_Value);
c_bar.Label.String = 'Oxygene';

% %% PLOT OF ZONES WITH GRADIENT OF pH
% for i=1:length(Zone) 
% pH_Zone(i)=X_Zone(i).pH;
% end
% 
% min_Temp=min(pH_Zone); % Min of the Temperature
% max_Temp=max(pH_Zone); % Max of the Temperature
% 
% Index_min=0; % Min Color
% Index_max=255;
% Index=round((Index_max/(max_Temp-min_Temp))*pH_Zone-(Index_max/(max_Temp-min_Temp))*min_Temp,0); % Computation of the Color depending of the Temperature at time k
% 
% index_color=(max(Index)-min(Index))/5/255; % Number of intervalles
% Ticks_color=min(Index)/255:index_color:max(Index)/255; % Number of Ticks
% Ticks_Value=round(Ticks_color*255*((max_Temp-min_Temp)/max(Index))+min_Temp,2); % Values of Ticks
% 
% Color_Range=colormap;
% 
% for i=1:length(Index) 
%     if isnan(Index(i))==1  % White Zone for Nan
%         lineColored(i,1)=1;
%         lineColored(i,2)=1;
%         lineColored(i,3)=1;
%     else
%         lineColored(i,:)=Color_Range(Index(i)+1,:);
%     end
% end
% 
%  figure(122);
%  for i=1:length(Zone) 
%   hold on
%   grid on
% plot(Zone(i),'color', lineColored(i,1:3),...
%     'FaceColor', 'no',...
% 	'EdgeColor', 'm',...
% 	'LineWidth', 1);
%  end
%  title('pH - Zone 2')
%  xlabel('E')
%  ylabel('N')
%    c_bar=colorbar('Ticks',Ticks_color,...
%          'TickLabels',Ticks_Value);
% c_bar.Label.String = 'pH';
% 
% 
%  %% PLOT OF ZONES WITH GRADIENT OF REDOX
% for i=1:length(Zone) 
% Redox_Zone(i)=X_Zone(i).Redox;
% end
% 
% min_Temp=min(Redox_Zone); % Min of the Temperature
% max_Temp=max(Redox_Zone); % Max of the Temperature
% 
% Index_min=0; % Min Color
% Index_max=255;
% Index=round((Index_max/(max_Temp-min_Temp))*Redox_Zone-(Index_max/(max_Temp-min_Temp))*min_Temp,0); % Computation of the Color depending of the Temperature at time k
% 
% index_color=(max(Index)-min(Index))/5/255; % Number of intervalles
% Ticks_color=min(Index)/255:index_color:max(Index)/255; % Number of Ticks
% Ticks_Value=round(Ticks_color*255*((max_Temp-min_Temp)/max(Index))+min_Temp,2); % Values of Ticks
% 
% Color_Range=colormap;
% 
% for i=1:length(Index) 
%     if isnan(Index(i))==1  % White Zone for Nan
%         lineColored(i,1)=1;
%         lineColored(i,2)=1;
%         lineColored(i,3)=1;
%     else
%         lineColored(i,:)=Color_Range(Index(i)+1,:);
%     end
% end
% 
%  figure(123);
%  for i=1:length(Zone) 
%   hold on
%   grid on
% plot(Zone(i),'color', lineColored(i,1:3),...
%     'FaceColor', 'no',...
% 	'EdgeColor', 'm',...
% 	'LineWidth', 1);
%  end
%  title('Redox - Zone 2')
%  xlabel('E')
%  ylabel('N')
%   c_bar=colorbar('Ticks',Ticks_color,...
%          'TickLabels',Ticks_Value);
% c_bar.Label.String = 'Redox';
% 
% %% PLOT OF ZONES WITH GRADIENT OF CONDUC
% for i=1:length(Zone) 
% Conduc_Zone(i)=X_Zone(i).Conduc;
% end
% 
% min_Temp=min(Conduc_Zone); % Min of the Temperature
% max_Temp=max(Conduc_Zone); % Max of the Temperature
% 
% Index_min=0; % Min Color
% Index_max=255;
% Index=round((Index_max/(max_Temp-min_Temp))*Conduc_Zone-(Index_max/(max_Temp-min_Temp))*min_Temp,0); % Computation of the Color depending of the Temperature at time k
% 
% index_color=(max(Index)-min(Index))/5/255; % Number of intervalles
% Ticks_color=min(Index)/255:index_color:max(Index)/255; % Number of Ticks
% Ticks_Value=round(Ticks_color*255*((max_Temp-min_Temp)/max(Index))+min_Temp,2); % Values of Ticks
% 
% Color_Range=colormap;
% 
% for i=1:length(Index) 
%     if isnan(Index(i))==1  % White Zone for Nan
%         lineColored(i,1)=1;
%         lineColored(i,2)=1;
%         lineColored(i,3)=1;
%     else
%         lineColored(i,:)=Color_Range(Index(i)+1,:);
%     end
% end
% 
%  figure(124);
%  for i=1:length(Zone) 
%   hold on
%   grid on
% plot(Zone(i),'color', lineColored(i,1:3),...
%     'FaceColor', 'no',...
% 	'EdgeColor', 'm',...
% 	'LineWidth', 1);
%  end
%  title('Conductivity - Zone 2')
%  xlabel('E')
%  ylabel('N')
%  c_bar=colorbar('Ticks',Ticks_color,...
%          'TickLabels',Ticks_Value);
% c_bar.Label.String = 'Conductivity';
% 
%  %% PLOT OF ZONES WITH GRADIENT OF TEMPERATURE
% for i=1:length(Zone) 
% Turbidity_Zone(i)=X_Zone(i).Turbidity;
% end
% 
% min_Temp=min(Turbidity_Zone); % Min of the Temperature
% max_Temp=max(Turbidity_Zone); % Max of the Temperature
% 
% Index_min=0; % Min Color
% Index_max=255;
% Index=round((Index_max/(max_Temp-min_Temp))*Turbidity_Zone-(Index_max/(max_Temp-min_Temp))*min_Temp,0); % Computation of the Color depending of the Temperature at time k
% 
% index_color=(max(Index)-min(Index))/5/255; % Number of intervalles
% Ticks_color=min(Index)/255:index_color:max(Index)/255; % Number of Ticks
% Ticks_Value=round(Ticks_color*255*((max_Temp-min_Temp)/max(Index))+min_Temp,2); % Values of Ticks
% 
% Color_Range=colormap;
% 
% for i=1:length(Index) 
%     if isnan(Index(i))==1  % White Zone for Nan
%         lineColored(i,1)=1;
%         lineColored(i,2)=1;
%         lineColored(i,3)=1;
%     else
%         lineColored(i,:)=Color_Range(Index(i)+1,:);
%     end
% end
% 
%  figure(125);
%  for i=1:length(Zone) 
%   hold on;
%   grid on;
% plot(Zone(i),'color', lineColored(i,1:3),...
%     'FaceColor', 'no',...
% 	'EdgeColor', 'm',...
% 	'LineWidth', 1);
%  end
%  title('Turbidity - Zone 2')
%  xlabel('E')
%  ylabel('N')
% c_bar=colorbar('Ticks',Ticks_color,...
%          'TickLabels',Ticks_Value);
% c_bar.Label.String = 'Turbidity';


%% ************************************************************************
%  ************************* Plot Variance Maps ***************************
compute_var=1;

if compute_var==1;
    run variance_plot(Zone,X_Zone)
end

%% ************************************************************************
%  ************************** Rebuild Data *********************************
%% ************************************************************************

switch Shape
    case 0  % Rectangle Meshing
        run neighbourhood
    case 1 % Hexagonal Meshing
        run neighbourhood_Hexa(Zone,X_Zone)

    otherwise % None
        a=pi*2
end

size(Zone)

 save Mapping_Lac_Heron_2021.mat