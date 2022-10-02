#include "ros/ros.h"
#include "tf/tf.h"
#include <iostream>
#include <cmath>
#include <time.h> 
#include <ctime>
#include "geometry_msgs/Point.h"
#include "geometry_msgs/Twist.h"
#include "visualization_msgs/Marker.h"
#include "std_msgs/Float64.h"
#include "math.h"
using namespace std;

//Paramètres
std::vector<double> X {33,29,0,0,0,0};
double input_g ;
double input_d ;
clock_t start = clock();
clock_t end ;
double dt = 1; //pour initialiser le départ


void calcul_stateCallback(const geometry_msgs::Point::ConstPtr& msg1)
{	
	//mis a jour des inputs 
	input_g = msg1 ->x;
	input_d = msg1 ->y;	
	
	//chrono entre deux reception d'input
	clock_t end = clock();
	
	//euler integration modele dynamique a partir de la thèse 2016 
	//"Design and Implementation of an Effective Communication and Coordination System forUnmanned Surface Vehicles (USV)"
	//de Yoann Hervagault
	
	X[5] = X[5] + dt*( ((0.127/2)*(input_g-input_d))/0.8 - (0.28/0.33)*X[5] );
	X[3] = X[3] + dt*( (input_g + input_d)/28 + (28/28)*X[4]*X[5] -(8.5/28)*X[3] );
	X[4] = X[4] + dt*( (-28/28)*X[3]*X[5] - (19/28)*X[4] );
	X[2] = X[2] + dt*( X[5]);
	X[0] = X[0] + dt*( X[3]*cos(X[2]) - X[4]*sin(X[2]) );
	X[1] = X[1] + dt*( X[3]*sin(X[2]) - X[4]*cos(X[2]) );
	dt =0.5;//  (end - start);
	clock_t start = end;
	

}



int main(int argc, char **argv){

		// Initialisation du node : le troisième argument est son nom
		ros::init(argc, argv, "data_sensors1");

		ros::NodeHandle n1;

		ros::Publisher chatter_sensors = n1.advertise<geometry_msgs::Twist>("/data_sensors", 1000);

		ros::Publisher chatter_time = n1.advertise<std_msgs::Float64>("/sample_time", 1000);
		
		ros::Publisher chatter_x = n1.advertise<std_msgs::Float64>("/state_x_bds", 1000);
		
		ros::Publisher chatter_y = n1.advertise<std_msgs::Float64>("/state_y_bds", 1000);
		
		ros::Subscriber sub_control = n1.subscribe("/input_vessel", 1000, calcul_stateCallback);

		ros::Publisher vis_pub = n1.advertise<visualization_msgs::Marker>( "visualization_marker", 0 );

		//ros::spin();
		ros::Rate loop_rate(1); //frequence lente pour les calculs de casadi

		while (ros::ok()){

			
			std_msgs::Float64 x_state;
			std_msgs::Float64 y_state;
			geometry_msgs::Twist msg;
			x_state.data = X[0];
			y_state.data = X[1];
			msg.linear.x = X[0];
			msg.linear.y = X[1];
			msg.linear.z = X[2];
			msg.angular.x = X[3];
			msg.angular.y = X[4];
			msg.angular.z = X[5];
			

			visualization_msgs::Marker marker;

			tf::Quaternion q;
			q.setRPY(0, 0, X[2]);
			
			marker.header.frame_id = "map";
			marker.header.stamp = ros::Time::now();
			std::string ns = ros::this_node::getNamespace();
			marker.ns = ns;
			marker.id = 0;
			marker.type = visualization_msgs::Marker::MESH_RESOURCE;
			marker.action = visualization_msgs::Marker::ADD;
			marker.pose.position.x = X[0];
 			marker.pose.position.y = X[1];
 			marker.pose.position.z = 0;
 			tf::quaternionTFToMsg(q, marker.pose.orientation);
			marker.scale.x = 0.5;
			marker.scale.y = 0.5;
			marker.scale.z = 0.5;
			marker.color.a = 1.0; // alpha = transparence
			marker.color.r = 1.0;
			marker.color.g = 1.0;
			marker.color.b = 1.0;
			marker.mesh_resource = "package://matlab_ros/mesh/boat.dae";

			std_msgs::Float64 time;
			time.data = dt;

			// publication du message
			chatter_sensors.publish(msg);
			chatter_time.publish(time);
			vis_pub.publish( marker );
			chatter_x.publish(x_state);
			chatter_y.publish(y_state);
			
			ros::spinOnce();
			//pause
			loop_rate.sleep();

		}
	}
