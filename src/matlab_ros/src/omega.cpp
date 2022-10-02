#include "ros/ros.h"
#include "tf/tf.h"
#include <iostream>
#include <cmath>
#include "geometry_msgs/PoseStamped.h"
#include "geometry_msgs/Twist.h"
#include "math.h"
#include "visualization_msgs/Marker.h"
using namespace std;


//Coordonnées du target set oméga 
std::vector<double> C {-15, 0}; 



int main(int argc, char **argv){
	    // Initialisation du node 
		ros::init(argc, argv, "omega");

		// Connexion au master et initialisation du NodeHandle qui permet d'avoir accès aux topics et services
		ros::NodeHandle n2;

		ros::Publisher chatter_pub1 = n2.advertise<geometry_msgs::PoseStamped>("omega", 1000);

		ros::Rate loop_rate(25);

		while (ros::ok()){
			
			ros::spinOnce();

			geometry_msgs::PoseStamped omega;

			omega.pose.position.x = C[0];
			omega.pose.position.y = C[1];
			omega.header.frame_id = "map";
			omega.header.stamp = ros::Time::now();

			chatter_pub1.publish(omega);

			loop_rate.sleep();

			}
	}
