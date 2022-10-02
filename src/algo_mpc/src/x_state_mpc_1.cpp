// Node pour un deuxième robot !!
#include "ros/ros.h"

#include <sensor_msgs/NavSatFix.h>
#include <sensor_msgs/Imu.h>
#include <geometry_msgs/Vector3Stamped.h>
#include <geometry_msgs/Twist.h>
#include <std_msgs/Float64.h>
#include <cmath>
#include "tf/tf.h"
#include <iostream>
using namespace std;


std::vector<double> X {0,0,0,0,0,0}; 			//état initial
double lata = -33.728122; //-33.7232726;		//latitude et longitude de l'origine du plan 2D
double lona = 150.670025; //150.6691994;
double re = 6378137;					//rayon de la terre
double dt = 1;
int cnt = 0;
double tot_x;
double tot_y;
double mean_x;
double mean_y;
double vel_x;
double vel_y;


void chatterCallback_bis(const sensor_msgs::NavSatFix& msg_bis)
{
    double lat = msg_bis.latitude;
    double lon = msg_bis.longitude;
    X[0] = (acos(sin(lata)*sin(lata) + cos(lata)*cos(lata)*cos(lon-lona)) * re) / 100;
    X[1] = (acos(sin(lata)*sin(lat) + cos(lata)*cos(lat)) * re) / 100;
}

void chatterCallback1_bis(const sensor_msgs::Imu& msg1_bis)
{
    double yaw = tf::getYaw(msg1_bis.orientation);
    double angular = msg1_bis.angular_velocity.z;
    double a_x = msg1_bis.linear_acceleration.x;
    double a_y = msg1_bis.linear_acceleration.y;
    cnt = cnt + 1;
    X[2] =  - yaw  ;
    X[5] = -angular;
    tot_x = tot_x + a_x;
    tot_y = tot_y + a_y;
    
    mean_x = tot_x/cnt;
    mean_y = tot_y/cnt;
    X[3] = (X[3] + dt *(a_x - 0.14)) / 30;						//euler pour la vitesse linéaire dans le repère du robot, attention : méthode très imprécise !!! 
    X[4] = (X[4] + dt *(a_y + 0.00059)) / 30;
    
    if ( sqrt(vel_x*vel_x + vel_y*vel_y) <1 or sqrt(vel_x*vel_x + vel_y*vel_y) >20 ) //remise à zéro pour éviter une accumulation d'erreur trop importante
    {
    	X[3] = 0;
    	X[4] = 0;
    }
}

void chatterCallback2_bis(const geometry_msgs::Vector3Stamped& msg2_bis)
{
    vel_x = msg2_bis.vector.x;
    vel_y = msg2_bis.vector.y;
    
}

int main(int argc, char **argv){

    	// talker
	ros::init(argc, argv, "x_state_mpc_1"); // initialise le node
	ros::NodeHandle n1; // connection master
    	ros::Subscriber sub_pos_1 = n1.subscribe("/wamv1/sensors/gps/gps/fix", 100, chatterCallback_bis);
    	ros::Subscriber sub_yaw_angular_1 = n1.subscribe("/wamv1/sensors/imu/imu/data", 100, chatterCallback1_bis);
    	ros::Subscriber sub_vel_1 = n1.subscribe("/wamv1/sensors/gps/gps/fix_velocity", 100, chatterCallback2_bis);
    	ros::Publisher chatter_state_1 = n1.advertise<geometry_msgs::Twist>("/vessel_state_1", 1000);
    	ros::Publisher chatter_x_1 = n1.advertise<std_msgs::Float64>("/x_state_vrx_1", 1000);
    	ros::Publisher chatter_y_1 = n1.advertise<std_msgs::Float64>("/y_state_vrx_1", 1000);
    	ros::Publisher chatter_yaw_1 = n1.advertise<std_msgs::Float64>("/yaw_state_vrx_1", 1000);
    	ros::Publisher chatter_x_vel_1 = n1.advertise<std_msgs::Float64>("/x_vel_state_vrx_1", 1000);
    	ros::Publisher chatter_y_vel_1 = n1.advertise<std_msgs::Float64>("/y_vel_state_vrx_1", 1000);
    	ros::Publisher chatter_ang_vel_1 = n1.advertise<std_msgs::Float64>("/ang_vel_state_vrx_1", 1000);
	
	
    ros::Rate loop_rate(0.7);

    while (ros::ok()){
		ros::spinOnce(); // regarde si nouveau message
		
		std_msgs::Float64 x_state_1;
		std_msgs::Float64 y_state_1;
		std_msgs::Float64 yaw_state_1;
		std_msgs::Float64 x_vel_state_1;
		std_msgs::Float64 y_vel_state_1;
		std_msgs::Float64 ang_vel_state_1;
		geometry_msgs::Twist state_1;
		
		state_1.linear.x = X[0];
		state_1.linear.y = X[1];
		state_1.linear.z = X[2];
		state_1.angular.x = X[3];
		state_1.angular.y = X[4];
		state_1.angular.z = X[5];
		chatter_state_1.publish(state_1);
		
		x_state_1.data = X[0];
		y_state_1.data = X[1];
		yaw_state_1.data = X[2];
		x_vel_state_1.data = X[3];
		y_vel_state_1.data = X[4];
		ang_vel_state_1.data = X[5];
		chatter_x_1.publish(x_state_1);
		chatter_y_1.publish(y_state_1);
		chatter_yaw_1.publish(yaw_state_1);
		chatter_x_vel_1.publish(x_vel_state_1);
		chatter_y_vel_1.publish(y_vel_state_1);
		chatter_ang_vel_1.publish(ang_vel_state_1);
		loop_rate.sleep();
	}
}

