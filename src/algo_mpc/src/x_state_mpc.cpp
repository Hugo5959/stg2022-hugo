
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


std::vector<double> X {0,0,0,0,0,0}; //état initial
double lata = -33.728122; //-33.7232726;		//latitude et longitude de l'origine du plan 2D
double lona = 150.670025; //150.6691994;
double re = 6378137;			//rayon de la terre
double dt = 1;
int cnt = 0;
double tot_x;
double tot_y;
double mean_x;
double mean_y;
double vel_x;
double vel_y;


void chatterCallback(const sensor_msgs::NavSatFix& msg)
{
    double lat = msg.latitude;
    double lon = msg.longitude;
    X[0] = (acos(sin(lata)*sin(lata) + cos(lata)*cos(lata)*cos(lon-lona)) * re) / 100;
    X[1] = (acos(sin(lata)*sin(lat) + cos(lata)*cos(lat)) * re) / 100;
}

void chatterCallback1(const sensor_msgs::Imu& msg1)
{
    double yaw = tf::getYaw(msg1.orientation);
    double angular = msg1.angular_velocity.z;
    double a_x = msg1.linear_acceleration.x;
    double a_y = msg1.linear_acceleration.y;
    cnt = cnt + 1;
    X[2] =  - yaw  ;
    X[5] = -angular;
    tot_x = tot_x + a_x;
    tot_y = tot_y + a_y;
    
    mean_x = tot_x/cnt;
    mean_y = tot_y/cnt;
    X[3] = (X[3] + dt *(a_x - 0.14)) / 30;		//euler pour la vitesse linéaire dans le repère du robot très imprécis !!! 
    X[4] = (X[4] + dt *(a_y + 0.00059)) / 30;
    
    if ( sqrt(vel_x*vel_x + vel_y*vel_y) <1 or sqrt(vel_x*vel_x + vel_y*vel_y) >20 ) //remise à zéro pour éviter une accumulation d'erreur trop importante
    {
    	X[3] = 0;
    	X[4] = 0;
    }
}

void chatterCallback2(const geometry_msgs::Vector3Stamped& msg2)
{
    vel_x = msg2.vector.x;
    vel_y = msg2.vector.y;
    
}

int main(int argc, char **argv){

    // talker
	ros::init(argc, argv, "x_state_mpc"); // innitialise la node
	ros::NodeHandle n; // connection master
    	ros::Subscriber sub_pos = n.subscribe("/wamv/sensors/gps/gps/fix", 100, chatterCallback);
    	ros::Subscriber sub_yaw_angular = n.subscribe("/wamv/sensors/imu/imu/data", 100, chatterCallback1);
    	ros::Subscriber sub_vel = n.subscribe("/wamv/sensors/gps/gps/fix_velocity", 100, chatterCallback2);
    	ros::Publisher chatter_state = n.advertise<geometry_msgs::Twist>("/vessel_state", 1000);
    	ros::Publisher chatter_x = n.advertise<std_msgs::Float64>("/x_state_vrx", 1000);
    	ros::Publisher chatter_y = n.advertise<std_msgs::Float64>("/y_state_vrx", 1000);
    	ros::Publisher chatter_yaw = n.advertise<std_msgs::Float64>("/yaw_state_vrx", 1000);
    	ros::Publisher chatter_x_vel = n.advertise<std_msgs::Float64>("/x_vel_state_vrx", 1000);
    	ros::Publisher chatter_y_vel = n.advertise<std_msgs::Float64>("/y_vel_state_vrx", 1000);
    	ros::Publisher chatter_ang_vel = n.advertise<std_msgs::Float64>("/ang_vel_state_vrx", 1000);
	
	
    ros::Rate loop_rate(0.7);

    while (ros::ok()){
		ros::spinOnce(); // regarde si nouveau message
		
		std_msgs::Float64 x_state;
		std_msgs::Float64 y_state;
		std_msgs::Float64 yaw_state;
		std_msgs::Float64 x_vel_state;
		std_msgs::Float64 y_vel_state;
		std_msgs::Float64 ang_vel_state;
		geometry_msgs::Twist state;
		
		state.linear.x = X[0];
		state.linear.y = X[1];
		state.linear.z = X[2];
		state.angular.x = X[3];
		state.angular.y = X[4];
		state.angular.z = X[5];
		chatter_state.publish(state);
		
		x_state.data = X[0];
		y_state.data = X[1];
		yaw_state.data = X[2];
		x_vel_state.data = X[3];
		y_vel_state.data = X[4];
		ang_vel_state.data = X[5];
		chatter_x.publish(x_state);
		chatter_y.publish(y_state);
		chatter_yaw.publish(yaw_state);
		chatter_x_vel.publish(x_vel_state);
		chatter_y_vel.publish(y_vel_state);
		chatter_ang_vel.publish(ang_vel_state);
		loop_rate.sleep();
	}
}

