# **USV MPC (Model Predictive Control) ** on VRX simulator


## Installation


On an ubuntu 20.04 + Matlab (version R2022a) (you need a license to use it) + ROS Noetic (version ros-noetic-desktop-full) :

``` bash
sudo apt update
sudo apt full-upgrade
```


``` bash
sudo apt install ros-noetic-geographic-msgs
sudo apt install ros-noetic-geodesy
sudo apt install -y build-essential cmake cppcheck curl git gnupg libeigen3-dev libgles2-mesa-dev lsb-release pkg-config protobuf-compiler qtbase5-dev python3-dbg python3-pip python3-venv ruby software-properties-common wget 
sudo sh -c 'echo "deb http://packages.ros.org/ros/ubuntu $(lsb_release -sc) main" > /etc/apt/sources.list.d/ros-latest.list'
sudo apt-key adv --keyserver 'hkp://keyserver.ubuntu.com:80' --recv-key C1CF6E31E6BADE8868B172B4F42ED6FBAB17C654
sudo sh -c 'echo "deb http://packages.osrfoundation.org/gazebo/ubuntu-stable `lsb_release -cs` main" > /etc/apt/sources.list.d/gazebo-stable.list'
wget http://packages.osrfoundation.org/gazebo.key -O - | sudo apt-key add -
sudo apt update
DIST=noetic
GAZ=gazebo11
sudo apt install ${GAZ} lib${GAZ}-dev ros-${DIST}-gazebo-plugins ros-${DIST}-gazebo-ros ros-${DIST}-hector-gazebo-plugins ros-${DIST}-joy ros-${DIST}-joy-teleop ros-${DIST}-key-teleop ros-${DIST}-robot-localization ros-${DIST}-robot-state-publisher ros-${DIST}-joint-state-publisher ros-${DIST}-rviz ros-${DIST}-ros-base ros-${DIST}-teleop-tools ros-${DIST}-teleop-twist-keyboard ros-${DIST}-velodyne-simulator ros-${DIST}-xacro ros-${DIST}-rqt ros-${DIST}-rqt-common-plugins
```


``` bash
cd $HOME
mkdir -p ~/workspaceRos_VRX
cd ~/workspaceRos_VRX
git clone https://github.com/Hugo5959/stg2022-hugo.git
cd $HOME/workspaceRos_VRX/stg2022-hugo
catkin_make
source devel/setup.bash
```
Some error messages will probably ask you to install additional ros-noetic packages, please install them.


## Launch VRX simulator

```bash
cd $HOME/workspaceRos_VRX/stg2022-hugo
roslaunch algo_mpc robot_stage_vrx.launch
```

## Launch Matlab MPC Node

- Open your Matlab Application and then open the file nodeMPC.m located in $HOME/workspaceRos_VRX/stg2022-hugo/matlab_node_mpc.
(new : you can open the file nodeMPC_lac_heron.m instead for the simulation of the new MPC on the Heron Lake with the oxygen level)

- You will have to adapt the path for the casadi library (in "Importations" section) (normally located at $HOME/workspaceRos_VRX/stg2022-hugo/matlab_node_mpc/casadi-linux-matlabR2014b-v3.5.5)  and the IP adress of your node Ros master (in "Ros Architecture" section) (you can check it with a "ifconfig" command) (execute "sudo apt install net-tools" in command line if the command is not recognized by your system).

- Also add the tbxmanager library in your Matlab path (normally located at $HOME/workspaceRos_VRX/stg2022-hugo/matlab_node_mpc/tbxmanager), it is used for the build of Polyhedron target set).

## Multi-robots (new)

(Please try with one robot before to better understand the simulation with two robots)

- To run a simulation with two robots on VRX:

```bash
cd $HOME/workspaceRos_VRX/stg2022-hugo
roslaunch algo_mpc robot_stage_vrx.launch
```
and in an other shell : 

```bash
cd $HOME/workspaceRos_VRX/stg2022-hugo
roslaunch algo_mpc robot_stage_vrx_1.launch
```

- Then you run in two different Matlab sessions : "nodeMPC_lac_heron.m" and "nodeMPC_lac_heron_robot2.m"

## External links

For more informations about VRX : https://github.com/osrf/vrx

about casadi installation : https://web.casadi.org/get/

about tbxmanager installation : https://www.tbxmanager.com

## Contributors

- Hugo Reubrecht
- Alejandro Anderson (MPC algorithm)
- Eric Duviella (Map generation)
