# **USV MPC (Model Predictive Control ** on VRX simulator


## Installation

On an ubuntu 20.04 + Matlab (version R2022a) (you need a license to use it) + ROS Noetic (version ros-noetic-desktop-full) :

``` bash
cd $HOME
mkdir -p ~/workspaceRos_VRX
cd ~/workspaceRos_VRX
git clone https://bitbucket.org/imt-mobisyst/stg2022-hugo.git
cd $HOME/workspaceRos_VRX
catkin_make

```
Some errors messages will probably ask you to install additional ros-noetic packages, please install them


## Launch VRX simulator

```bash
cd $HOME/workspaceRos_VRX/stg2022-hugo
roslaunch algo_mpc robot_stage_vrx.launch
```

## Launch Matlab MPC Node

Open your Matlab Application and then launch the file nodeMPC.m locate in $HOME/workspaceRos_VRX/stg2022-hugo/matlab_node_mpc

`
## Contributor

- Hugo Reubrecht

