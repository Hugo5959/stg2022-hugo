
rm ~/workspaceRos/src/larm1_slam/maps/map.bag
# rostopic echo -p /map > ~/workspaceRos/src/larm1_slam/bags/map.bag

# rosbag record --duration=10 --output-name=~/workspaceRos/src/larm1_slam/bags/map.map /clock

# cd ~/workspaceRos/src/larm1_slam/bags
# rosbag record -O map /base_scan /tf

cd ~/workspaceRos/src/larm1_slam/maps
rosrun map_server map_saver -f map