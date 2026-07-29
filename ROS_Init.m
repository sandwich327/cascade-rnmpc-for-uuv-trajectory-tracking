% 这里的ip换成自己的ip
% ip='192.168.3.87';
ip = '192.168.126.129';
setenv('ROS_MASTER_URI','http://192.168.126.129:11311');
% setenv('ROS_IP','192.168.223.2');
rosinit(ip);
%机器人实际力输出看 /rexrov/thrusters/0/thrust_wrench
