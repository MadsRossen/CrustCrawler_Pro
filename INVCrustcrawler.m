%Inverse kinematics for Crust crawler
%Gruppe 364 3. semester 2020
%NOTE: That to use this code, the files TDH.m and eulerZYX.m should be in 
%the same workfolder as the file with this code
%% Setup
clear
syms theta1 theta2 theta3 theta4 theta5 theta6 

%% Links measurements for Crust crawler in Fusion
%Fusion craustcrawler
l_1_fusion = 173.6000;   %mm
l_2_fusion = 219.80;   %mm
l_3_fusion = 277.8;    %mm
%Real life crustcrawler
l_1 = 194;   %mm
l_2 = 219;   %mm
l_3 = 261;    %mm

%% Setup of data
%Enter end effector coordinates and orientation.
   X        =  -362.938;
   Y        =  -209.542;
   Z        =  341.977;
   
   Roll     =  -180;
   Pitch    =  30;
   Yaw      =  90;

% Matrix setup.
TBEEtarget = eulerZYX2T(X,Y,Z,Roll*pi/180,Pitch*pi/180,Yaw*pi/180);

% Transalation from baseframe to frame0 along Z-axis with l_1
TB0= [  1       0       0       0;
        0       1       0       0;
        0       0       1     l_1;
        0       0       0       1];

%% Compute theta1:
theta1_1    = pi+atan2(TBEEtarget(2,4),TBEEtarget(1,4));
theta1_2    = atan2(TBEEtarget(2,4),TBEEtarget(1,4));

%% Compute theta2:
T01=TDH(0       ,   0    ,    0     ,theta1_1     );
T1EEtarget  = inv(TB0)*inv(T01)*TBEEtarget;
l_4         = sqrt((-T1EEtarget(1,4))^2+(T1EEtarget(3,4))^2);
phi_1       = acos((l_4^2+l_2^2-l_3^2)/(2*l_4*l_2));
phi_2       = atan2(T1EEtarget(3,4),-T1EEtarget(1,4));
theta2_1    = pi/2-phi_1-phi_2;
theta2_2    = -theta2_1;

%% Compute theta3:
phi_3       = acos((l_2^2+l_3^2-l_4^2)/(2*l_2*l_3));
theta3_1    = pi-phi_3;
theta3_2    = -theta3_1;

%% Solutions
% Solution1
Theta1_1_deg   = round(theta1_1*180/pi,5)
%Theta1_2_deg   = round(theta1_2*180/pi,5)
theta2_1_deg   = round(theta2_1*180/pi,5)
%theta2_2_deg   = round(theta2_2*180/pi,5)
theta3_1_deg   = round(theta3_1*180/pi,5)
%theta3_2_deg   = round(theta3_2*180/pi,5)








