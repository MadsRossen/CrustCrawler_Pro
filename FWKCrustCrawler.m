%Forward Kinematics for Crustcrawler
%Gruppe 364 3. semester 2020
%NOTE: That to use this code, the files TDH.m and eulerZYX.m should be in 
%the same workfolder as the file with this code
%% Input theta values
syms PI theta1 theta2 theta3 theta4 theta5 theta6 theta7
%Input theta values instead of t1, t2, t3

theta1=deg2rad(30);
theta2=deg2rad(40);
theta3=deg2rad(50);

%% Links measurements for Crust crawler in Fusion
%Fusion craustcrawler
l_1_fusion = 173.6000;   %mm
l_2_fusion = 219.80;   %mm
l_3_fusion = 277.8;    %mm
%Real life crustcrawler
l_1 = 194;   %mm
l_2 = 219;   %mm
l_3 = 261;    %mm

%% DH parametre
T01=TDH(0       ,   0    ,    0     ,theta1     );
T12=TDH(0       ,   0    , pi/2     ,theta2+pi/2);
T23=TDH(l_2     ,   0    ,    0     ,theta3     );

%% The T03 matrix
T03 = T01*T12*T23;

%% The TBW matrix
% Transalation from baseframe to frame0 along Z-axis with l_1
TB0_tz= [   1       0       0       0;
            0       1       0       0;
            0       0       1     l_1;
            0       0       0       1];

TB0 = TB0_tz;

% Transalation of frame3 to EE along X-axis with l_4

T3EE_tx= [  1       0       0     l_3;
            0       1       0       0;
            0       0       1       0;
            0       0       0       1];
        
T3EE = T3EE_tx;

TBEE = TB0*T03*T3EE;

%% Transform TBW the to Euler angles and ZYX positions and euler angles
R = eulerZYX(TBEE);
RPY = R*180/pi;
POS=[TBEE(1,4),TBEE(2,4),TBEE(3,4),R*180/pi]


x_pos = TBEE(1,4)
y_pos = TBEE(2,4)
z_pos = TBEE(3,4)
Roll_x  = RPY(1)
Pitch_y = RPY(2)
Yaw_z   = RPY(3)
