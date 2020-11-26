%Dynamics for crustcrawler
%Gruppe 364 3. semester 2020
%% Setting up symbolic values (for symbolic calculation)
syms g I_xx_1 I_xy_1 I_xz_1 I_xy_1 I_yy_1 I_yz_1 I_xz_1 I_yz_1 I_zz_1
syms I_xx_2 I_xy_2 I_xz_2 I_xy_2 I_yy_2 I_yz_2 I_xz_2 I_yz_2 I_zz_2
syms I_xx_3 I_xy_3 I_xz_3 I_xy_3 I_yy_3 I_yz_3 I_xz_3 I_yz_3 I_zz_3
syms l1 lc1 l2 lc2 l3 lc3
syms m1 m2 m3
syms th1 th2 th3
syms thd1 thd2 thd3
syms thdd1 thdd2 thdd3
syms thd12 thd13 thd23
syms thd1sq thd2sq thd3sq

%% Gravity vector

g_vec = [
    0
    0
    g];

%% Inertia tensors
IC_1 = [
    I_xx_1  I_xy_1  I_xz_1
    I_xy_1  I_yy_1  I_yz_1
    I_xz_1  I_yz_1  I_zz_1];

IC_2 = [
    I_xx_2  I_xy_2  I_xz_2
    I_xy_2  I_yy_2  I_yz_2
    I_xz_2  I_yz_2  I_zz_2];

IC_3 = [
    I_xx_3  I_xy_3  I_xz_3
    I_xy_3  I_yy_3  I_yz_3
    I_xz_3  I_yz_3  I_zz_3];

%% Ratation matrices;
R_01 = [
    cos(th1)  -sin(th1)  0
    sin(th1)  cos(th1)   0
    0          0         1];

R_12 = [
    -sin(th2)  -cos(th2)  0
    0          0        -1
    cos(th2)   -sin(th2)  0];

R_23 = [
    cos(th3)   -sin(th3)  0
    sin(th3)   cos(th3)   0
    0              0      1];

R_02 = R_01 * R_12;

R_03 = R_01 * R_12 * R_23;

%% Unit vectors
e1 = [
    1
    0
    0];

e3 = [
    0
    0
    1];

%% Positions vectors
Pv1 = R_01 * (l1 * e3);
Pv2 = R_02 * (l2 * e1);

PvC1 = R_01 * lc1 * e3;
PvC2 = R_02 * lc2 * e1;
PvC3 = R_03 * lc3 * e1;

PC1 = PvC1;
PC2 = PvC2 + Pv1;
PC3 = PvC3 + Pv1 + Pv2;

%% Angular velocities
Om11 = thd1 * e3;
Om22 = inv(R_12) * Om11 + thd2 * e3;
Om33 = inv(R_23) * Om22 + thd3 * e3;

Om01 = R_01 * thd1 * e3;
Om02 = R_01 * Om11 + R_02 * thd2 * e3;
Om03 = R_02 * Om22 + R_03 * thd3 * e3;

%% Linear velocities
Vc1 = cross(Om01,PC1);
Vc2 = cross(Om02,PC2);
Vc3 = cross(Om03,PC3);

%% Kinetic energies
T1 = 1/2 * m1 * transpose(Vc1) * Vc1 + 1/2 * transpose(Om11) * IC_1 ...
    * Om11;
T2 = 1/2 * m2 * transpose(Vc2) * Vc2 + 1/2 * transpose(Om22) * IC_2 ...
    * Om22;
T3 = 1/2 * m3 * transpose(Vc3) * Vc3 + 1/2 * transpose(Om33) * IC_3 ...
    * Om33;

%% potential energies
V1 = -m1 * transpose(g_vec) * PC1;
V2 = -m2 * transpose(g_vec) * PC2;
V3 = -m3 * transpose(g_vec) * PC3;

%% The lagrangian
T_total = simplify(T1 + T2 + T3);
V_total = simplify(V1 + V2 + V3);

Langragian_system = T_total-V_total;

%% Torques
%tau 1
diff_L_th1    = diff(Langragian_system,th1);
diff_L_thd1   = diff(Langragian_system,thd1);
diff_time_diff_L_th1_d...
    = diff(diff_L_thd1,th1) * thd1      ...
    + diff(diff_L_thd1,th2) * thd2      ...
    + diff(diff_L_thd1,th3) * thd3      ...
    + diff(diff_L_thd1,thd1) * thdd1    ...
    + diff(diff_L_thd1,thd2) * thdd2    ...
    + diff(diff_L_thd1,thd3) * thdd3;

tau_1 = simplify(diff_time_diff_L_th1_d - diff_L_th1)

%tau 2
diff_L_th2      = diff(Langragian_system,th2);
diff_L_thd2   = diff(Langragian_system,thd2);
diff_time_diff_L_th2_d...
    = diff(diff_L_thd2,th1) * thd1      ...
    + diff(diff_L_thd2,th2) * thd2      ...
    + diff(diff_L_thd2,th3) * thd3      ...
    + diff(diff_L_thd2,thd1) * thdd1    ...
    + diff(diff_L_thd2,thd2) * thdd2    ...
    + diff(diff_L_thd2,thd3) * thdd3;

tau_2 = simplify(diff_time_diff_L_th2_d - diff_L_th2)

%tau 3
diff_L_th3      = diff(Langragian_system,th3);
diff_L_thd3   = diff(Langragian_system,thd3);
diff_time_diff_L_th3_d...
    = diff(diff_L_thd3,th1) * thd1      ...
    + diff(diff_L_thd3,th2) * thd2      ...
    + diff(diff_L_thd3,th3) * thd3      ...
    + diff(diff_L_thd3,thd1) * thdd1    ...
    + diff(diff_L_thd3,thd2) * thdd2    ...
    + diff(diff_L_thd3,thd3) * thdd3;

tau_3 = simplify(diff_time_diff_L_th3_d - diff_L_th3)

%% Statespace equation terms

% Mass Matrix
MassMatrix = [
    simplify(diff(tau_1,thdd1)),...
    simplify(diff(tau_1,thdd2)),...
    simplify(diff(tau_1,thdd3));...
    ...
    simplify(diff(tau_2,thdd1)),...
    simplify(diff(tau_2,thdd2)),...
    simplify(diff(tau_2,thdd3));...
    ...
    simplify(diff(tau_3,thdd1)),...
    simplify(diff(tau_3,thdd2)),...
    simplify(diff(tau_3,thdd3))
    ]

% subs function is used to substitute ex thd1*thd2 to thd12, for partial
% differentiation
t1_thdthd_sub=...
    subs(tau_1,[thd1*thd2,thd1*thd2,thd2*thd3],[thd12,thd13,thd23]);
t2_thdthd_sub=...
    subs(tau_2,[thd1*thd2,thd1*thd2,thd2*thd3],[thd12,thd13,thd23]);
t3_thdthd_sub=...
    subs(tau_3,[thd1*thd2,thd1*thd2,thd2*thd3],[thd12,thd13,thd23]);

% Coriolis coeficients
CoriolisMatrix = [
    simplify(diff(t1_thdthd_sub,thd12)),...
    simplify(diff(t1_thdthd_sub,thd13)),...
    simplify(diff(t1_thdthd_sub,thd23));...
    ...
    simplify(diff(t2_thdthd_sub,thd12)),...
    simplify(diff(t2_thdthd_sub,thd13)),...
    simplify(diff(t2_thdthd_sub,thd23));...
    ...
    simplify(diff(t3_thdthd_sub,thd12)),...
    simplify(diff(t3_thdthd_sub,thd13)),...
    simplify(diff(t3_thdthd_sub,thd23))
    ]

% subs function is used to substitute ex thd1^2 to thd1sq, for partial
% differentiation
t1_thdsq_sub=subs(tau_1,[thd1^2,thd2^2,thd3^2],[thd1sq,thd2sq,thd3sq]);
t2_thdsq_sub=subs(tau_2,[thd1^2,thd2^2,thd3^2],[thd1sq,thd2sq,thd3sq]);
t3_thdsq_sub=subs(tau_3,[thd1^2,thd2^2,thd3^2],[thd1sq,thd2sq,thd3sq]);

% centrifugal coeficients
CentrifugalMatrix = [
    simplify(diff(t1_thdsq_sub,thd1sq)),...
    simplify(diff(t1_thdsq_sub,thd2sq)),...
    simplify(diff(t1_thdsq_sub,thd3sq));...
    ...
    simplify(diff(t2_thdsq_sub,thd1sq)),...
    simplify(diff(t2_thdsq_sub,thd2sq)),...
    simplify(diff(t2_thdsq_sub,thd3sq));...
    ...
    simplify(diff(t3_thdsq_sub,thd1sq))...
    simplify(diff(t3_thdsq_sub,thd2sq))...
    simplify(diff(t3_thdsq_sub,thd3sq))
    ]

thdthdVector = [thd1*thd2;thd1*thd3;thd2*thd3];
thdsqVector = [thd1^2;thd2^2;thd3^2];

% VelocityVector
VelocityVector = CoriolisMatrix*thdthdVector+CentrifugalMatrix*thdsqVector

% GravityVector
GravityVector_no_g  = [
    simplify(diff(tau_1,g));
    simplify(diff(tau_2,g));
    simplify(diff(tau_3,g))
    ];
GravityVector = g*GravityVector_no_g

%% Statespace equation
jointAngleAccVector = [
    thdd1;
    thdd2;
    thdd3
    ];

%Inverse Dynamics Equation
tauVector = simplify(MassMatrix*jointAngleAccVector+VelocityVector...
    +GravityVector)

%Forward Dynamics Equation
jointAngleAccVector = simplify(inv(MassMatrix)*(tauVector-...
    (VelocityVector)+GravityVector))




